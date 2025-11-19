! (C) British Crown Copyright 2025 Met Office
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module for reconstructed radiance observation operator

module ufo_reconradop_mod

  use fckit_log_module, only : fckit_log
  use iso_c_binding, only: c_int, c_double
  use kinds, only : kind_real ! from oops
  use rttov_types
  use ufo_radiancerttov_utils_mod
  use ufo_rttovonedvarcheck_ob_mod
  use ufo_utils_mod, only : ufo_utils_iogetfreeunit

  implicit none
  private

  integer, parameter                    :: max_string=800

  !> Fortran derived type for the observation type
  type, public :: ufo_reconradop
    real(kind_real), allocatable     :: Cmatrix_bias(:)         ! C matrix bias vector used to simulate reconstructed radiances.
    real(kind_real), allocatable     :: Cmatrix(:,:)            ! C matrix used to simulate reconstructed radiances.
    real(kind_real), allocatable     :: cwn(:)
    real(kind_real), allocatable     :: rad(:)
    real(kind_real), allocatable     :: bt_arr(:)
    real(kind_real), allocatable     :: dBT_dRad_arr(:)
    real(kind_real), allocatable     :: profile_k_rr(:,:,:)
    real(kind_real), allocatable     :: profile_k_rr_surf(:,:)
  contains
    procedure :: GetCmatrix => ufo_reconradop_GetCmatrix
    procedure :: setup  => ufo_reconradop_setup
    procedure :: allocate => ufo_reconradop_allocate
    procedure :: hofx_jac_calc => ufo_reconradop_hofx_jac_calc
    procedure :: delete => ufo_reconradop_delete
  end type ufo_reconradop

contains

  !! Read the C matrix that multiplies the obs operator for reconstructed radiances from file
  !! Based on ufo_rttovonedvarcheck_GetEmisEigenVec
  subroutine ufo_reconradop_GetCmatrix(self, filepath)

    implicit none

    ! Subroutine arguments:
    class(ufo_reconradop), target, intent(inout)     :: self
    character(len=*), intent(in)                     :: filepath

    ! Local declarations:
    character(len=*), parameter :: RoutineName = "ufo_reconradop_GetCmatrix"
    integer :: nchans_first, nchans_second, nchans, ich
    integer :: readstatus
    integer :: fileunit
    character(len=max_string)   :: message

    ! Open file for reading
    fileunit = ufo_utils_iogetfreeunit()
    open(unit = fileunit, file = trim(filepath))

    read(fileunit, *, iostat = readstatus) nchans

    allocate(self % Cmatrix_bias(nchans))
    allocate(self % Cmatrix(nchans,nchans))

    read (fileunit, *, iostat = readstatus) self % Cmatrix_bias(:)
    do ich = 1, nchans
      read (fileunit, *, iostat = readstatus) self % Cmatrix(:,ich)
    end do

    if (readstatus /= 0) then
      write(message,*) RoutineName,  &
      ' Problem reading C matrix elements '
      call abor1_ftn(message)
    end if

    write(*, '(A,I0,A)') 'Finished reading the C matrix and bias vector for ', nchans, ' reconstructed radiances.'

    close(unit = fileunit)

  end subroutine ufo_reconradop_GetCmatrix

  ! ------------------------------------------------------------------------------
  subroutine ufo_reconradop_setup(self, rtprof, conf)

    implicit none
    class(ufo_reconradop), intent(inout) :: self
    type(ufo_rttov_io),    intent(inout) :: rtprof
    type(rttov_conf),      intent(inout) :: conf
    character(len=max_string)            :: message

    call self % getCmatrix(conf % Cmatrix_path)
    write(*,*) "Cmatrix and bias read in ufo_radiancerttov_setup"
    if (size(self % Cmatrix_bias) /=  rtprof % nchan_inst &
        .or. size(self % Cmatrix,1) /= rtprof % nchan_inst &
        .or. size(self % Cmatrix,2) /= rtprof % nchan_inst) then
        write(*,*) "C matrix or bias vector size does not match expected size"
        write(*,*) "C matrix bias size: ", size(self % Cmatrix_bias)
        write(*,*) "C matrix size: ", size(self % Cmatrix,1), "x", size(self % Cmatrix,2)
        write(*,*) "Channels size: ", rtprof % nchan_inst
        call abor1_ftn("C matrix size does not match channel size")
    end if
    if (conf % RTTOV_switchrad) then
      ! the C matrix needs to be applied to jacobians in radiance units
      write(*,*) 'RTTOV jacobians now set to be given in radiance units'
      conf % rttov_opts % rt_all % switchrad = .false.
    end if

    message = 'Finished reading C matrix and bias in ufo_reconradop_setup'
    call fckit_log%info(message)

  end subroutine ufo_reconradop_setup
  ! ------------------------------------------------------------------------------

  subroutine ufo_reconradop_allocate(self, rtprof)

    implicit none
    class(ufo_reconradop), intent(inout) :: self
    type(ufo_rttov_io),    intent(inout) :: rtprof

    integer                              :: nlevels      !Number of levels in the profile 

    nlevels = size(rtprof % profiles(1) % p)
    allocate(self % cwn(rtprof % nchan_inst))
    allocate(self % rad(rtprof % nchan_inst))
    allocate(self % bt_arr(rtprof % nchan_inst))
    allocate(self % dBT_dRad_arr(rtprof % nchan_inst))
    allocate(self % profile_k_rr(rtprof % nchan_inst, nlevels,2))
    allocate(self % profile_k_rr_surf(rtprof % nchan_inst,9))

  end subroutine ufo_reconradop_allocate
  ! ------------------------------------------------------------------------------

  subroutine ufo_reconradop_hofx_jac_calc(self, rtprof, conf, chanprof, prof_start, jacobian_needed, hofx, ob_info)

    implicit none
    class(ufo_reconradop),   intent(inout) :: self
    type(ufo_rttov_io),      intent(inout) :: rtprof
    type(rttov_conf),        intent(in)    :: conf
    type(rttov_chanprof),    intent(in)    :: chanprof(:)
    integer,                 intent(in)    :: prof_start
    logical,                 intent(in)    :: jacobian_needed
    real(c_double),optional, intent(inout) :: hofx(:,:)
    type(ufo_rttovonedvarcheck_ob), optional, intent(inout) :: ob_info 
    
    character(len=max_string)                   :: message

    integer                              :: ichan, nchan_sim, iprof, localchan, nlevels

    nlevels = size(rtprof % profiles(1) % p)
    nchan_sim = size(chanprof)
    do ichan = 1, nchan_sim, rtprof % nchan_inst
      iprof = prof_start + chanprof(ichan) % prof - 1          
      self % cwn = rtprof % radiance % total(ichan:ichan+rtprof % nchan_inst-1) ! temporary array, in mW/cm-1/sr/m2
      ! Cmat_bias is in W/m-1/sr/m2
      self % rad(1:rtprof % nchan_inst) = matmul(self % Cmatrix, 1e-5*self % cwn) + self % Cmatrix_bias(1:rtprof % nchan_inst) ! simulated reconstructed radiance from mW/cm-1/sr/m2 to W/m-1/sr/m2
      ! make sure radiance is not negative
      if (any(self % rad(1:rtprof % nchan_inst) <= 0.0)) then
        write(*,*) 'Warning: radiance <= 0:', pack(self % rad(1:rtprof % nchan_inst), self % rad(1:rtprof % nchan_inst) <= 0.0), &
          ' chan:', pack(chanprof(ichan:ichan+rtprof % nchan_inst-1) % chan, self % rad(1:rtprof % nchan_inst) <= 0.0), ' iprof:', iprof
        where (self % rad(1:rtprof % nchan_inst) <= 0.0)
          ! set negative simulated reconstructed radiances equal to simulated radiances
          self % rad(1:rtprof % nchan_inst) = 1e-5*self % cwn(1:rtprof % nchan_inst) ! radiance from mW/cm-1/sr/m2 to W/m-1/sr/m2
        end where
        ! convert wavenumber from cm-1 to m-1
        ! note we are assuming we are not using cut-down coefficients so that chanprof(ichan)%chan is the actual channel number
        self % cwn = conf % rttov_coef_array(1) % coef % ff_cwn(chanprof(ichan:ichan+rtprof % nchan_inst-1) % chan) * 1e2
        ! convert radiance to brightness temperature 
        call rtprof % rad_to_bt(self % cwn(1:rtprof % nchan_inst), self % rad(1:rtprof % nchan_inst), self % bt_arr(1:rtprof % nchan_inst))              
        if (present(hofx)) then
          hofx(1:rtprof % nchan_inst,iprof) = self % bt_arr(1:rtprof % nchan_inst)
          ! this is optional, to use BT directly from RTTOV when rad is not positive
          ! (perhaps more precise when the conversion depends of the channel SRF?)
          where (self % rad(1:rtprof % nchan_inst) <= 0.0)
            ! set bt_arr equal to simulated BT from RTTOV
            hofx(1:rtprof % nchan_inst,iprof) = rtprof % radiance % bt(ichan:ichan+rtprof % nchan_inst-1)
          end where
        end if
      else
        ! convert wavenumber from cm-1 to m-1
        ! note we are assuming we are not using cut-down coefficients so that chanprof(ichan)%chan is the actual channel number
        self % cwn = conf % rttov_coef_array(1) % coef % ff_cwn(chanprof(ichan:ichan+rtprof % nchan_inst-1) % chan) * 1e2
        ! convert radiance to brightness temperature 
        call rtprof % rad_to_bt(self % cwn(1:rtprof % nchan_inst), self % rad(1:rtprof % nchan_inst), self % bt_arr(1:rtprof % nchan_inst))              
        if (present(hofx)) &
          hofx(1:rtprof % nchan_inst,iprof) = self % bt_arr(1:rtprof % nchan_inst)  
      end if              
      ! now redefines the RTprof quantities that are going to end up in hofxdiags
      ! check ufo_rttov_populate_hofxdiags to make sure everything is covered
      if (present(hofx)) then
        rtprof % radiance % bt(ichan:ichan+rtprof % nchan_inst-1) = hofx(1:rtprof % nchan_inst,iprof)
      else
        ! redefine the BT when radiance is positive
        ! and use BT unchanged from RTTOV when rad is not positive 
        ! (perhaps more precise when the conversion depends of the channel SRF?)
        where (self % rad(1:rtprof % nchan_inst) > 0.0)
          rtprof % radiance % bt(ichan:ichan+rtprof % nchan_inst-1) = self % bt_arr(1:rtprof % nchan_inst)
        end where
      end if
    
      !store transmittance if ob_info present in call and transmittance part of structure
      if (present(ob_info)) then
        if (allocated(ob_info % transmittance)) then
          if (conf % do_mw_scatt) then
            ob_info % transmittance(1:rtprof % nchan_inst) = &
              rtprof % mw_scatt % emis_retrieval % tau_clr(ichan:ichan+rtprof % nchan_inst-1)
          else
            ob_info % transmittance(1:rtprof % nchan_inst) = &
              rtprof % transmission % tau_total(ichan:ichan+rtprof % nchan_inst-1)
          end if
        end if
      end if

      if (jacobian_needed) then
        ! Now apply the C matrix to the Jacobian
        do localchan = 1, rtprof % nchan_inst
          ! var_ts
          self % profile_k_rr(localchan,:,1) = rtprof % profiles_k(ichan + localchan - 1) % t(:)
          ! var_q
          self % profile_k_rr(localchan,:,2) = rtprof % profiles_k(ichan + localchan - 1) % q(:)
          ! var_sfc_t2m
          self % profile_k_rr_surf(localchan,1) = rtprof % profiles_k(ichan + localchan - 1) % s2m % t
          ! var_sfc_q2m
          self % profile_k_rr_surf(localchan,2) = rtprof % profiles_k(ichan + localchan - 1) % s2m % q
          ! var_sfc_u10
          self % profile_k_rr_surf(localchan,3) = rtprof % profiles_k(ichan + localchan - 1) % s2m % u
          ! var_sfc_u10
          self % profile_k_rr_surf(localchan,4) = rtprof % profiles_k(ichan + localchan - 1) % s2m % v
          ! surf press
          self % profile_k_rr_surf(localchan,5) = rtprof % profiles_k(ichan + localchan - 1) % s2m % p
          ! var_sfc_tskin
          self % profile_k_rr_surf(localchan,6) = rtprof % profiles_k(ichan + localchan - 1) % skin % t
          ! var_sfc_emiss
          self % profile_k_rr_surf(localchan,7) = rtprof % emissivity_k(ichan + localchan - 1) % emis_in
          ! ctp
          self % profile_k_rr_surf(localchan,8) = rtprof % profiles_k(ichan + localchan - 1) % ctp
          ! cfraction
          self % profile_k_rr_surf(localchan,9) = rtprof % profiles_k(ichan + localchan - 1) % cfraction
        end do
        self % profile_k_rr(:,:,1) = matmul(self % Cmatrix,self % profile_k_rr(:,:,1))
        self % profile_k_rr(:,:,2) = matmul(self % Cmatrix,self % profile_k_rr(:,:,2))
        self % profile_k_rr_surf(:,1:9) = matmul(self % Cmatrix,self % profile_k_rr_surf(:,1:9))
        ! Calculate the Jacobian of BT (K) with respect to radiance (W/m-1/sr/m2)
        call rtprof % dBT_dRad(self % cwn, self % rad, self % dBT_dRad_arr)
        ! Convert the jacobians back to BT over field variable units 
        do localchan = 1, rtprof % nchan_inst
          ! var_ts
          rtprof % profiles_k(ichan + localchan - 1) % t(:) = &
            self % dBT_dRad_arr(localchan)*1e-5*self % profile_k_rr(localchan,:,1)
          ! var_q      
          rtprof % profiles_k(ichan + localchan - 1) % q(:) = &
            self % dBT_dRad_arr(localchan)*1e-5*self % profile_k_rr(localchan,:,2)
          ! var_sfc_t2m
          rtprof % profiles_k(ichan + localchan - 1) % s2m % t = &
            self % dBT_dRad_arr(localchan)*1e-5*self % profile_k_rr_surf(localchan,1)
          ! var_sfc_q2m
          rtprof % profiles_k(ichan + localchan - 1) % s2m % q = &
            self % dBT_dRad_arr(localchan)*1e-5*self % profile_k_rr_surf(localchan,2)
          ! var_sfc_u10
          rtprof % profiles_k(ichan + localchan - 1) % s2m % u = &
            self % dBT_dRad_arr(localchan)*1e-5*self % profile_k_rr_surf(localchan,3)
          ! var_sfc_u10
          rtprof % profiles_k(ichan + localchan - 1) % s2m % v = &
            self % dBT_dRad_arr(localchan)*1e-5*self % profile_k_rr_surf(localchan,4)
          ! surf press
          rtprof % profiles_k(ichan + localchan - 1) % s2m % p = &
            self % dBT_dRad_arr(localchan)*1e-5*self % profile_k_rr_surf(localchan,5)
          ! var_sfc_tskin
          rtprof % profiles_k(ichan + localchan - 1) % skin % t = &
            self % dBT_dRad_arr(localchan)*1e-5*self % profile_k_rr_surf(localchan,6)
          ! var_sfc_emiss
          rtprof % emissivity_k(ichan + localchan - 1) % emis_in = &
            self % dBT_dRad_arr(localchan)*1e-5*self % profile_k_rr_surf(localchan,7)
          ! ctp
          rtprof % profiles_k(ichan + localchan - 1) % ctp = &
            self % dBT_dRad_arr(localchan)*1e-5*self % profile_k_rr_surf(localchan,8)
          ! cfraction
          rtprof % profiles_k(ichan + localchan - 1) % cfraction = &
            self % dBT_dRad_arr(localchan)*1e-5*self % profile_k_rr_surf(localchan,9)
        end do
      end if 
    enddo

  end subroutine ufo_reconradop_hofx_jac_calc
  ! ------------------------------------------------------------------------------
  
  subroutine ufo_reconradop_delete(self)

    implicit none
    class(ufo_reconradop), intent(inout) :: self

    if (allocated(self % cwn)) deallocate(self % cwn)
    if (allocated(self % rad)) deallocate(self % rad)
    if (allocated(self % bt_arr)) deallocate (self % bt_arr)
    if (allocated(self % dBT_dRad_arr)) deallocate (self % dBT_dRad_arr)
    if (allocated(self % profile_k_rr)) deallocate (self % profile_k_rr)
    if (allocated(self % profile_k_rr_surf)) deallocate (self % profile_k_rr_surf)

  end subroutine ufo_reconradop_delete

end module ufo_reconradop_mod
