! (C) Copyright 2025 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module for CRTM reconstructed radiance observation operator

module ufo_reconradop_crtm_mod

 use crtm_module
 use CRTM_Planck_Functions
 use fckit_log_module, only : fckit_log
 use, intrinsic :: iso_c_binding, only: c_int, c_double
 use kinds, only : kind_real ! from oops
 use ufo_crtm_utils_mod
 use ufo_utils_mod, only : ufo_utils_iogetfreeunit

 implicit none
 private

 !> Fortran derived type for cmatrix
 type, public :: ufo_reconradop_crtm
 contains
   procedure :: GetCmatrix => ufo_reconradop_crtm_GetCmatrix
   procedure :: setup  => ufo_reconradop_crtm_setup
   procedure :: apply => ufo_reconradop_crtm_apply
   procedure :: delete => ufo_reconradop_crtm_delete
 end type ufo_reconradop_crtm

contains


! ------------------------------------------------------------------------------

! same file format as rttov cmatrix
subroutine ufo_reconradop_crtm_GetCmatrix(self, filepath, Cmatrix_bias, Cmatrix)

implicit none

! Subroutine arguments:
class(ufo_reconradop_crtm),   intent(in)       :: self
character(len=*), intent(in)                   :: filepath
real(kind_real), intent(inout),allocatable     :: Cmatrix_bias(:)         ! C matrix bias vector used to simulate reconstructed radiances.
real(kind_real), intent(inout),allocatable     :: Cmatrix(:,:)            ! C matrix used to simulate reconstructed radiances.
! Local declarations:
character(len=*), parameter :: RoutineName = "ufo_reconradop_crtm_GetCmatrix"
integer :: nchans, ich
integer :: readstatus
integer :: fileunit
character(len=max_string)   :: message
character(len=4) :: char_chan
! Open file for reading
 fileunit = ufo_utils_iogetfreeunit()
 open(unit = fileunit, file = trim(filepath))

 read(fileunit, *, iostat = readstatus) nchans
 if (readstatus /= 0) then
   message = RoutineName//" Problem reading nchans."
   call abor1_ftn(message)
 end if

 allocate( Cmatrix_bias(nchans))
 allocate( Cmatrix(nchans,nchans))

 read (fileunit, *, iostat = readstatus)  Cmatrix_bias(:)
 if (readstatus /= 0) then
   message = RoutineName//" Problem reading C matrix bias."
   call abor1_ftn(message)
 end if
 do ich = 1, nchans
   read (fileunit, *, iostat = readstatus)  Cmatrix(:,ich)
   if (readstatus /= 0) then
     message = RoutineName//" Problem reading C matrix elements."
     call abor1_ftn(message)
   end if
 end do
 close(unit = fileunit)

 write(char_chan,"(I4)") nchans
 message = "Finished reading the C matrix and bias vector for "//trim(char_chan)//" reconstructed radiances."
 call fckit_log%info(message)

end subroutine ufo_reconradop_crtm_GetCmatrix

! ------------------------------------------------------------------------------

subroutine ufo_reconradop_crtm_setup(self, filepath, channels, Cmatrix_bias, Cmatrix)

implicit none
class(ufo_reconradop_crtm),   intent(in)        :: self
character(len=*), intent(in)                    :: filepath
character(len=max_string)                       :: message
integer(c_int),            intent(in)           :: channels(:)
real(kind_real), allocatable, intent(inout)     :: Cmatrix_bias(:)         ! C matrix bias vector used to simulate reconstructed radiances.
real(kind_real), allocatable, intent(inout)     :: Cmatrix(:,:)            ! C matrix used to simulate reconstructed radiances.
integer:: nchan
character(len=4) :: size_char, size_char2

 call self%GetCmatrix(filepath, Cmatrix_bias, Cmatrix)
 nchan = size(channels)

 message = "Cmatrix and bias read in ufo_radiancecrtm"
 call fckit_log%info(message)

 if (size(Cmatrix_bias) /=  nchan  .or. size(Cmatrix,1) /=  nchan .or. size(Cmatrix,2) /=  nchan) then


   message =  "C matrix or bias vector size does not match expected size"
   call fckit_log%info(message)

   write(size_char,"(I4)") size(Cmatrix_bias)
   message = "C matrix bias size: "//trim(size_char)//"."
   call fckit_log%info(message)

   write(size_char,"(I4)") size(Cmatrix,1)
   write(size_char2,"(I4)") size(Cmatrix,2)
   message = "C matrix size: "//trim(size_char)//"x"//trim(size_char2)//"."
   call fckit_log%info(message)

   write(size_char,"(I4)") nchan
   message = "Channels size: "//trim(size_char)
   call fckit_log%info(message)

   message = "C matrix size does not match channel size"
   call abor1_ftn(message)
 end if

 message = "Finished reading C matrix and bias in ufo_reconradop_crtm_setup"
 call fckit_log%info(message)

end subroutine ufo_reconradop_crtm_setup

! ------------------------------------------------------------------------------

subroutine ufo_reconradop_crtm_apply(self, filepath,  rts, n_Sensor, n_Profiles, n_Channels, n_Absorbers, n_Levels, channels, rts_K, atm_K, sfc_K, jacobian_needed)

implicit none
class(ufo_reconradop_crtm),   intent(in) :: self
character(len=*), intent(in)                     :: filepath
type(CRTM_RTSolution_type), intent(inout) :: rts(:,:)
integer, intent(in) :: n_Sensor, n_Profiles, n_Channels, n_Absorbers
integer, intent(in) :: n_Levels
integer, intent(in) :: channels(:)
type(CRTM_RTSolution_type), optional, intent(inout):: rts_K(:,:)
type(CRTM_Atmosphere_type), optional,intent(inout) :: atm_K(:,:)
type(CRTM_Surface_type), optional,intent(inout) :: sfc_K(:,:)
logical, intent(in) :: jacobian_needed
character(len=max_string)            :: message
real(c_double)  :: rad_tmp(n_Channels)
real(c_double),allocatable :: jac_t_tmp(:,:)
real(c_double),allocatable :: jac_a_tmp(:, :, :)
real(c_double),allocatable :: jac_s_tmp(:,:)
real(kind_real), allocatable     :: Cmatrix_bias(:)         ! C matrix bias vector used to simulate reconstructed radiances.
real(kind_real), allocatable     :: Cmatrix(:,:)            ! C matrix used to simulate reconstructed radiances.
integer :: m, l, jabs

 message = "Applying reconstructed radiance operator"
 call fckit_log%info(message)
 call  self%setup(filepath, channels, Cmatrix_bias, CMatrix)
 message = "jacobian_needed="//merge("T","F",jacobian_needed)
 call fckit_log%info(message)

 if(jacobian_needed) then
   allocate(jac_t_tmp(n_Channels,n_Levels),&
            jac_a_tmp(n_Channels, n_Levels, n_Absorbers ),&
            jac_s_tmp(n_Channels,8) )
 end if

 profile_loop: do m = 1, n_Profiles

   rad_tmp(:) = rts(:,m)%Radiance
   ! need factor of 1e5 for bias MTG-IRS/IASI products (W/m-1/sr/m2) to CRTM/Planck Functions (mW/cm-1/sr/m2)
   ! Apply reconstructed radiance operator matrix and add bias term
   rad_tmp = matmul(cmatrix, rad_tmp) + 1e5*cmatrix_bias

   if(jacobian_needed) then

    ! fill temporary jacobian arrays for profile
    call  ufo_reconradop_crtm_rr_jacobian_tmp_populate(m, n_Channels, n_Levels, n_Absorbers, sfc_K, rts_K, atm_K, jac_s_tmp, jac_t_tmp, jac_a_tmp)

     ! apply matrix to temperature profile Jacobian
     jac_t_tmp = matmul(cmatrix, jac_t_tmp)
     ! apply matrix to surface Jacobians
     jac_s_tmp = matmul(cmatrix, jac_s_tmp)
     ! apply matrix to constituent (H2O, O3, etc)  profile Jacobians
     do jabs = 1, n_Absorbers
       jac_a_tmp(:,:,jabs) = matmul(cmatrix, jac_a_tmp(:,:,jabs))
     end do
   end if

   ! Fill in each channel with brightness temperatures (recall we have set radiance output)
   channel_loop: do l= 1, n_Channels
   ! only change everything if RR is positive.
     if(rad_tmp(l) > 0) then

       ! replace Tb with reconstructed radiance Tb (as long as there's a valid result for RR operator)
       ! Note there isn't another Planck temperature call for when rad_tmp<0,
       ! because rts(l,m)%Brightness_Temperature will already have CRTM Tb

       call CRTM_Planck_Temperature(n_Sensor, channels(l), rad_tmp(l), rts(l,m)%Brightness_Temperature  )
       if(jacobian_needed) then
         call ufo_reconradop_crtm_rr_jacobian_convert(l, m, n_Sensor,n_Levels, n_Absorbers, channels, rad_tmp(l), jac_s_tmp, jac_t_tmp, jac_a_tmp, sfc_K, rts_K, atm_K)
       end if
       ! if we have a negative radiance, convert the jacobians that are set to be dR/dx to dTb/dx
     else if (jacobian_needed) then
       ! fill in original CRTM Tb Jacobian
       call ufo_reconradop_crtm_normal_jacobian_convert(l, m, n_Sensor, n_Levels, n_Absorbers, channels, rts(l,m)%Radiance, sfc_K, rts_K, atm_K)
     end if

   end do channel_loop
 end do profile_loop

 ! deallocate temporary arrays
 if(allocated(jac_t_tmp)) deallocate(jac_t_tmp)
 if(allocated(jac_a_tmp)) deallocate(jac_a_tmp)
 if(allocated(jac_s_tmp)) deallocate(jac_s_tmp)

 call self%delete(Cmatrix_bias, Cmatrix)

end subroutine ufo_reconradop_crtm_apply

! ------------------------------------------------------------------------------
subroutine ufo_reconradop_crtm_rr_jacobian_tmp_populate(m, n_Channels, n_Levels, n_Absorbers, sfc_K, rts_K, atm_K, jac_s_tmp, jac_t_tmp, jac_a_tmp)
implicit none
integer, intent(in) :: m, n_Channels, n_Levels, n_Absorbers
type(CRTM_RTSolution_type), intent(in):: rts_K(:,:)
type(CRTM_Atmosphere_type), intent(in) :: atm_K(:,:)
type(CRTM_Surface_type),  intent(in) :: sfc_K(:,:)
real(c_double),intent(inout) :: jac_t_tmp(:,:)
real(c_double),intent(inout) :: jac_a_tmp(:, :, :)
real(c_double),intent(inout) :: jac_s_tmp(:,:)
integer :: l, jlevel, jabs

 jac_s_tmp(:,1) = sfc_K(:,m) % Water_Temperature
 jac_s_tmp(:,2) = sfc_K(:,m) % Land_Temperature
 jac_s_tmp(:,3) = sfc_K(:,m) % Ice_Temperature
 jac_s_tmp(:,4) = sfc_K(:,m) % Snow_Temperature
 jac_s_tmp(:,5) = rts_K(:,m) % Surface_Emissivity
 jac_s_tmp(:,6) = sfc_K(:,m) % Wind_Speed
 jac_s_tmp(:,7) = sfc_K(:,m) % Wind_Direction
 jac_s_tmp(:,8) = sfc_K(:,m) % Salinity
 do l = 1, n_Channels
   do jlevel = 1, n_Levels
     jac_t_tmp(l,jlevel) = atm_K(l,m)%Temperature(jlevel)
     do jabs = 1, n_Absorbers
       jac_a_tmp(l,jlevel,jabs) = atm_K(l,m)%Absorber(jlevel,jabs)
     end do
   end do
 end do

end subroutine ufo_reconradop_crtm_rr_jacobian_tmp_populate

! ------------------------------------------------------------------------------
subroutine ufo_reconradop_crtm_normal_jacobian_convert(l, m, n_Sensor, n_Levels, n_Absorbers, channels, rad, sfc_K, rts_K, atm_K)
implicit none
integer, intent(in) :: l ,m, n_Sensor, channels(:), n_Absorbers, n_levels
real(c_double),intent(in) :: rad
type(CRTM_RTSolution_type), intent(inout):: rts_K(:,:)
type(CRTM_Atmosphere_type), intent(inout) :: atm_K(:,:)
type(CRTM_Surface_type),  intent(inout) :: sfc_K(:,:)
integer:: jlevel, jabs
 ! fill in RR operator Tb surface Jacobians

 call CRTM_Planck_Temperature_TL(n_Sensor, channels(l), rad, sfc_K(l,m) % Water_Temperature, sfc_K(l,m) % Water_Temperature)
 call CRTM_Planck_Temperature_TL(n_Sensor, channels(l), rad, sfc_K(l,m) % Land_Temperature, sfc_K(l,m) % Land_Temperature)
 call CRTM_Planck_Temperature_TL(n_Sensor, channels(l), rad, sfc_K(l,m) % Ice_Temperature, sfc_K(l,m) % Ice_Temperature)
 call CRTM_Planck_Temperature_TL(n_Sensor, channels(l), rad, sfc_K(l,m) % Snow_Temperature , sfc_K(l,m) % Snow_Temperature)
 call CRTM_Planck_Temperature_TL(n_Sensor, channels(l), rad, rts_K(l,m) % Surface_Emissivity, rts_K(l,m) % Surface_Emissivity)
 call CRTM_Planck_Temperature_TL(n_Sensor, channels(l), rad, sfc_K(l,m) % Wind_Speed, sfc_K(l,m) % Wind_Speed)
 call CRTM_Planck_Temperature_TL(n_Sensor, channels(l), rad, sfc_K(l,m) % Wind_Direction, sfc_K(l,m) % Wind_Direction)
 call CRTM_Planck_Temperature_TL(n_Sensor, channels(l), rad, sfc_K(l,m) % Salinity, sfc_K(l,m) % Salinity)

 do jlevel = 1,n_Levels
   ! fill RR operator for Tb air temperature Jacobians
   call CRTM_Planck_Temperature_TL(n_Sensor, channels(l), rad, atm_K(l,m) % Temperature(jlevel) , atm_K(l,m) % Temperature(jlevel))
   ! fill RR operator for Tb constituent profile Jacobians
   do jabs = 1,n_Absorbers
     call CRTM_Planck_Temperature_TL(n_Sensor, channels(l), rad, atm_K(l,m) % Absorber(jlevel,jabs), atm_K(l,m) % Absorber(jlevel,jabs))
   end do
 end do
end subroutine ufo_reconradop_crtm_normal_jacobian_convert
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
subroutine ufo_reconradop_crtm_rr_jacobian_convert(l, m, n_Sensor, n_Levels, n_Absorbers,  channels, rad, jac_s_tmp, jac_t_tmp, jac_a_tmp, sfc_K, rts_K, atm_K)
implicit none
integer, intent(in) :: l ,m, n_Sensor, channels(:), n_Levels, n_Absorbers
real(c_double),intent(in) :: rad, jac_s_tmp(:,:), jac_t_tmp(:,:), jac_a_tmp(:,:,:)
type(CRTM_RTSolution_type), intent(inout):: rts_K(:,:)
type(CRTM_Atmosphere_type), intent(inout) :: atm_K(:,:)
type(CRTM_Surface_type),  intent(inout) :: sfc_K(:,:)
integer:: jlevel, jabs
 ! fill in RR operator Tb surface Jacobians

 call CRTM_Planck_Temperature_TL(n_Sensor, channels(l), rad, jac_s_tmp(l,1), sfc_K(l,m) % Water_Temperature)
 call CRTM_Planck_Temperature_TL(n_Sensor, channels(l), rad, jac_s_tmp(l,2), sfc_K(l,m) % Land_Temperature)
 call CRTM_Planck_Temperature_TL(n_Sensor, channels(l), rad, jac_s_tmp(l,3), sfc_K(l,m) % Ice_Temperature)
 call CRTM_Planck_Temperature_TL(n_Sensor, channels(l), rad, jac_s_tmp(l,4), sfc_K(l,m) % Snow_Temperature)
 call CRTM_Planck_Temperature_TL(n_Sensor, channels(l), rad, jac_s_tmp(l,5), rts_K(l,m) % Surface_Emissivity)
 call CRTM_Planck_Temperature_TL(n_Sensor, channels(l), rad, jac_s_tmp(l,6), sfc_K(l,m) % Wind_Speed)
 call CRTM_Planck_Temperature_TL(n_Sensor, channels(l), rad, jac_s_tmp(l,7), sfc_K(l,m) % Wind_Direction)
 call CRTM_Planck_Temperature_TL(n_Sensor, channels(l), rad, jac_s_tmp(l,8), sfc_K(l,m) % Salinity)

 do jlevel = 1,n_Levels
   ! fill RR operator for Tb air temperature Jacobians
   call CRTM_Planck_Temperature_TL(n_Sensor, channels(l), rad, jac_t_tmp(l,jlevel) , atm_K(l,m) % Temperature(jlevel))
   ! fill RR operator for Tb constituent profile Jacobians
   do jabs = 1,n_Absorbers
     call CRTM_Planck_Temperature_TL(n_Sensor, channels(l), rad, jac_a_tmp(l,jlevel,jabs), atm_K(l,m) % Absorber(jlevel,jabs))
   end do
 end do
end subroutine ufo_reconradop_crtm_rr_jacobian_convert
! ------------------------------------------------------------------------------



subroutine ufo_reconradop_crtm_delete(self, Cmatrix_bias, Cmatrix)
implicit none

class(ufo_reconradop_crtm),   intent(in)     :: self
real(kind_real), allocatable, intent(inout)  :: Cmatrix_bias(:)         ! C matrix bias vector used to simulate reconstructed radiances.
real(kind_real), allocatable, intent(inout)  :: Cmatrix(:,:)            ! C matrix used to simulate reconstructed radiances.

 if(allocated(Cmatrix_bias)) deallocate(Cmatrix_bias)
 if(allocated(Cmatrix)) deallocate(Cmatrix)

end subroutine ufo_reconradop_crtm_delete

end module ufo_reconradop_crtm_mod
