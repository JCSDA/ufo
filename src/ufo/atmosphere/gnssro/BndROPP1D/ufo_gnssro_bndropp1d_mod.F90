! (C) Copyright 2017-2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Fortran module to handle gnssro bending angle observations following 
!> the ROPP (2018 Aug) implementation

module ufo_gnssro_bndropp1d_mod
  use iso_c_binding
  use ufo_vars_mod
  use ufo_geovals_mod
  use ufo_geovals_mod_c, only: ufo_geovals_registry
  use ufo_basis_mod,     only: ufo_basis
  use vert_interp_mod
  use lag_interp_mod, only: lag_interp_const, lag_interp_smthWeights
  use obsspace_mod  

  use kinds
  implicit none
  public             :: ufo_gnssro_bndropp1d
  private

  !> Fortran derived type for gnssro trajectory
  type, extends(ufo_basis) :: ufo_gnssro_BndROPP1D
  contains
    procedure :: simobs    => ufo_gnssro_bndropp1d_simobs
  end type ufo_gnssro_BndROPP1D

  contains
! ------------------------------------------------------------------------------
  subroutine ufo_gnssro_bndropp1d_simobs(self, geovals, hofx, obss)
    use ropp_fm_types, only: State1dFM
    use ropp_fm_types, only: Obs1dBangle
    use datetimetypes, only: dp
      implicit none
      class(ufo_gnssro_BndROPP1D), intent(in):: self
      type(ufo_geovals), intent(in)          :: geovals
      real(kind_real),   intent(inout)       :: hofx(:)
      type(c_ptr), value, intent(in)         :: obss

      type(State1dFM)                 :: x
      type(Obs1dBangle)               :: y

      character(len=*), parameter     :: myname_="ufo_gnssro_bndropp1d_simobs"
      real(kind=dp)                   :: ob_time
      integer, parameter              :: max_string = 800

      character(max_string)           :: err_msg
      character(len=250)              :: record
      integer                         :: iobs
      integer                         :: nlev, nobs
      integer                         :: nvprof
      integer, allocatable, dimension(:)      :: ichk
      type(ufo_geoval), pointer       :: t, q, prs, z, z_sfc
      real(kind_real), allocatable    :: obsLat(:), obsLon(:), obsImpP(:), obsLocR(:), obsGeoid(:)
      real(kind_real), allocatable    :: obsYYYY(:), obsMM(:), obsDD(:), obsHH(:), obsMN(:), obsSS(:)
      integer :: obss_nobs

      write(*,*) "TRACE: ufo_gnssro_bndropp1d_simobs: begin"
      ! check if nobs is consistent in geovals & hofx
      if (geovals%nobs /= size(hofx)) then
        write(err_msg,*) myname_, ' error: nobs inconsistent!'
        call abor1_ftn(err_msg)
      endif
      ! check if prs (pressure at model levels)  variable is in geovals and get it
      call ufo_geovals_get_var(geovals, var_prs, prs)
      ! check if specific humidity variable is in geovals and get it
      call ufo_geovals_get_var(geovals, var_q, q)
      ! check if geopotential height variable is in geovals and get it
      call ufo_geovals_get_var(geovals, var_z, z)
      ! check if surface geopotential height variable is in geovals and get it
      call ufo_geovals_get_var(geovals, var_sfc_z, z_sfc)
      ! check if t variable is in geovals and get it
      call ufo_geovals_get_var(geovals, var_t, t)

      nlev  = q%nval ! number of model levels
      nobs  = geovals%nobs ! number of observations

      ! read observation vectors
      obss_nobs = obsspace_get_nobs(obss)
      !allocate(obsYYYY(obss_nobs))
      !allocate(obsMM(obss_nobs))
      !allocate(obsDD(obss_nobs))
      !allocate(obsHH(obss_nobs))
      !allocate(obsMN(obss_nobs))
      !allocate(obsSS(obss_nobs))
      allocate(obsLon(obss_nobs))
      allocate(obsLat(obss_nobs))
      allocate(obsImpP(obss_nobs))
      allocate(obsLocR(obss_nobs))
      allocate(obsGeoid(obss_nobs))

     !call obsspace_get_db(obss, "Metadata", "year", obsYYYY) !Note: variable name not consistent with BUFR table
     !call obsspace_get_db(obss, "Metadata", "month", obsMM) !Note: variable name not consistent with BUFR table
     !call obsspace_get_db(obss, "Metadata", "day", obsDD) !Note: variable name not consistent with BUFR table
     !call obsspace_get_db(obss, "Metadata", "hour", obsHH) !Note: variable name not consistent with BUFR table
     !call obsspace_get_db(obss, "Metadata", "minute", obsMN) !Note: variable name not consistent with BUFR table
     !call obsspace_get_db(obss, "Metadata", "second", obsSS) !Note: variable name not consistent with BUFR table
      call obsspace_get_db(obss, "Metadata", "Longitude", obsLon) !Note: variable name not consistent with BUFR table
      call obsspace_get_db(obss, "Metadata", "Latitude", obsLat) !Note: variable name not consistent with BUFR table
      call obsspace_get_db(obss, "Metadata", "IMPP", obsImpP) !observed impact parameter
      call obsspace_get_db(obss, "Metadata", "ELRC", obsLocR) !local radius of earth. Note: need add to test data
      call obsspace_get_db(obss, "Metadata", "GEODU", obsGeoid) !Geoid. Note: need add to test data

      nvprof = 1 ! number of vertical profiles (occultation points)
      allocate(ichk(nvprof))
      ichk(:) = 0   ! this will hold QC values for observation from QC flags

      write(record,*) "DEBUG: ufo_gnssro_bndropp1d_simobs: begin observation loop  ", nobs
      obs_loop: do iobs = 1, nobs 

!         call init_ob_time(int(obsYYYY(iobs)),   &
!                           int(obsMM(iobs)),     &
!                           int(obsDD(iobs)),     &
!                           int(obsHH(iobs)),     &
!                           int(obsMN(iobs)),     &
!                           int(obsSS(iobs)),     &
!                           ob_time)
                            ob_time = 0.0
         ! alternatively you can use the analysis time and the observation value of time
         call init_ropp_1d_statevec(ob_time,               &
                                    obsLon(iobs),   &
                                    obsLat(iobs),   &
                                    t%vals(:,iobs),        &
                                    q%vals(:,iobs),        &
                                    prs%vals(:,iobs),      &
                                    z%vals(:,iobs),        &
                                    nlev,                  &
                                    z_sfc%vals(1,iobs),    &
                                    x)
         call init_ropp_1d_obvec(nvprof,                &
                                 obsImpP(iobs),  &
                                 ichk, ob_time,         &
                                 obsLat(iobs),   &
                                 obsLon(iobs),   &
                                 obsLocR(iobs),  &  
                                 obsGeoid(iobs), &
                                 y)

         call ropp_fm_bangle_1d(x,y)

         hofx(iobs) = y%bangle(nvprof)  ! nvprof is just one point

      end do obs_loop
      
      deallocate(obsLat) !Note: to be removed
      deallocate(obsLon)
      deallocate(obsImpP)
      deallocate(obsLocR)
      deallocate(obsGeoid)
      !deallocate(obsYYYY)
      !deallocate(obsMM)
      !deallocate(obsDD)
      !deallocate(obsHH)
      !deallocate(obsMN)
      !deallocate(obsSS)
      write(*,*) "TRACE: ufo_gnssro_bndropp1d_simobs: completed"

   end subroutine ufo_gnssro_bndropp1d_simobs
! ------------------------------------------------------------------------------

      subroutine init_ropp_1d_statevec(step_time,rlon,rlat, &
             temp,shum,pres,phi,lm,phi_sfc,x)

!  Description:
!     subroutine to fill a ROPP state vector structure with
!     model background fields: Temperature, pressure, specific
!     humidity at full-levels, and surface geopotential height.
!
!  Inputs:
!     temp   background temperature
!     shum   background specific humidity
!     pres   background pressure
!     phi    geopotential height
!
!     phi_sfc  terrain geopotential of background field
!     lm       number of vertical levels in  the background
!
!  Outputs:
!     x      Forward model state vector
!
! ###############################################################

! For ROPP data type and library subroutine

      use typesizes, only: wp => EightByteReal
      use datetimetypes, only: dp
      use kinds, only: kind_real
      use ropp_fm_types, only: State1dFM
      use geodesy,   only: gravity, R_eff, geometric2geopotential
      use arrays, only: callocate

      implicit none

! Function arguments

!     Output state vector
      type(State1dFM),      intent(out)   :: x
      real(kind=dp),        intent(in)    :: step_time
      real(kind=kind_real), intent(in)    :: rlat, rlon
      real(kind=kind_real), intent(in)    :: phi_sfc
      integer,              intent(in)    :: lm
      real(kind=kind_real), dimension(lm), intent(in)  :: temp,shum,pres,phi

! Local variables
      character(len=250)                 :: record
      real(kind=kind_real)               :: rlon_local
      integer n,i,j,k

!-------------------------------------------------------------------------
      x%state_ok = .TRUE.
      x%new_bangle_op = .TRUE.     ! activate ROPP v8 new interpolation scheme

! ROPP Longitude value is -180.0 to 180.0

      x%lat = real(rlat,kind=wp)
      rlon_local = rlon
      if (rlon_local .gt. 180) rlon_local = rlon_local - 360.
      x%lon = real(rlon_local,kind=wp)
      x%time = real(step_time,kind=wp)


! Number of levels in background profile.  What about (lm+1) field ?

      x%n_lev = lm

!--------------------------------------------------------
! allocate arrays for temperature, specific humidity, pressure
! and geopotential height data
!---------------------------------------------------------
      allocate(x%temp(x%n_lev))
      allocate(x%shum(x%n_lev))
      allocate(x%pres(x%n_lev))
      allocate(x%geop(x%n_lev))
!----------------------------------------------------
! ROPP FM requires vertical height profile to be of the ascending order.
! (see ropp_io_ascend ( ROdata )). So we need to flip the data.
!----------------------------------------------------
      n = lm
      write(record,'(4a9,a11)') 'lvl','temp','shum','pres','geop'
      do k = 1, lm
         x%temp(n) = real(temp(k),kind=wp)
         x%shum(n) = real(shum(k),kind=wp)
         x%pres(n) = real(pres(k)*100.,kind=wp)
         x%geop(n) = real(phi(k),kind=wp)
         write(record,'(5x,i4,f9.2,f9.4,f9.1,f15.1)') &
              n, x%temp(n), x%shum(n), x%pres(n), x%geop(n)
         n = n - 1
      end do

! sufrace geopotential height value

      x%geop_sfc = real(phi_sfc,kind=wp)
      write(record,'("geop_sfc",f15.2)') x%geop_sfc

!------------------------------------------------
! covariance matrix, is this used by ROPP FM?
!------------------------------------------------
      x%cov_ok = .TRUE.
 
! Allocate memory
! For ECMWF example, Covariance matrix for temperature sigma and
! specific humidity sigma, and surface pressure. There the
! size of covariance matrix is 2 * nlevel + 1.
 
      n = (2*x%n_lev)+1  ! Number of elements in the state vector

      if (associated(x%cov%d)) deallocate(x%cov%d)
      call callocate(x%cov%d, n*(n+1)/2)    ! From ROPP utility library
!     allocate(x%cov%d((n*(n+1))/2))        ! or just use standard Fortran call

      do i = 1, x%n_lev
         x%cov%d(i + i*(i-1)/2) = 1.0_wp
      end do

      do i = 1, x%n_lev
         j = x%n_lev + i
         x%cov%d(j + j*(j-1)/2) = 1.0_wp
      enddo

      x%cov%d(n + n*(n-1)/2) = 1.0_wp

!-------------------------------------------------------------
! Rest of the covariance marix
!--------------------------------------------------------------

      if (associated(x%cov%e)) deallocate(x%cov%e)
      if (associated(x%cov%f)) deallocate(x%cov%f)
      if (associated(x%cov%s)) deallocate(x%cov%s)

      x%cov%fact_chol = .FALSE.
      x%cov%equi_chol = 'N'

      return
      end subroutine init_ropp_1d_statevec

      subroutine init_ropp_1d_obvec(nvprof,obs_impact,  &
                         ichk,ob_time,rlat,rlon,roc,undulat,y)

!  Description:
!     subroutine to fill a ROPP observation vector structure
!     observation provides the inputs of:
!     impact parameter, lat, lon, time, radius of curvature
!
!     forward model will provide the simulation of bending angle
!
!  Inputs:
!     obs_impact    impact parameter
!     ob_time       time of the observation
!     rlat          latitude
!     rlon          longitude
!     roc           radius of curvature
!     undulat       undulation correction for radius of curvature
!
!  Outputs:
!     y:     Partially filled Forward model observation vector
!
! ###############################################################

! For ROPP data type

      use typesizes, only: wp => EightByteReal
      use kinds, only: kind_real
      use datetimetypes, only: dp
      use ropp_fm_types, only: Obs1dBangle
      use geodesy,   only: gravity, R_eff, geopotential2geometric

      implicit none


!     Output state vector
      type(Obs1dBangle),          intent(out)   :: y

      integer,                    intent(in)    :: nvprof
      integer, dimension(nvprof), intent(in)    :: ichk
      real(kind=kind_real), dimension(nvprof), intent(in)    :: obs_impact
      real(kind=kind_real),       intent(in)    :: rlat, rlon
      real(kind=kind_real),       intent(in)    :: roc, undulat
      real(kind=dp),              intent(in)    :: ob_time

      real(kind=wp)                             :: r8lat
      real(kind=kind_real)                      :: rlon_local
      character(len=250)                        :: record

      integer                                   :: i

!-------------------------------------------------------------------------

      y%time = real(ob_time,kind=wp)
      r8lat = real(rlat,kind=wp)
      y%lat = r8lat
      rlon_local = rlon
      if (rlon_local .gt. 180) rlon_local = rlon_local - 360.
      y%lon = real(rlon_local,kind=wp)
      y%nobs   = nvprof
      y%g_sfc = gravity(r8lat, 0.0_wp)          ! 2nd argument is height(m) above the sfc
      y%r_curve = real(roc,kind=wp)             ! ROPP rad of curve (m)
      y%undulation = real(undulat,kind=wp)      ! ROPP undulation corr for rad of curve (m)
      y%r_earth = r_eff(r8lat)

!--------------------------------------------------------
! allocate bending angle, impact parameter & weights
!---------------------------------------------------------
      if (associated(y%bangle)) then
        deallocate(y%bangle)
        deallocate(y%impact)
        deallocate(y%weights)
        nullify(y%bangle)
        nullify(y%impact)
        nullify(y%weights)
      end if

      allocate(y%bangle(1:nvprof))                     ! value computed in fwd model
      allocate(y%impact(1:nvprof))
      allocate(y%weights(1:nvprof))                    ! value set in fwd model

      do i=1,nvprof
        y%impact(i) = real(obs_impact(i),kind=wp)  ! ROPP expects impact parameter in meters
        if (ichk(i) .le. 0) then
          y%weights(i) = 1.0_wp                    ! following t_fascod example
        else
          y%weights(i) = 0.0_wp
        end if
      end do
      y%bangle(:)  = 0.0_wp                    ! following t_fascod example

      write(record,'(a9,2a11,3a15)') 'ROPPyvec:','lat', 'lon',  &
                        'g_sfc', 'roc', 'r_earth_eff'
      write(record,'(9x,2f11.2,f15.6,2f15.2)') y%lat, y%lon,    &
                        y%g_sfc, y%r_curve, y%r_earth

!------------------------------------------------
! covariance matrix, is this used by ROPP FM?
!------------------------------------------------
      y%obs_ok = .TRUE.

      return
      end subroutine init_ropp_1d_obvec

      subroutine init_ob_time(yyyy, mm, dd, hh, mn, ss, ob_time)

      use datetimetypes, only: dp

      integer, intent(in)             :: yyyy, mm, dd, hh, mn, ss
      real(dp), intent(out)           :: ob_time

      integer, dimension(8)           :: dt8

!---------------------------------------------------------------
! Compute Julian seconds from YYYYMMDDHH information of anal time
! (ropp_utils-6.0/datetime/timesince.f90.
!---------------------------------------------------------------
! EXAMPLES
!   1) Current time in seconds since midnight, 1-Jan-2000
!      USE DateTime
!      INTEGER  :: CDT(8)
!      REAL(dp) :: JDF
!      CALL Date_and_Time_UTC ( Values=CDT )
!      --> CDT = 2010,4,9,15,17,14,234
!      CALL TimeSince ( CDT, JDF, 1, "JS2000" )
!      --> JDF = 324141434.23402011
!---------------------------------------------------------------

      dt8(1) = yyyy
      dt8(2) = mm
      dt8(3) = dd
      dt8(4) = 0     ! time zone offset
      dt8(5) = hh
      dt8(6) = mn    ! minute
      dt8(7) = ss    ! second
      dt8(8) = 0     ! millisecond
      ! in the call to timesince the first argument CDT has dimension(:) which
      !  causes an issue with the addressing in gfortran
      !  the time is not used in the operator at this time but this should
      !  be remedied
      !call timesince ( dt8, ob_time, 1, "JS2000" )
      ob_time = 0.

      end subroutine init_ob_time

end module ufo_gnssro_bndropp1d_mod
