! (C) Copyright 2017-2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Fortran module for ground-based gnss refractivity ropp forward operator
!> following the ROPP (2018 Aug) implementation

module ufo_groundgnss_ropp_mod

use iso_c_binding
use kinds
use ufo_vars_mod
use ufo_geovals_mod
use ufo_geovals_mod_c, only: ufo_geovals_registry
use ufo_basis_mod,     only: ufo_basis
use obsspace_mod  
use missing_values_mod
use ufo_groundgnss_ropp_utils_mod
use fckit_log_module,  only : fckit_log

implicit none
public             :: ufo_groundgnss_ropp
private

  !> Fortran derived type for ground based gnss trajectory
type, extends(ufo_basis) :: ufo_groundgnss_ROPP
  contains
    procedure :: simobs    => ufo_groundgnss_ropp_simobs
end type ufo_groundgnss_ROPP

contains

! ------------------------------------------------------------------------------
!  Description:
!    simulate the zenith total delay point observation 
!       this is done by integrating refractivities 
!       computed by the ROPP forward operator
!
!  Inputs:
!     obsLat        latitude
!     obsLon        longitude
!     obsHeight     station height
!
!     t             temperature
!     q             specific humidity
!     prs           pressure
!     gph           geopotential height
!     z_sfc         surface geometric height
!     
!  Outputs:
!     model_ztd     simulated zenith total delay (put into HofX)
!
! ------------------------------------------------------------------------------
subroutine ufo_groundgnss_ropp_simobs(self, geovals, hofx, obss)
  use ropp_fm_types, only: State1dFM
  use ropp_fm_types, only: Obs1dRefrac
  use typesizes,     only: wp => EightByteReal
  use datetimetypes, only: dp

  implicit none
  class(ufo_groundgnss_ROPP),  intent(in)    :: self
  type(ufo_geovals),           intent(in)    :: geovals
  real(kind_real),             intent(inout) :: hofx(:)
  type(c_ptr), value,          intent(in)    :: obss
  real(c_double)                             :: missingDouble

  type(State1dFM)                    :: x                    ! ROPP 1d state vector data type
  type(Obs1dRefrac)                  :: y                    ! ROPP 1d Refractivity "observations" vector data type

  character(len=*), parameter        :: myname_="ufo_groundgnss_ropp_simobs"
  real(kind=dp)                      :: ob_time
  integer, parameter                 :: max_string = 800

  character(max_string)              :: err_msg                   ! error message 
  integer                            :: nlev, nobs, iobs          ! number of model level, observations, and counter
  integer, allocatable, dimension(:) :: ichk                      ! quality indicator for observation
  type(ufo_geoval), pointer          :: t, q, prs, gph, z_sfc     ! model state variables
  real(kind_real), allocatable       :: obsLat(:), obsLon(:)      ! observation latitude and longitude
  real(kind_real), allocatable       :: obsHeight(:), obsValue(:) ! observation station height and value
  real(kind_real), allocatable       :: model_z(:)                ! model geometric height
  real(kind_real)                    :: station_phi               ! observation geopotential height
  real(kind_real)                    :: model_ztd                 ! simulated observation (zenith total delay)
  integer                            :: iflip                     ! monotonic increasing height
  logical                            :: l_linear                  ! indicator of refractivity interpolation method
  write(err_msg,*) "TRACE: ufo_groundgnss_ropp_simobs: begin"
  call fckit_log%info(err_msg)

! check if nlocs is consistent in geovals & hofx
  if (geovals%nlocs /= size(hofx)) then
      write(err_msg,*) myname_, ' error: nlocs inconsistent!'
      call abor1_ftn(err_msg)
  endif

! get variables from geovals
  call ufo_geovals_get_var(geovals, var_ts,    t)             ! temperature
  call ufo_geovals_get_var(geovals, var_q,     q)             ! specific humidity
  call ufo_geovals_get_var(geovals, var_prs,   prs)           ! pressure
  call ufo_geovals_get_var(geovals, var_z,     gph)           ! geopotential height
  call ufo_geovals_get_var(geovals, var_sfc_geomz, z_sfc)     ! surface geometric height

  missingDouble = missing_value(missingDouble)

  nlev  = t%nval ! number of model levels
  nobs  = obsspace_get_nlocs(obss)

  if (nobs > 0) then 
     iflip = 0
     if (prs%vals(1,1) .lt. prs%vals(prs%nval,1) ) then
       iflip = 1 
       write(err_msg,'(a)') '  ufo_groundgnss_ropp_simobs:'//new_line('a')//                         &
                            '  Model vertical height profile is in descending order,'//new_line('a')// &
                            '  but ROPP requires it to be ascending order, need flip'
       call fckit_log%info(err_msg)
     end if

   ! set obs space struture
     allocate(obsLon(nobs))
     allocate(obsLat(nobs))
     allocate(obsHeight(nobs))
     allocate(obsValue(nobs))

   ! create array for model geometric heights
     allocate(model_z(nlev))
     model_z(:) = 0.0

     call obsspace_get_db(obss, "MetaData", "longitude",        obsLon) 
     call obsspace_get_db(obss, "MetaData", "latitude",         obsLat) 
     call obsspace_get_db(obss, "MetaData", "stationElevation", obsHeight) 
     call obsspace_get_db(obss, "ObsValue", "zenithTotalDelay", obsValue) 
!    obsValue = 0.0

     allocate(ichk(nlev))
     ichk(:) = 0   ! this will hold QC values for observation from QC flags

     write(err_msg,*) "TRACE: ufo_groundgnss_ropp_simobs: begin observation loop, nobs =  ", nobs
     call fckit_log%info(err_msg)   ! always print

     obs_loop: do iobs = 1, nobs 

       ob_time = 0.0
       l_linear = .False.
       call init_ropp_1d_statevec(ob_time,              &
                                  obsLon(iobs),         &
                                  obsLat(iobs),         &
                                  t%vals(:,iobs),       &
                                  q%vals(:,iobs),       &
                                  prs%vals(:,iobs),     &
                                  gph%vals(:,iobs),     &
                                  nlev,                 &
                                  z_sfc%vals(1,iobs),   &
                                  x, iflip)

       call calc_station_phi(obsLat(iobs), obsHeight(iobs), station_phi)
       call init_ropp_1d_obvec(nlev,            &
                               ichk, ob_time,   &
                               obsLat(iobs),    &
                               obsLon(iobs),    &
                               station_phi,     &
                               x, y)

       call ropp_fm_refrac_1d(x,y)

       call calc_model_z(nlev, obsLat(iobs), y%geop, model_z)
       ! note the scaling of the refractivity by 1.e-6 is done in subroutine after integral
       call gnss_ztd_integral(nlev, y%refrac, model_z, obsHeight(iobs), model_ztd, l_linear)
       ! add error trap ! model_ztd initialized to 0 in integral if 0 is returned something very wrong
       if ( model_ztd == 0.0 ) then
         model_ztd = missingDouble
       end if

       ! matching a print format used in initialization of obvec
       if ( iobs == 1 ) then
         write(err_msg,'(a9,2a11,2a16,2a15)') 'ROPPsim: ', 'lat', 'lon',  &
                           'ob', 'bk', 'station_hgt',  'model_terr'
         call fckit_log%debug(err_msg)   ! print when MAIN_DEBUG=1
       end if
       write(err_msg,'(9x,2f11.3,2f16.6,2f15.3)')                    &
                                           obsLat(iobs),             &
                                           obsLon(iobs),             &
                                         obsValue(iobs),             &
                                              model_ztd,             &
                                        obsHeight(iobs),             &
                                             model_z(1)
       call fckit_log%debug(err_msg)   ! print when MAIN_DEBUG=1

       hofx(iobs) = model_ztd

       ! hack -- handling ropp missing value 
       call ropp_tidy_up_1d(x,y)

     end do obs_loop
 
     deallocate(ichk)
     deallocate(obsLat) 
     deallocate(obsLon)
     deallocate(obsHeight)
     deallocate(obsValue)
  end if ! nobs > 0

  write(err_msg,*) "TRACE: ufo_groundgnss_ropp_simobs: completed"
  call fckit_log%info(err_msg)

end subroutine ufo_groundgnss_ropp_simobs
! ------------------------------------------------------------------------------

end module ufo_groundgnss_ropp_mod
