! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Fortran module to handle adt observations

module ufo_adt_mod
  
use ioda_obs_adt_mod
use ioda_obs_vectors
use ufo_vars_mod
use ioda_locs_mod
use ufo_geovals_mod
use kinds
use ncd_kinds, only:  i_kind,r_single,r_kind,r_double
  
implicit none
public :: ufo_adt
public :: ufo_adt_eqv
private
integer, parameter :: max_string=800

!> Fortran derived type for adt observation operator
type :: ufo_adt
end type ufo_adt

  type diag_adt_header
    character(3),    dimension(:), allocatable :: ObsType
    integer(i_kind)                            :: n_ObsType
    integer(i_kind)                            :: n_Observations_Tracer
    integer(i_kind)                            :: date
 end type diag_adt_header

  type diag_adt_tracer
    integer      :: Station_ID
    real(r_kind) :: Observation_Type
    real(r_kind) :: Latitude
    real(r_kind) :: Longitude
    real(r_kind) :: Station_Depth
    real(r_kind) :: Time
    real(r_kind) :: Observation
    real(r_kind) :: Obs_Minus_Forecast_adjusted
    real(r_kind) :: Obs_Minus_Forecast_unadjusted
 end type diag_adt_tracer



! ------------------------------------------------------------------------------

contains
 
! ------------------------------------------------------------------------------

subroutine ufo_adt_eqv(self, geovals, hofx, obs_adt)
implicit none
type(ufo_adt), intent(in) :: self
type(ufo_geovals), intent(in)    :: geovals
type(ioda_obs_adt), intent(in) :: obs_adt     !< adt observations
type(obs_vector),  intent(inout) :: hofx

character(len=*), parameter :: myname_="ufo_adt_eqv"
character(max_string) :: err_msg

! nc_diag stuff
logical :: append
character(len=120) :: filename !< name of outpu file for omf, lon, lat, ...
type(diag_adt_tracer), allocatable :: adt_out(:)
type(diag_adt_header) :: adt_out_header


integer :: iobs 
real :: sum_obs,sum_hofx
type(ufo_geoval), pointer :: adt 

! check if nobs is consistent in geovals & hofx
if (geovals%nobs /= hofx%nobs) then
  write(err_msg,*) myname_, ' error: nobs inconsistent!',geovals%nobs,hofx%nobs
  call abor1_ftn(err_msg)
endif

! check if adt variable is in geovals and get it
if (.not. ufo_geovals_get_var(geovals, var_abs_topo, adt)) then
  write(err_msg,*) myname_, trim(var_abs_topo), ' doesnt exist'
  call abor1_ftn(err_msg)
endif

! Information for temporary output file ---------------------------------------!
allocate(adt_out(hofx%nobs))
allocate(adt_out_header%ObsType(1))
adt_out_header%ObsType(1)            = 'aaa'
adt_out_header%n_ObsType             = 1
adt_out_header%n_Observations_Tracer = hofx%nobs
adt_out_header%date                  = 20120729
! -----------------------------------------------------------------------------!


hofx%values = 0.0
sum_hofx=sum(adt%vals(1,:))
sum_obs=sum(obs_adt%adt(:))
print *,'offset',(sum_obs-sum_hofx)/hofx%nobs
! adt obs operator
do iobs = 1, hofx%nobs
!   hofx%values(iobs) = adt%vals(1,iobs)
! remove offset from hofx
   hofx%values(iobs) = adt%vals(1,iobs)+(sum_obs-sum_hofx)/hofx%nobs
   write(302,*)hofx%values(iobs)

  ! Output information:
  adt_out(iobs)%Station_ID                     = 1 
  adt_out(iobs)%Observation_Type               = 1.0
  adt_out(iobs)%Latitude                       = obs_adt%lat(iobs)
  adt_out(iobs)%Longitude                      = obs_adt%lon(iobs)
  adt_out(iobs)%Station_Depth                  = 0 
  adt_out(iobs)%Time                           = 1.0
  adt_out(iobs)%Observation                    = obs_adt%adt(iobs)
  adt_out(iobs)%Obs_Minus_Forecast_adjusted    = obs_adt%adt(iobs) - hofx%values(iobs)
  adt_out(iobs)%Obs_Minus_Forecast_unadjusted  = obs_adt%adt(iobs) - hofx%values(iobs)

!   print *,'hofx is',hofx%values(iobs)
enddo

filename='test.nc'
call tonc(filename,adt_out)

deallocate(adt_out)


end subroutine ufo_adt_eqv

! ------------------------------------------------------------------------------
subroutine tonc(filename,adt_out)
! ------------------------------------------------------------------------------
use netcdf
use ufo_vars_mod

implicit none

character(len=120) :: filename                            !< name of output file for omf, lon, lat, ...
type(diag_adt_tracer) :: adt_out(:)                         !< profile observations

integer :: iobs

!netcdf stuff
integer(kind=4) :: iNcid,i,iStationNo
integer(kind=4) :: iDimStation_ID,iDimTime_ID,iDimLenStringName_ID, iDimLev_ID
integer(kind=4) :: iVarLON_ID, iVarLAT_ID, iVarLev_ID, iVarObs_ID
integer(kind=4) :: iVarOMF_ID, iVarOMA_ID, iVarSIGO_ID, iVarOBSID_ID

integer, allocatable :: obsid(:)
integer :: nlev,nobs

! Create file.
call check(nf90_create(filename, NF90_CLOBBER, iNcid))
! Define the dimensions. The Station-record dimension is defined to have
call check(nf90_def_dim(iNcid, "nobs", NF90_UNLIMITED, iDimStation_ID))

! Define of variables.
! Obs space
call check( nf90_def_var(iNcid, "OBSID", NF90_INT, (/ iDimStation_ID /), iVarOBSID_ID) )
call check( nf90_def_var(iNcid, "lon", NF90_REAL, (/ iDimStation_ID /), iVarLON_ID) )
call check( nf90_def_var(iNcid, "lat", NF90_REAL, (/ iDimStation_ID /), iVarLAT_ID) )
call check( nf90_def_var(iNcid, "lev", NF90_REAL, (/ iDimStation_ID /), iVarLev_ID) )
call check( nf90_def_var(iNcid, "obs", NF90_REAL, (/ iDimStation_ID /), iVarObs_ID) )
call check( nf90_def_var(iNcid, "sigo", NF90_REAL, (/ iDimStation_ID /), iVarSIGO_ID) )
call check( nf90_def_var(iNcid, "omf", NF90_REAL, (/ iDimStation_ID /), iVarOMF_ID) )
call check( nf90_def_var(iNcid, "oma", NF90_REAL, (/ iDimStation_ID /), iVarOMA_ID) )

! End define mode.
call check(nf90_enddef(iNcid))

! Writing
allocate(obsid(size(adt_out(:),1)))
obsid=1073
call check(nf90_put_var(iNcid, iVarLON_ID , obsid(:)))!adt_out(:)%Station_ID))    
call check(nf90_put_var(iNcid, iVarLON_ID , adt_out(:)%Longitude))
call check(nf90_put_var(iNcid, iVarLAT_ID , adt_out(:)%Latitude))
call check(nf90_put_var(iNcid, iVarLev_ID , adt_out(:)%Station_Depth))
call check(nf90_put_var(iNcid, iVarObs_ID , adt_out(:)%Observation))
call check(nf90_put_var(iNcid, iVarOMF_ID , adt_out(:)%Obs_Minus_Forecast_adjusted))
call check(nf90_put_var(iNcid, iVarOMA_ID , adt_out(:)%Obs_Minus_Forecast_adjusted))   ! fix later

! Close file.
call check(nf90_close(iNcid))

end subroutine tonc

! ------------------------------------------------------------------------------
subroutine check(status)
! ------------------------------------------------------------------------------
use netcdf
IMPLICIT NONE
integer(4), intent ( in) :: status
if(status /= nf90_noerr) then
  print *, trim(nf90_strerror(status))
  stop "Stopped"
end if
end subroutine check

end module ufo_adt_mod
