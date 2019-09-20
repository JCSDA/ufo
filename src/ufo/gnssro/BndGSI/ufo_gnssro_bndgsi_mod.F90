!  
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Fortran module to handle gnssro bending angle observations following 
!> the NCEP/GSI (2018 Aug) implementation

module ufo_gnssro_bndgsi_mod
  use, intrinsic::  iso_c_binding
  use kinds
  use ufo_vars_mod
  use ufo_geovals_mod
  use ufo_geovals_mod_c,   only: ufo_geovals_registry
  use ufo_basis_mod,      only: ufo_basis
  use obsspace_mod
  use missing_values_mod
  use gnssro_mod_conf
  use gnssro_mod_constants
  use fckit_log_module,  only : fckit_log
  use ufo_gnssro_bndgsi_util_mod

  implicit none
  public             :: ufo_gnssro_BndGSI
  private

  !> Fortran derived type for gnssro trajectory
  type, extends(ufo_basis) :: ufo_gnssro_BndGSI
   type(gnssro_conf) :: roconf
  contains
   procedure :: setup     => ufo_gnssro_bndgsi_setup
   procedure :: simobs    => ufo_gnssro_bndgsi_simobs
  end type ufo_gnssro_BndGSI

  contains
! ------------------------------------------------------------------------------
subroutine ufo_gnssro_bndgsi_setup(self, f_conf)

use fckit_configuration_module, only: fckit_configuration 
implicit none
class(ufo_gnssro_BndGSI), intent(inout) :: self
type(fckit_configuration), intent(in)   :: f_conf

call gnssro_conf_setup(self%roconf,f_conf)

end subroutine ufo_gnssro_bndgsi_setup

subroutine ufo_gnssro_bndgsi_simobs(self, geovals, hofx, obss)
  implicit none
  class(ufo_gnssro_bndGSI), intent(in)    :: self
  type(ufo_geovals),        intent(in)    :: geovals
  real(kind_real),          intent(inout) :: hofx(:)
  type(c_ptr), value,       intent(in)    :: obss
  character(len=*), parameter             :: myname  = "ufo_gnssro_bndgsi_simobs"
  character(max_string)                   :: err_msg
  integer                                 :: nrecs
  integer                                 :: nlocs
  integer, parameter                      :: nlevAdd = 13 !num of additional levels on top of exsiting model levels
  integer, parameter                      :: ngrd    = 80 !num of new veritcal grids for bending angle computation
  integer                                 :: iobs, k, igrd, irec, icount
  integer                                 :: nlev, nlev1, nlevExt, nlevCheck
  real(kind_real)                         :: rnlevExt
  type(ufo_geoval), pointer               :: t, q, gph, prs
  real(kind_real), allocatable            :: gesT(:,:), gesZ(:,:), gesP(:,:), gesQ(:,:), gesTv(:,:)
  real(kind_real), allocatable            :: obsLat(:)
  real(kind_real), allocatable            :: obsImpP(:)
  real(kind_real), allocatable            :: obsLocR(:)
  real(kind_real), allocatable            :: obsGeoid(:)
  integer(c_size_t), allocatable          :: obsRecnum(:)
  real(kind_real), allocatable            :: temperature(:)
  real(kind_real)                         :: grids(ngrd)
  real(kind_real)                         :: bendingAngle
  integer                                 :: iflip
  integer,allocatable                     :: nlocs_begin(:)
  integer,allocatable                     :: nlocs_end(:)
  real(c_double)                          :: missing

  write(err_msg,*) myname, ": begin"
  call fckit_log%info(err_msg)

  nlocs   = obsspace_get_nlocs(obss) ! number of observations
  nrecs   = obsspace_get_nrecs(obss) ! number of observations
  write(err_msg,*) myname, ': nlocs from gelvals and hofx, nrecs', geovals%nlocs, nlocs, nrecs
  call fckit_log%info(err_msg)

! check if nobs is consistent in geovals & hofx
  if (geovals%nlocs /= size(hofx)) then
    write(err_msg,*) myname, ': nlocs inconsistent!', geovals%nlocs, size(hofx)
    call abor1_ftn(err_msg)
  endif

  missing = missing_value(missing)

! get variables from geovals
  call ufo_geovals_get_var(geovals, var_ts,  t)         ! air temperature
  call ufo_geovals_get_var(geovals, var_q,   q)         ! specific humidity
  if (self%roconf%vertlayer .eq. "mass") then
    call ufo_geovals_get_var(geovals, var_prs,   prs)       ! pressure
    call ufo_geovals_get_var(geovals, var_z,     gph)       ! geopotential height
  else if (self%roconf%vertlayer .eq. "full") then
    call ufo_geovals_get_var(geovals, var_prsi,  prs)       ! pressure
    call ufo_geovals_get_var(geovals, var_zi,    gph)       ! geopotential height
  else
    call ufo_geovals_get_var(geovals, var_prsi,  prs)       ! pressure
    call ufo_geovals_get_var(geovals, var_zi,    gph)      
    write(err_msg,*) myname,': vertlayer has to be mass of full, '//new_line('a')// &
                            '  will use full layer anyway'
    call fckit_log%info(err_msg)
  end if

  nlev  = t%nval   ! number of model mass levels
  nlev1 = prs%nval ! number of model pressure/height levels 

  allocate(gesP(nlev1,nlocs)) 
  allocate(gesZ(nlev1,nlocs)) 
  allocate(gesT(nlev,nlocs)) 
  allocate(gesTv(nlev,nlocs))
  allocate(gesQ(nlev,nlocs)) 

! copy geovals to local background arrays
  iflip = 0
  if (prs%vals(1,1) .lt. prs%vals(prs%nval,1) ) then
     iflip = 1
     write(err_msg,'(a)')'  ufo_gnssro_bndgsi_simobs:'//new_line('a')//                         &
                         '  Model vertical height profile is in descending order,'//new_line('a')// &
                         '  but bndGSI requires it to be ascending order, need flip'
    call fckit_log%info(err_msg)
    do k=1, nlev
       gesT(k,:) = t%vals(nlev-k+1,:)
       gesQ(k,:) = q%vals(nlev-k+1,:)
       gesTv(k,:)= gesT(k,:)*(1+gesQ(k,:)*(rv_over_rd-1))
    enddo
    do k=1, nlev1
       gesP(k,:) = prs%vals(nlev1-k+1,:)
       gesZ(k,:) = gph%vals(nlev1-k+1,:)
    enddo
  else  ! not flipping
    do k=1, nlev
       gesT(k,:)  = t%vals(k,:)
       gesQ(k,:)  = q%vals(k,:)
       gesTv(k,:) = gesT(k,:)*(1+gesQ(k,:)*(rv_over_rd-1))
    enddo

    do k=1, nlev1
       gesP(k,:) = prs%vals(k,:)
       gesZ(k,:) = gph%vals(k,:)
    enddo
  end if

! if background t and q are on mass layers, 
!    while p and z are on interface layers, take the mean of t and q
!       -- gsi manner
  if ( nlev1 /= nlev ) then  
     do k = nlev, 2, -1
        gesQ(k,:) = half* (gesQ(k,:) + gesQ(k-1,:))
        gesTv(k,:) = half* (gesTv(k,:) + gesTv(k-1,:))
!       PLEASE KEEP this COMMENT:
!       to exactly reproduce gsi, t is converted to tv, tv mean is calcualted,
!       then tv mean is converted to t mean
        gesT(k,:) = gesTv(k,:)/(1+ gesQ(k,:)*(rv_over_rd-1))
!       gesT(k,:) = half* (gesT(k,:) + gesT(k-1,:))
     enddo
  end if

! set obs space struture
  allocate(obsLat(nlocs))
  allocate(obsImpP(nlocs))
  allocate(obsLocR(nlocs))
  allocate(obsGeoid(nlocs))
  allocate(obsRecnum(nlocs))
  allocate(nlocs_begin(nrecs))
  allocate(nlocs_end(nrecs))

  call obsspace_get_db(obss, "MetaData", "latitude",         obsLat)
  call obsspace_get_db(obss, "MetaData", "impact_parameter", obsImpP)
  call obsspace_get_db(obss, "MetaData", "earth_radius_of_curvature", obsLocR)
  call obsspace_get_db(obss, "MetaData", "geoid_height_above_reference_ellipsoid", obsGeoid)
  call obsspace_get_recnum(obss, obsRecnum)

  nlocs_begin=1
  nlocs_end=1
  icount = 1
  do iobs = 1, nlocs-1
    if (obsRecnum(iobs+1) /= obsRecnum(iobs)) then
      icount = icount +1  !counting number of records
      nlocs_end(icount-1)= iobs
      nlocs_begin(icount) = iobs+1
    end if
  end do
  nlocs_end(nrecs)= nlocs
  if (icount /= nrecs) then
    write(err_msg,*) "record number is not consistent :", icount, nrecs
    call fckit_log%info(err_msg)
  end if

  nlevExt   = nlev + nlevAdd
  nlevCheck = min(23, nlev)   !number of levels to check super refaction

! define new integration grids
  do igrd = 0, ngrd-1
     grids(igrd+1) = igrd * ds
  end do 

! bending angle forward model starts

  allocate(temperature(nlocs))
  temperature = missing

  iobs = 0
  rec_loop: do irec = 1, nrecs
    obs_loop: do icount = nlocs_begin(irec), nlocs_end(irec)

      iobs = iobs + 1
      call ufo_gnssro_bndgsi_simobs_single(   &
           obsLat(iobs), obsGeoid(iobs), obsLocR(iobs), obsImpP(iobs), &
           gesZ(:,iobs), gesT(:,iobs), gesQ(:,iobs), gesP(:,iobs), &
           grids, self%roconf%use_compress, &
           nlev, nlev1, nlevExt, nlevAdd, nlevCheck, ngrd, &
           temperature(iobs), bendingAngle) 

      hofx(iobs) = bendingAngle

    end do obs_loop
  end do rec_loop

  if (iobs /= nlocs) then
  write(err_msg,*) myname, ": number of obs are not consistent before and after grouping", nlocs, iobs
  call abor1_ftn(err_msg)
  end if

  write(err_msg,*) myname, ": complete"
  call fckit_log%info(err_msg)

! putting temeprature at obs location to obs space for BackgroundCheck ROGSI
  call obsspace_put_db(obss, "MetaData", "temperature", temperature)

  deallocate(obsLat)
  deallocate(obsImpP)
  deallocate(obsLocR)
  deallocate(obsGeoid)
  deallocate(gesP) 
  deallocate(gesZ) 
  deallocate(gesT) 
  deallocate(gesTv) 
  deallocate(gesQ) 
  deallocate(temperature)
  deallocate(obsRecnum)
  deallocate(nlocs_begin)
  deallocate(nlocs_end)

  write(err_msg,*) myname, ": complete"
  call fckit_log%info(err_msg)
end subroutine ufo_gnssro_bndgsi_simobs
! ------------------------------------------------------------------------------
end module ufo_gnssro_bndgsi_mod
