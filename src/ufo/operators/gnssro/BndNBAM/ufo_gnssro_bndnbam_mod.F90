!  
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
!> Fortran module of Gnssro NBAM (NCEP's Bending Angle Method)
!> nonlinear operator

module ufo_gnssro_bndnbam_mod
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
  use gnssro_mod_transform, only: geop2geometric, compute_refractivity
  use gnssro_mod_grids,  only : get_coordinate_value
  use fckit_log_module,  only : fckit_log
  use ufo_gnssro_bndnbam_util_mod
  use ufo_utils_mod, only: cmp_strings 
  use ufo_constants_mod, only: zero, half, one, two, grav, rd, rv_over_rd

  implicit none
  public             :: ufo_gnssro_BndNBAM
  private

  !> Fortran derived type for gnssro trajectory
  type, extends(ufo_basis) :: ufo_gnssro_BndNBAM
   type(gnssro_conf) :: roconf
  contains
   procedure :: setup     => ufo_gnssro_bndnbam_setup
   procedure :: simobs    => ufo_gnssro_bndnbam_simobs
  end type ufo_gnssro_BndNBAM

  contains
! ------------------------------------------------------------------------------
subroutine ufo_gnssro_bndnbam_setup(self, f_conf)

use fckit_configuration_module, only: fckit_configuration 
implicit none
class(ufo_gnssro_BndNBAM), intent(inout) :: self
type(fckit_configuration), intent(in)   :: f_conf

call gnssro_conf_setup(self%roconf,f_conf)

end subroutine ufo_gnssro_bndnbam_setup

subroutine ufo_gnssro_bndnbam_simobs(self, geovals, hofx, obss)
  implicit none
  class(ufo_gnssro_BndNBAM), intent(in)    :: self
  type(ufo_geovals),        intent(in)    :: geovals
  real(kind_real),          intent(inout) :: hofx(:)
  type(c_ptr), value,       intent(in)    :: obss
  character(len=*), parameter             :: myname  = "ufo_gnssro_bndnbam_simobs"
  character(max_string)                   :: err_msg
  integer                                 :: nrecs, nlocs
  integer, parameter                      :: nlevAdd = 13 !num of additional levels on top of exsiting model levels
  integer, parameter                      :: ngrd    = 80 !num of new veritcal grids for bending angle computation
  integer                                 :: iobs, k, igrd, irec, icount, kk
  integer                                 :: nlev, nlev1, nlevExt, nlevCheck
  type(ufo_geoval), pointer               :: t, q, gph, prs, zs
  real(kind_real), allocatable            :: gesT(:,:), gesZ(:,:), gesP(:,:), gesQ(:,:), gesTv(:,:), gesZs(:)
  real(kind_real), allocatable            :: obsLat(:), obsImpP(:),obsLocR(:), obsGeoid(:), obsValue(:)
  integer(c_size_t), allocatable          :: obsRecnum(:)
  real(kind_real), allocatable            :: temperature(:)
  real(kind_real), allocatable            :: humidity(:),refractivity(:),pressure(:)
  real(kind_real)                         :: temp, geop
  real(kind_real)                         :: wf
  integer                                 :: wi, wi2
  real(kind_real)                         :: grids(ngrd)
  real(kind_real), allocatable            :: refIndex(:), refXrad(:), geomz(:)
  real(kind_real), allocatable            :: ref(:), radius(:)
  real(kind_real)                         :: sIndx
  integer                                 :: indx
  integer                                 :: iflip
  integer,allocatable                     :: nlocs_begin(:)
  integer,allocatable                     :: nlocs_end(:)
  real(c_double)                          :: missing
  integer,          allocatable           :: super_refraction_flag(:), super(:), obs_max(:)
  real(kind_real),  allocatable           :: toss_max(:)
  integer                                 :: sr_hgt_idx
  real(kind_real)                         :: gradRef, obsImpH
  integer,          allocatable           :: LayerIdx(:)

  write(err_msg,*) myname, ": begin"
  call fckit_log%debug(err_msg)

  nlocs   = obsspace_get_nlocs(obss) ! number of observations
  nrecs   = obsspace_get_nrecs(obss) ! number of records/profiles
  write(err_msg,*) myname, ': nlocs from gelvals and hofx, nrecs', geovals%nlocs, nlocs, nrecs
  call fckit_log%debug(err_msg)
  missing = missing_value(missing)

  allocate(temperature(nlocs))
  temperature = missing
  if (trim(self%roconf%output_diags) .eq. "true") then
     allocate(humidity(nlocs))      ! at obs location
     allocate(pressure(nlocs))      ! at obs location
     allocate(refractivity(nlocs))  ! at obs location
     humidity = missing
     pressure = missing
     refractivity = missing
  endif
  allocate(super_refraction_flag(nlocs))
  super_refraction_flag = 0
  allocate(LayerIdx(nlocs))
  LayerIdx = 0

  if (nlocs > 0) then ! check if ZERO OBS

! check if nobs is consistent in geovals & hofx
  if (geovals%nlocs /= size(hofx)) then
    write(err_msg,*) myname, ': nlocs inconsistent!', geovals%nlocs, size(hofx)
    call abor1_ftn(err_msg)
  endif

! get variables from geovals
  call ufo_geovals_get_var(geovals, var_ts,  t)         ! air temperature
  call ufo_geovals_get_var(geovals, var_q,   q)         ! specific humidity
  call ufo_geovals_get_var(geovals, var_sfc_geomz, zs)      ! surface geopotential height/surface altitude

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
  allocate(gesZs(nlocs))

! copy geovals to local background arrays
  iflip = 0
  if (prs%vals(1,1) .lt. prs%vals(prs%nval,1) ) then
     iflip = 1
     write(err_msg,'(a)')'  ufo_gnssro_bndnbam_simobs:'//new_line('a')//                         &
                         '  Model vertical height profile is in descending order,'//new_line('a')// &
                         '  but bndNBAM requires it to be ascending order, need flip'
    call fckit_log%debug(err_msg)
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
       gesZs(:) = zs%vals(1,:)

! if background t and q are on mass layers, 
!    while p and z are on interface layers, take the mean of t and q
!       -- NBAM manner
  if ( nlev1 /= nlev ) then  
     do k = nlev, 2, -1
        gesQ(k,:) = half* (gesQ(k,:) + gesQ(k-1,:))
        gesTv(k,:) = half* (gesTv(k,:) + gesTv(k-1,:))
!       PLEASE KEEP this COMMENT:
!       to exactly reproduce nbam, t is converted to tv, tv mean is calcualted,
!       then tv mean is converted to t mean
        gesT(k,:) = gesTv(k,:)/(1+ gesQ(k,:)*(rv_over_rd-1))
     enddo
  end if

! set obs space struture
  allocate(obsLat(nlocs))
  allocate(obsImpP(nlocs))
  allocate(obsLocR(nlocs))
  allocate(obsGeoid(nlocs))
  allocate(obsValue(nlocs))
  allocate(obsRecnum(nlocs))
  allocate(nlocs_begin(nrecs))
  allocate(nlocs_end(nrecs))

  call obsspace_get_db(obss, "MetaData", "latitude",         obsLat)
  call obsspace_get_db(obss, "MetaData", "impactParameterRO", obsImpP)
  call obsspace_get_db(obss, "MetaData", "earthRadiusCurvature", obsLocR)
  call obsspace_get_db(obss, "MetaData", "geoidUndulation", obsGeoid)
  call obsspace_get_db(obss, "ObsValue", "bendingAngle", obsValue)
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
    write(err_msg,*) myname, ': record number is not consistent :', icount, nrecs
    call abor1_ftn(err_msg)
  end if

  nlevExt   = nlev + nlevAdd
  nlevCheck = min(23, nlev)   !number of levels to check super refaction

! define new integration grids
  do igrd = 0, ngrd-1
     grids(igrd+1) = igrd * ds
  end do 

! bending angle forward model starts
  allocate(geomz(nlev))    ! geometric height
  allocate(radius(nlev))   ! tangent point radisu to earth center
  allocate(ref(nlevExt))   ! refractivity
  allocate(refIndex(nlev))              !refactivity index n
  allocate(refXrad(0:nlevExt+1))        !x=nr, model conuterpart impact parameter

  allocate(super(nlocs))
  allocate(toss_max(nrecs))
  allocate(obs_max(nrecs))

  iobs = 0
  hofx =  missing
  super = 0
  obs_max  = 0
  toss_max = 0

  rec_loop: do irec = 1, nrecs

    obs_loop: do icount = nlocs_begin(irec), nlocs_end(irec)

      iobs = icount
      do k = 1, nlev
!        compute guess geometric height from geopotential height
         call geop2geometric(obsLat(iobs), gesZ(k,iobs)-gesZs(iobs), geomz(k))
         radius(k) = geomz(k) + gesZs(iobs) + obsGeoid(iobs) + obsLocR(iobs)   ! radius r
!        guess refactivity, refactivity index,  and impact parameter
         call compute_refractivity(gesT(k,iobs), gesQ(k,iobs), gesP(k,iobs),   &
                                ref(k), self%roconf%use_compress)
         refIndex(k) = one + (r1em6*ref(k))
         refXrad(k)  = refIndex(k) * radius(k)
      end do

!     Data rejection based on model background !
!     (1) skip data beyond model levels
      call get_coordinate_value(obsImpP(iobs),sIndx,refXrad(1),nlev,"increasing")
      if (sIndx < one .or. sIndx > float(nlev))  cycle obs_loop

!     save the obs vertical location index (unit: model layer)
      LayerIdx(iobs) = min(max(1, int(sIndx)), nlev)

!     calculating virtual temperature at obs location to obs space for BackgroundCheck RONBAM
      indx=sIndx
      wi=min(max(1,indx),nlev)
      wi2=max(1,min(indx+1,nlev))
      wf=sIndx-float(wi)
      wf=max(zero,min(wf,one))
      temperature(iobs)=gesTv(wi,iobs)*(one-wf)+gesTv(wi2,iobs)*wf
      if (trim(self%roconf%output_diags) .eq. "true") then
         humidity(iobs)= gesQ(wi,iobs)*(one-wf)+gesQ(wi2,iobs)*wf 
         temp          = gesT(wi,iobs)*(one-wf)+gesT(wi2,iobs)*wf 
         geop          = gesZ(wi,iobs)*(one-wf)+gesZ(wi2,iobs)*wf
         pressure(iobs)= gesP(wi,iobs)/exp(two*grav*(geop-gesZ(wi,iobs))/(rd*(temperature(iobs)+gesTv(wi,iobs))))
         call compute_refractivity(temp, humidity(iobs), pressure(iobs),   &
                                 refractivity(iobs), self%roconf%use_compress)
      end if

!     (2) super-refaction
!     (2.1) GSI style super refraction check
      if(cmp_strings(self%roconf%super_ref_qc, "NBAM")) then

        obsImpH = (obsImpP(iobs) - obsLocR(iobs)) * r1em3 !impact heigt: a-r_earth

        if (obsImpH <= six) then
           do k = nlevCheck, 1, -1

!             N gradient
              gradRef = 1000.0 * (ref(k+1)-ref(k))/(radius(k+1)-radius(k))
!             check for model SR layer
              if (abs(gradRef) >= 0.75*crit_gradRefr .and. obsImpP(iobs) <= refXrad(k+2)) then
                  super_refraction_flag(iobs) = 1
                  cycle obs_loop
              endif
           end do

!          relax to close-to-SR conditions, and check if obs is inside model SR layer
           
           if (self%roconf%sr_steps > 1                 &
              .and. obsValue(iobs) >= 0.03 ) then

               do k = nlevCheck, 1, -1
                  gradRef = 1000.0 * (ref(k+1)-ref(k))/(radius(k+1)-radius(k))
                  if (abs(gradRef) >= half*crit_gradRefr & 
                     .and. super(iobs) == 0                   &
                     .and. toss_max(irec) <= obsValue(iobs)) then
                      obs_max(irec) = iobs
                      toss_max(irec)= max(toss_max(irec), obsValue(iobs))
                      super(iobs) = 1
                  end if
              end do ! k
 
           end if   ! end if(self%roconf%sr_steps > 1 
        end if ! obsImpH <= six

!    ROPP style super refraction check
     else if(cmp_strings(self%roconf%super_ref_qc, "ECMWF")) then

       sr_hgt_idx = 1
       do k = nlev, 2, -1
          if (refXrad(k) - refXrad(k-1) < 10.0) THEN
             sr_hgt_idx = k
             exit
          end if
       end do

       if (obsImpP(iobs) < refXrad(sr_hgt_idx)) then
          super_refraction_flag(iobs) = 1
          cycle obs_loop
       end if

     else
       write(err_msg,*) myname, ': super refraction method has to be NBAM or ECMWF!'
       call abor1_ftn(err_msg)
     end if

     if (super_refraction_flag(iobs) .eq. 0) then
     call ufo_gnssro_bndnbam_simobs_single(   &
               obsLat(iobs), obsGeoid(iobs), obsLocR(iobs), obsImpP(iobs), &
               grids, ngrd, &
               nlev, nlevExt, nlevAdd, nlevCheck, &
               radius(1:nlev),ref(1:nlevExt),refIndex(1:nlev),refXrad(0:nlevExt),  &
               hofx(iobs))

     end if
    end do obs_loop
  end do rec_loop

  if (cmp_strings(self%roconf%super_ref_qc, "NBAM") .and. self%roconf%sr_steps > 1 ) then
     rec_loop2: do irec = 1, nrecs

       if (obs_max(irec) > 0 ) then

          obs_loop2: do k = nlocs_begin(irec), nlocs_end(irec)
            obsImpH = (obsImpP(k) - obsLocR(k)) * r1em3
            if (obsImpH<=six .and. obsImpP(k)<=obsImpP(obs_max(irec)) .and.  &
                hofx(k)/=missing .and. super_refraction_flag(k)==0) then
                super_refraction_flag(k)=2
            end if
          end do obs_loop2

       end if

     end do rec_loop2
  end if

  deallocate(obsLat)
  deallocate(obsImpP)
  deallocate(obsLocR)
  deallocate(obsGeoid)
  deallocate(gesP) 
  deallocate(gesZ) 
  deallocate(gesT) 
  deallocate(gesTv) 
  deallocate(gesQ)
  deallocate(gesZs) 
  deallocate(ref)
  deallocate(refIndex)
  deallocate(refXrad)
  deallocate(geomz)
  deallocate(radius)
  deallocate(obsRecnum)
  deallocate(nlocs_begin)
  deallocate(nlocs_end)
  deallocate(super)
  deallocate(toss_max)
  deallocate(obs_max)

  write(err_msg,*) myname, ": complete"
  call fckit_log%debug(err_msg)
  end if ! end check if ZERO OBS

! putting virtual temeprature at obs location to obs space for BackgroundCheck RONBAM
  call obsspace_put_db(obss, "MetaData", "virtual_temperature", temperature)
! putting super refraction flag to obs space 
  call obsspace_put_db(obss, "SRflag",   "bending_angle", super_refraction_flag)
! saving obs vertical model layer postion for later
  call obsspace_put_db(obss, "LayerIdx",   "bending_angle", LayerIdx)
  if (trim(self%roconf%output_diags) .eq. "true") then
      call obsspace_put_db(obss, "ObsDiag", "specific_humidity", humidity)
      call obsspace_put_db(obss, "ObsDiag", "refractivity", refractivity)
      call obsspace_put_db(obss, "ObsDiag", "pressure", pressure)
      deallocate(humidity)
      deallocate(pressure)
      deallocate(refractivity)
  end if
  deallocate(super_refraction_flag)
  deallocate(temperature)
  deallocate(LayerIdx)

end subroutine ufo_gnssro_bndnbam_simobs
! ------------------------------------------------------------------------------
end module ufo_gnssro_bndnbam_mod
