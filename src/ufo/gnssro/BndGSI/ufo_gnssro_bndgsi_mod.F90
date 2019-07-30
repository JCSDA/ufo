!  
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Fortran module to handle gnssro bending angle observations following 
!> the NCEP/GSI (2018 Aug) implementation

module ufo_gnssro_bndgsi_mod
  use fckit_configuration_module, only: fckit_configuration 
  use iso_c_binding
  use kinds
  use ufo_vars_mod
  use ufo_geovals_mod
  use ufo_geovals_mod_c,   only: ufo_geovals_registry
  use vert_interp_mod
  use ufo_basis_mod,      only: ufo_basis
  use lag_interp_mod,     only: lag_interp_const, lag_interp_smthWeights
  use obsspace_mod
  use missing_values_mod
  use gnssro_mod_conf
  use gnssro_mod_constants
  use fckit_log_module,  only : fckit_log

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

implicit none
class(ufo_gnssro_BndGSI), intent(inout) :: self
type(fckit_configuration), intent(in)   :: f_conf

call gnssro_conf_setup(self%roconf,f_conf)

end subroutine ufo_gnssro_bndgsi_setup

subroutine ufo_gnssro_bndgsi_simobs(self, geovals, hofx, obss)
  use gnssro_mod_transform
  use gnssro_mod_grids, only: get_coordinate_value
  implicit none
  class(ufo_gnssro_bndGSI), intent(in)    :: self
  type(ufo_geovals),        intent(in)    :: geovals
  real(kind_real),          intent(inout) :: hofx(:)
  type(c_ptr), value,       intent(in)    :: obss
  integer                                 ::  nlocs
  character(len=*), parameter     :: myname  = "ufo_gnssro_bndgsi_simobs"
  character(max_string)           :: err_msg
  integer, parameter              :: nlevAdd = 13 !num of additional levels on top of exsiting model levels
  integer, parameter              :: newAdd  = 20 !num of additional levels on top of extended levels
  integer, parameter              :: ngrd    = 80 !num of new veritcal grids for bending angle computation
  integer                         :: iobs, k, igrd
  integer                         :: nlev, nlev1, nlevExt, nlevCheck
  real(kind_real)                 :: rnlevExt
  real(kind_real)                 :: w4(4), dw4(4)
  type(ufo_geoval), pointer       :: t, q, gph, prs
  real(kind_real), allocatable    :: gesT(:,:), gesZ(:,:), gesP(:,:), gesQ(:,:), gesTv(:,:)
  real(kind_real), allocatable    :: ref(:), radius(:)
  real(kind_real), allocatable    :: refIndex(:), refXrad(:), geomz(:), refXrad_new(:)
  real(kind_real), allocatable    :: lagConst(:,:)
  real(kind_real), allocatable    :: obsLat(:), obsImpP(:), obsLocR(:), obsGeoid(:), obsValue(:)
  real(kind_real), allocatable    :: temperature(:)
  real(kind_real)                 :: wf
  integer                         :: wi, wi2
  real(kind_real)                 :: d_refXrad
  real(kind_real)                 :: derivRef_s(ngrd),grids(ngrd),refXrad_s(ngrd)
  real(kind_real)                 :: sIndx, sIndxExt
  integer                         :: indx
  real(kind_real)                 :: bndIntgd, bendingAngle, obsImpH, gradRef
  logical                         :: obs_check, qc_layer_SR
  integer                         :: nobs_outIntgl
  integer                         :: count_SR, top_layer_SP, top_layer_SR, bot_layer_SR !for super refaction
  integer                         :: count_rejection
  integer                         :: iflip
  real(c_double)                  :: missing

  write(err_msg,*) myname, ": begin"
  call fckit_log%info(err_msg)

! check if nobs is consistent in geovals & hofx
  if (geovals%nlocs /= size(hofx)) then
    write(err_msg,*) myname, ': nlocs inconsistent!'
    call abor1_ftn(err_msg)
  endif

  missing = missing_value(missing)
  nlocs   = obsspace_get_nlocs(obss) ! number of observations

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
  allocate(ObsValue(nlocs))

  call obsspace_get_db(obss, "MetaData", "latitude",         obsLat)
  call obsspace_get_db(obss, "MetaData", "impact_parameter", obsImpP)
  call obsspace_get_db(obss, "MetaData", "earth_radius_of_curvature", obsLocR)
  call obsspace_get_db(obss, "MetaData", "geoid_height_above_reference_ellipsoid", obsGeoid)
  call obsspace_get_db(obss, "ObsValue", "bending_angle", obsValue)

  sIndxExt  = one
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
  allocate(lagConst(3,nlevExt))         !
  allocate(refXrad_new(nlevExt+newAdd)) !

  allocate(temperature(nlocs))

  nobs_outIntgl = 0 !initialize count of observations out of integral grids  
  count_rejection = 0


  obs_loop: do iobs = 1, nlocs

    hofx(iobs) =  missing
    temperature(iobs) =  missing

    do k = 1, nlev
!     compute guess geometric height from geopotential height
      call geop2geometric(obsLat(iobs), gesZ(k,iobs), geomz(k))
      radius(k) = geomz(k) + obsGeoid(iobs) + obsLocR(iobs)   ! radius r
!     guess refactivity, refactivity index,  and impact parameter
      call compute_refractivity(gesT(k,iobs), gesQ(k,iobs), gesP(k,iobs),   &
                                ref(k), self%roconf%use_compress)
      refIndex(k) = one + (r1em6*ref(k))
      refXrad(k)  = refIndex(k) * radius(k)
    end do 

!   data rejection based on model background !
!   (1) skip data beyond model levels
    call get_coordinate_value(obsImpP(iobs),sIndx,refXrad(1),nlev,"increasing")
    if (sIndx < one .or. sIndx > float(nlev))    cycle obs_loop

!   calculating temeprature at obs location to obs space for BackgroundCheck ROGSI
    indx=sIndx
    wi=min(max(1,indx),nlev)
    wi2=max(1,min(indx+1,nlev))
    wf=indx-float(wi)
    wf=max(zero,min(wf,one))
    temperature(iobs)=gesT(wi,iobs)*(one-wf)+gesT(wi2,iobs)*wf

!   (2) super-refaction
    qc_layer_SR  = .false.
    count_SR     = 0
    top_layer_SR = 0
    bot_layer_SR = 0

    obsImpH = (obsImpP(iobs) - obsLocR(iobs)) * r1em3 !impact heigt: a-r_earth
    if (obsImpH <= six) then
       do k = nlevCheck, 1, -1

!         check for model SR layer
          gradRef = 1000.0_kind_real * (ref(k+1)-ref(k)) /       &
                                    (radius(k+1)-radius(k))
!         PLEASE KEEP this COMMENT:
!         this check needs RO profile, which was done with MPI reduce in GSI
!         not applied here yet
!         only check once - SR-likely layer detected
          if (.not.qc_layer_SR .and. abs(gradRef)>= half*crit_gradRefr) then
             qc_layer_SR=.true. 
          endif

!         relax to close-to-SR conditions
          if (abs(gradRef) >= 0.75_kind_real*crit_gradRefr) then
             count_SR=count_SR+1        ! layers of SR
             if (count_SR > 1 ) then
                bot_layer_SR=k
             else
                top_layer_SR=k
                bot_layer_SR=top_layer_SR
             endif
          endif
       end do
     
!      obs inside model SR layer
       if (top_layer_SR >= 1 .and. obsImpP(iobs) <= refXrad(top_layer_SR+2)) then
          cycle obs_loop
       end if

    end if ! obsImpH <= six

!  Extend atmosphere above interface level nlev
    d_refXrad = refXrad(nlev) - refXrad(nlev-1)
    do k = 1, nlevAdd
      refXrad(nlev+k)=refXrad(nlev)+ k*d_refXrad    ! extended x_i
      ref(nlev+k)=ref(nlev+k-1)**2/ref(nlev+k-2) ! exended N_i
    end do

    refXrad(0)=refXrad(3)
    refXrad(nlevExt+1)=refXrad(nlevExt-2)
    do k = 1,nlevExt
       call lag_interp_const(lagConst(:,k),refXrad(k-1:k+1),3)
    enddo

!   integrate on a new set of equally-spaced vertical grid 
    derivRef_s = zero
    grids_loop: do igrd =1,ngrd
      refXrad_s(igrd)=sqrt(grids(igrd)**2 + obsImpP(iobs)**2) !x_s^2=s^2+a^2
      call get_coordinate_value(refXrad_s(igrd), sIndx,refXrad(1:nlevExt),nlevExt,"increasing")
      rnlevExt = float(nlevExt)

      if (sIndx > zero .and. sIndx < rnlevExt) then  !obs inside the new grid
        indx=sIndx
!       Compute derivative at new grids (dN/ds) using Lagrange interpolators
        call lag_interp_smthWeights(refXrad(indx-1:indx+2),refXrad_s(igrd),&
                        lagConst(:,indx),lagConst(:,indx+1),&
                        w4,dw4,4)
        if (indx==1) then
          w4(4)=w4(4)+w4(1); w4(1:3)=w4(2:4);w4(4)=zero
          dw4(4)=dw4(4)+dw4(1);dw4(1:3)=dw4(2:4);dw4(4)=zero
          indx=indx+1
        endif
        if (indx==nlevExt-1) then
          w4(1)=w4(1)+w4(4); w4(2:4)=w4(1:3);w4(1)=zero
          dw4(1)=dw4(1)+dw4(4); dw4(2:4)=dw4(1:3);dw4(1)=zero
          indx=indx-1
        endif

        derivRef_s(igrd)=dot_product(dw4,ref(indx-1:indx+2)) !derivative dN/dx_s
        derivRef_s(igrd)=max(zero,abs(derivRef_s(igrd)))

      else
        cycle  obs_loop
      endif !obs in new grid
    end do grids_loop

!   bending angle (radians)
    bendingAngle = ds*derivRef_s(1)/refXrad_s(1)
    do igrd = 2,ngrd
       bndIntgd     = ds*derivRef_s(igrd)/refXrad_s(igrd)
       bendingAngle = bendingAngle + two*bndIntgd
    end do
    bendingAngle=r1em6 * obsImpP(iobs) * bendingAngle
    hofx(iobs) = bendingAngle

  end do obs_loop

! putting temeprature at obs location to obs space tor BackgroundCheck ROGSI
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
  deallocate(ref) 
  deallocate(refIndex) 
  deallocate(refXrad) 
  deallocate(geomz) 
  deallocate(radius) 
  deallocate(lagConst) 
  deallocate(refXrad_new) 
  deallocate(temperature)

  write(err_msg,*) myname, ": complete"
  call fckit_log%info(err_msg)
end subroutine ufo_gnssro_bndgsi_simobs
! ------------------------------------------------------------------------------
end module ufo_gnssro_bndgsi_mod
