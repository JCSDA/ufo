!  
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Fortran module of Gnssro NBAM (NCEP's Bending Angle Method) operator

module ufo_gnssro_bndnbam_util_mod
  use fckit_log_module, only: fckit_log
  use kinds
  use missing_values_mod
  use gnssro_mod_constants
  use vert_interp_mod
  use lag_interp_mod,     only: lag_interp_const, lag_interp_smthWeights
  use gnssro_mod_transform
  use gnssro_mod_grids, only: get_coordinate_value

  implicit none
  public             :: ufo_gnssro_bndnbam_simobs_single
  private

  contains
! ------------------------------------------------------------------------------
subroutine ufo_gnssro_bndnbam_simobs_single( &
           obsLat, obsGeoid, obsLocR, obsImpP, &
           grids, ngrd, &
           nlev, nlevExt, nlevAdd, nlevCheck, &
           radius,ref,refIndex,refXrad, &
           bendingAngle)
! -------------------------------------------------------------------------------
  character(len=*), parameter    :: myname  = "ufo_gnssro_bndnbam_simobs_single"

  real(kind_real), intent(out)   :: bendingAngle
 
  integer, intent(in)            :: nlev
  integer, intent(in)            :: nlev1
  integer, intent(in)            :: nlevExt
  integer, intent(in)            :: nlevAdd
  integer, intent(in)            :: nlevCheck
  integer, intent(in)            :: ngrd
  integer, intent(in)            :: use_compress
  real(kind_real), intent(in)    :: obsLat
  real(kind_real), intent(in)    :: obsGeoid
  real(kind_real), intent(in)    :: obsLocR
  real(kind_real), intent(in)    :: obsImpP
 
  real(kind_real), intent(in)    :: gesZ(nlev1)
  real(kind_real), intent(in)    :: gesT(nlev)
  real(kind_real), intent(in)    :: gesQ(nlev)
  real(kind_real), intent(in)    :: gesP(nlev1)
  real(kind_real), intent(in)    :: grids(ngrd)
 
  real(kind_real)                :: radius(nlev)
  real(kind_real)                :: ref(nlevExt)
  real(kind_real)                :: refIndex(nlev)
  real(kind_real)                :: refXrad(0:nlevExt+1)
  real(kind_real)                :: lagConst(3,nlevExt)
 
  real(kind_real)                 :: sIndx
  real(kind_real)                 :: obsImpH
  real(kind_real)                 :: geomz
  integer                         :: klev
  integer                         :: igrd
  integer                         :: indx
  integer                         :: wi
  integer                         :: wi2
  real(kind_real)                 :: wf
  real(kind_real)                 :: temperature
  logical                         :: qc_layer_SR
  integer                         :: count_SR
  integer                         :: top_layer_SR
  integer                         :: bot_layer_SR 
  real(kind_real)                 :: gradRef
  real(kind_real)                 :: d_refXrad
  real(kind_real)                 :: w4(4), dw4(4)
  real(kind_real)                 :: bndIntgd
  real(kind_real)                 :: rnlevExt
  real(kind_real)                 :: derivRef_s(ngrd)
  real(kind_real)                 :: refXrad_s(ngrd)

!------------------------------------------------------------

  do klev = 1, nlev
!    compute guess geometric height from geopotential height
     call geop2geometric(obsLat, gesZ(klev), geomz)
     radius(klev) = geomz + obsGeoid + obsLocR   ! radius r

!    guess refactivity, refactivity index,  and impact parameter
     call compute_refractivity(gesT(klev), gesQ(klev), gesP(klev),   &
                                ref(klev), use_compress)
     refIndex(klev) = one + (r1em6*ref(klev))
     refXrad(klev)  = refIndex(klev) * radius(klev)
  end do 

! data rejection based on model background !
!  (1) skip data beyond model levels
   call get_coordinate_value(obsImpP,sIndx,refXrad(1),nlev,"increasing")
   if (sIndx < one .or. sIndx > float(nlev))  then
      return
   endif

!   calculating temeprature at obs location to obs space for BackgroundCheck RONBAM
   indx=sIndx
   wi=min(max(1,indx),nlev)
   wi2=max(1,min(indx+1,nlev))
   wf=indx-float(wi)
   wf=max(zero,min(wf,one))
   temperature=gesT(wi)*(one-wf)+gesT(wi2)*wf

!  (2) super-refaction
   qc_layer_SR  = .false.
   count_SR     = 0
   top_layer_SR = 0
   bot_layer_SR = 0

   obsImpH = (obsImpP - obsLocR) * r1em3 !impact heigt: a-r_earth
   if (obsImpH <= six) then
      do klev = nlevCheck, 1, -1

!         check for model SR layer
          gradRef = 1000.0_kind_real * (ref(klev+1)-ref(klev)) /       &
                                    (radius(klev+1)-radius(klev))
         if (.not.qc_layer_SR .and. abs(gradRef)>= half*crit_gradRefr) then
            qc_layer_SR=.true. 
         endif

!         relax to close-to-SR conditions
         if (abs(gradRef) >= 0.75_kind_real*crit_gradRefr) then
            count_SR=count_SR+1        ! layers of SR
            if (count_SR > 1 ) then
               bot_layer_SR=klev
            else
               top_layer_SR=klev
               bot_layer_SR=top_layer_SR
            endif
        endif
      end do
     
!      obs inside model SR layer
       if (top_layer_SR >= 1 .and. obsImpP <= refXrad(top_layer_SR+2)) then
          return
       end if

    end if ! obsImpH <= six

!  Extend atmosphere above interface level nlev
    d_refXrad = refXrad(nlev) - refXrad(nlev-1)
    do klev = 1, nlevAdd
      refXrad(nlev+klev)=refXrad(nlev)+ klev*d_refXrad    ! extended x_i
      ref(nlev+klev)=ref(nlev+klev-1)**2/ref(nlev+klev-2) ! exended N_i
    end do

    refXrad(0)=refXrad(3)
    refXrad(nlevExt+1)=refXrad(nlevExt-2)
    do klev = 1,nlevExt
       call lag_interp_const(lagConst(:,klev),refXrad(klev-1:klev+1),3)
    enddo

!   integrate on a new set of equally-spaced vertical grid 
    derivRef_s = zero
    grids_loop: do igrd =1,ngrd
      refXrad_s(igrd)=sqrt(grids(igrd)**2 + obsImpP**2) !x_s^2=s^2+a^2
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
        return
      endif !obs in new grid
    end do grids_loop

!   bending angle (radians)
    bendingAngle = ds*derivRef_s(1)/refXrad_s(1)
    do igrd = 2,ngrd
       bndIntgd     = ds*derivRef_s(igrd)/refXrad_s(igrd)
       bendingAngle = bendingAngle + two*bndIntgd
    end do
    bendingAngle=r1em6 * obsImpP * bendingAngle

end subroutine ufo_gnssro_bndnbam_simobs_single
! ------------------------------------------------------------------------------
end module ufo_gnssro_bndnbam_util_mod
