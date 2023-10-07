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
  use gnssro_mod_grids,   only: get_coordinate_value
  use ufo_constants_mod,  only: two, zero

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
           bendingAngle,super_refraction_flag,top_layer_SR,super_ref_GEOS)
! -------------------------------------------------------------------------------
  character(len=*), parameter    :: myname  = "ufo_gnssro_bndnbam_simobs_single"

  real(kind_real), intent(out)   :: bendingAngle
 
  integer, intent(in)            :: nlev
  integer, intent(in)            :: nlevExt
  integer, intent(in)            :: nlevAdd
  integer, intent(in)            :: nlevCheck
  integer, intent(in)            :: ngrd
  integer, intent(in)            :: top_layer_SR
  logical, intent(in)            :: super_ref_GEOS
  real(kind_real), intent(in)    :: obsLat, obsGeoid, obsLocR, obsImpP ! obsspace
  real(kind_real), intent(in)    :: grids(ngrd)
 
  real(kind_real), intent(in)    :: radius(nlev)
  real(kind_real), intent(in)    :: refIndex(nlev)
  real(kind_real), intent(inout) :: ref(nlevExt)
  real(kind_real), intent(inout) :: refXrad(0:nlevExt+1)
  real(kind_real)                :: lagConst(3,nlevExt)

  real(kind_real)                 :: sIndx
  real(kind_real)                 :: obsImpH
  integer                         :: k, igrd, indx
  real(kind_real)                 :: d_refXrad
  real(kind_real)                 :: w4(4), dw4(4)
  real(kind_real)                 :: bndIntgd
  real(kind_real)                 :: rnlevExt
  real(kind_real)                 :: derivRef_s(ngrd)
  real(kind_real)                 :: refXrad_s(ngrd)
  integer, intent(inout) :: super_refraction_flag
!------------------------------------------------------------
! Extend atmosphere above interface level nlev
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

! integrate on a new set of equally-spaced vertical grid 
  derivRef_s = zero
  grids_loop: do igrd =1,ngrd
    refXrad_s(igrd)=sqrt(grids(igrd)**2 + obsImpP**2) !x_s^2=s^2+a^2
    call get_coordinate_value(refXrad_s(igrd), sIndx,refXrad(1:nlevExt),nlevExt,"increasing")
    if (top_layer_SR > 0 .and. sIndx < float(top_layer_SR+1)) then
       call get_coordinate_value(refXrad_s(igrd), sIndx,refXrad((top_layer_SR+1):nlevExt),&
            nlevExt-top_layer_SR-1,"increasing")
       sIndx = sIndx+top_layer_SR
    endif
    rnlevExt = float(nlevExt)

    if (sIndx > zero .and. sIndx < rnlevExt) then  !obs inside the new grid
       indx=sIndx
!      Compute derivative at new grids (dN/ds) using Lagrange interpolators
       call lag_interp_smthWeights(refXrad(indx-1:indx+2),refXrad_s(igrd),&
                                   lagConst(:,indx),lagConst(:,indx+1),   &
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
       if (super_ref_GEOS) then
          derivRef_s(igrd)=max(zero,abs(derivRef_s(igrd)))
       else
          if (derivRef_s(igrd).gt.zero) then
             super_refraction_flag = 3
             return
          end if
       end if

    else
      return
    endif !obs in new grid
  end do grids_loop

! bending angle (radians)
  bendingAngle = ds*derivRef_s(1)/refXrad_s(1)
  do igrd = 2,ngrd
     bndIntgd     = ds*derivRef_s(igrd)/refXrad_s(igrd)
     bendingAngle = bendingAngle + two*bndIntgd
  end do
  if (super_ref_GEOS) then
     bendingAngle=r1em6 * obsImpP * bendingAngle
  else
     bendingAngle = (-1)* r1em6 * obsImpP * bendingAngle
  end if

end subroutine ufo_gnssro_bndnbam_simobs_single
! ------------------------------------------------------------------------------
end module ufo_gnssro_bndnbam_util_mod
