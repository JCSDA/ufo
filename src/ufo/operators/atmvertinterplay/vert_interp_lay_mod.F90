! (C) Copyright 2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Fortran module to perform linear interpolation

module vert_interp_lay_mod

use, intrinsic :: iso_c_binding
use kinds, only: kind_real
use missing_values_mod

implicit none
public

contains

subroutine get_integral_limits(airpressure, botpressure, toppressure, modelpressure, nlevs, nlocs, nsig) 
implicit none
integer,intent(in)  :: nlevs, nlocs, nsig
real(kind_real) , dimension(nlocs),intent(in) ::  airpressure
real(kind_real), dimension(nlocs),intent(inout)  ::  botpressure,toppressure
real(kind_real), dimension(nsig+1,nlocs),intent(in) :: modelpressure
! local
integer :: nprofs, iobs, iprof, kk, k1, k2
if (nlevs == 1) then ! total column ozone
  do iobs = 1, nlocs
    toppressure(iobs) = modelpressure(nsig+1,iobs)
    botpressure(iobs) = modelpressure(1,iobs)
  enddo
else
  !Obs pressures read in as Pa
  nprofs = nlocs/nlevs
  iobs = 0
  do iprof = 1, nprofs
    do kk = 1, nlevs
      k1 = kk
      k2 = kk - 1
      if (k2 == 0) k2 = 1
      if (kk == nlevs) then
        k1 = nlevs - 1
        k2 = 1
      endif
      iobs = iobs+1
      toppressure(iobs) = airpressure(k2)
      botpressure(iobs) = airpressure(k1)
      if( kk == 1 ) then
        toppressure(iobs) = modelpressure(nsig+1, iobs) 
        botpressure(iobs) = airpressure(k1)
        if(botpressure(iobs) < modelpressure(nsig+1, iobs)) then 
          botpressure(iobs) = modelpressure(nsig+1, iobs)
        endif
      else if( kk == nlevs) then
        toppressure(iobs) = modelpressure(nsig+1, iobs)  
        botpressure(iobs) = modelpressure(1, iobs) 
      endif
    enddo
  enddo
endif

end subroutine get_integral_limits

real function pindex(nx, press, obspressure)
! This routine handles press (pressure) in decendent order (bottom2top)
! press(1): highest pressure value
! press(nx): lowerest pressure value

implicit none

integer :: ix, k, nx
real(kind_real) :: ozp, obspressure, psi
real(kind_real), dimension(nx) :: press

psi = 1.0_kind_real/press(1)
if(obspressure*psi < 1.) then
  ozp = obspressure
else
  ozp = press(1)
endif
if( ozp >= press(1)) then
  ix = 1
else
  ix = 0
  do k = 1, nx-1
    if(ozp >= press(k)) then
      ix = k
      exit
    endif
  enddo
  if(ix == 0) ix = nx
  if(ix > 1)ix = ix -1
endif
ozp = float(ix) + &
    (ozp-press(ix))/(press(ix+1)-press(ix))
pindex = ozp
return
end function pindex

subroutine apply_layer_integral(coefficient, modelozone, modelpressure, botpressure, toppressure, nsig, layer_oz)
implicit none
integer,intent(in)  :: nsig
real, intent(in)  :: coefficient
real(kind_real),intent(in)  :: botpressure, toppressure
real(kind_real), dimension(:), intent(in) :: modelpressure, modelozone
real(kind_real), intent(out) :: layer_oz
! local
integer :: kk, iz1, iz2
real(kind_real) :: pob,delz,g,delp4,dz1
real(kind_real) :: topozp, botozp
topozp = pindex(nsig+1, modelpressure, toppressure)
botozp = pindex(nsig+1, modelpressure, botpressure)
pob = botozp
iz1 = topozp
if (iz1>nsig) iz1=nsig
iz2 = pob
layer_oz = 0._kind_real
dz1 = topozp
do kk=iz1,iz2,-1
  delz = 1.0_kind_real
  if(kk == iz1) delz = dz1 - iz1
  if (kk == iz2) delz = delz - pob + iz2
  delp4 = modelpressure(kk)-modelpressure(kk+1)  ! [Pa]
  layer_oz = layer_oz + modelozone(kk)*coefficient*(delz*delp4)
enddo

end subroutine apply_layer_integral

subroutine undo_layer_integral(coefficient, modelozone, modelpressure, botpressure, toppressure, nsig, layer_oz)
implicit none
integer,intent(in) :: nsig
real,intent(in) :: coefficient
real(kind_real),intent(in) :: botpressure, toppressure
real(kind_real), dimension(:),intent(in) :: modelpressure
real(kind_real), dimension(:),intent(out):: modelozone
real(kind_real), intent(in)  :: layer_oz
! local
integer :: kk, iz1, iz2
real(kind_real) :: pob,delz,g,delp4,dz1
real(kind_real) :: topozp, botozp
topozp = pindex(nsig+1, modelpressure, toppressure)
botozp = pindex(nsig+1, modelpressure, botpressure)
pob = botozp
iz1 = topozp
if (iz1>nsig) iz1=nsig
iz2 = pob
dz1 = topozp
modelozone = 0.0_kind_real
do kk=iz1,iz2,-1
  delz = 1.0_kind_real
  if(kk == iz1) delz = dz1 - iz1
  if (kk == iz2) delz = delz - pob + iz2
  delp4 = modelpressure(kk)-modelpressure(kk+1)  ! [Pa]
  modelozone(kk) = modelozone(kk) + layer_oz*coefficient*(delz*delp4)
enddo

end subroutine undo_layer_integral


subroutine vert_interp_lay_apply_tl(modelozoned, layer_ozd, coefficient,  modelpressure, botpressure, toppressure, nsig)
real(kind_real), intent(in)  ::modelozoned(:)
real(kind_real),intent(in):: botpressure, toppressure
real(kind_real), intent(in), dimension(nsig+1) :: modelpressure 
real(kind_real), intent(out) :: layer_ozd
integer,intent(in) :: nsig
real,intent(in) :: coefficient
call apply_layer_integral(coefficient, modelozoned, modelpressure, botpressure, toppressure, nsig, layer_ozd)
end subroutine vert_interp_lay_apply_tl 

subroutine vert_interp_lay_apply_ad(modelozoneb, layer_ozb, coefficient,  modelpressure, botpressure, toppressure, nsig)
real(kind_real), intent(out)  ::modelozoneb(:)
real(kind_real),intent(in):: botpressure, toppressure
real(kind_real), intent(in), dimension(nsig+1) :: modelpressure 
real(kind_real), intent(in) :: layer_ozb
integer,intent(in) :: nsig
real,intent(in) :: coefficient
call undo_layer_integral(coefficient, modelozoneb, modelpressure, botpressure, toppressure, nsig, layer_ozb)
end subroutine vert_interp_lay_apply_ad


end module vert_interp_lay_mod
