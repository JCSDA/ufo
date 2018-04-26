! (C) Copyright 2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Fortran module to perform linear interpolation

module vert_interp_mod

use kinds, only: kind_real

implicit none
public

contains

! ------------------------------------------------------------------------------

subroutine vert_interp_weights(nlev,obl,vec,wi,wf)

implicit none
integer,         intent(in ) :: nlev       !Number of model levels
real(kind_real), intent(in ) :: obl        !Observation location
real(kind_real), intent(in ) :: vec(nlev)  !Structured vector of grid points
integer,         intent(out) :: wi         !Index for interpolation
real(kind_real), intent(out) :: wf         !Weight for interpolation

integer :: k

if (vec(1) < vec(nlev)) then !Pressure increases with index

  if (obl < vec(1)) then
     wi = 1
  elseif (obl > vec(nlev)) then
     wi = nlev - 1
  else
     do k = 1,nlev-1
        if (obl >= vec(k) .and. obl <= vec(k+1)) then
           wi = k
        endif
     enddo
  endif

else !Pressure decreases with index

  if (obl > vec(1)) then
     wi = 1
  elseif (obl < vec(nlev)) then
     wi = nlev - 1
  else
     do k = 1,nlev-1
        if (obl >= vec(k+1) .and. obl <= vec(k)) then
           wi = k
        endif
     enddo
  endif

endif

wf = (vec(wi+1) - obl)/(vec(wi+1) - vec(wi))

end subroutine vert_interp_weights

! ------------------------------------------------------------------------------

subroutine vert_interp_apply(nlev, fvec, f, wi, wf) 

implicit none
integer,         intent(in ) :: nlev        !Number of model levels
real(kind_real), intent(in ) :: fvec(nlev)  !Field at grid points
integer,         intent(in ) :: wi          !Index for interpolation
real(kind_real), intent(in ) :: wf          !Weight for interpolation
real(kind_real), intent(out) :: f           !Output at obs location using linear interp

f = fvec(wi)*wf + fvec(wi+1)*(1.0-wf)

end subroutine vert_interp_apply

! ------------------------------------------------------------------------------

subroutine vert_interp_apply_tl(nlev, fvec_tl, f_tl, wi, wf) 

implicit none
integer,         intent(in)  :: nlev
real(kind_real), intent(in)  :: fvec_tl(nlev)
integer,         intent(in)  :: wi
real(kind_real), intent(in)  :: wf
real(kind_real), intent(out) :: f_tl

f_tl = fvec_tl(wi)*wf + fvec_tl(wi+1)*(1.0_kind_real-wf)

end subroutine vert_interp_apply_tl

! ------------------------------------------------------------------------------

subroutine vert_interp_apply_ad(nlev, fvec_ad, f_ad, wi, wf) 

implicit none
integer,         intent(in)    :: nlev
real(kind_real), intent(inout) :: fvec_ad(nlev)
integer,         intent(in)    :: wi
real(kind_real), intent(in)    :: wf
real(kind_real), intent(inout) :: f_ad

fvec_ad(wi  ) = fvec_ad(wi  ) + f_ad*wf
fvec_ad(wi+1) = fvec_ad(wi+1) + f_ad*(1.0_kind_real-wf)

end subroutine vert_interp_apply_ad

! ------------------------------------------------------------------------------

end module vert_interp_mod
