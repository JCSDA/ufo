!==========================================================================
module gnssro_mod_obserror
!==========================================================================

use kinds
use gnssro_mod_constants

contains
subroutine bending_angle_obserr_ROPP(obsImpH, obsValue, nobs,  obsErr, QCflags, missing)
implicit none
integer,                         intent(in)  :: nobs
real(kind_real), dimension(nobs),intent(in)  :: obsImpH, obsValue
integer(c_int),  dimension(nobs),intent(in)  :: QCflags(:)
real(kind_real), dimension(nobs),intent(out) :: obsErr
real(kind_real)                  :: H_km, missing
integer :: i

obsErr = missing

do i = 1, nobs
   
if (QCflags(i) .eq. 0) then

    H_km  = obsImpH(i)/1000.0_kind_real

    if ( H_km <= 10.0 ) then
       obsErr(i) = (H_km*1.25 + (10-H_km)*20)/10.0
       obsErr(i) = obsErr(i)/100.0*obsValue(i)
    else if ( H_km > 10.0 .and. H_km <= 32.0 ) then
      obsErr(i) = 1.25/100.0*obsValue(i) 
    else
      obsErr(i) = 3.0*1e-6
    end if
end if

end do

end subroutine bending_angle_obserr_ROPP
!---------------------------------------

subroutine  bending_angle_obserr_NBAM(obsLat, obsImpH, obsSaid, nobs, obsErr, QCflags, missing)
implicit none
integer,                         intent(in)  :: nobs
real(kind_real), dimension(nobs),intent(in)  :: obsImpH, obsLat
integer(c_int),  dimension(nobs),intent(in)  :: obsSaid, QCflags(:)
real(kind_real), dimension(nobs),intent(out) ::  obsErr
real(kind_real)                 :: H_km, missing

integer :: i

obsErr = missing

do i = 1, nobs

if (QCflags(i) .eq. 0) then

   H_km  = obsImpH(i)/1000.0_kind_real
   if( (ObsSaid(i)==41).or.(ObsSaid(i)==722).or.(ObsSaid(i)==723).or.   &
       (ObsSaid(i)==4).or.(ObsSaid(i)==42).or.(ObsSaid(i)==3).or.       &
       (ObsSaid(i)==5).or.(ObsSaid(i)==821.or.(ObsSaid(i)==421)).or.    &
       (ObsSaid(i)==440).or.(ObsSaid(i)==43)) then
       if( abs(obsLat(i))>= 40.00 ) then
         if(H_km>12.0) then
           obsErr(i)=0.19032 +0.287535 *H_km-0.00260813*H_km**2
         else
           obsErr(i)=-3.20978 +1.26964 *H_km-0.0622538 *H_km**2
         endif
       else
         if(H_km>18.) then
           obsErr(i)=-1.87788 +0.354718 *H_km-0.00313189 *H_km**2
         else
           obsErr(i)=-2.41024 +0.806594 *H_km-0.027257 *H_km**2
         endif
       endif

   else !!!! CDAAC processing
     if( abs(obsLat(i))>= 40.00 ) then
       if ( H_km>12.00 ) then
          obsErr(i)=-0.685627 +0.377174 *H_km-0.00421934 *H_km**2
       else
          obsErr(i)=-3.27737 +1.20003 *H_km-0.0558024 *H_km**2
       endif
     else
       if( H_km>18.00 ) then
          obsErr(i)=-2.73867 +0.447663 *H_km-0.00475603 *H_km**2
       else
          obsErr(i)=-3.45303 +0.908216 *H_km-0.0293331 *H_km**2
       endif
     endif
     obsErr(i) = 0.001 /abs(exp(obsErr(i)))

  endif

end if

end do

end subroutine bending_angle_obserr_NBAM
!---------------------------------------

subroutine refractivity_obserr_NBAM(obsLat, obsZ, nobs, obsErr, QCflags,missing)
implicit none
integer,                         intent(in)  :: nobs
real(kind_real), dimension(nobs),intent(in)  :: obsLat, obsZ
real(kind_real), dimension(nobs),intent(out) :: obsErr
integer(c_int),  dimension(nobs),intent(in)  :: QCflags(:)
real(kind_real)                   :: H_km, missing

integer :: i

obsErr = missing

do i = 1, nobs

  if (QCflags(i) .eq. 0) then
     H_km  = obsZ(i)/1000.0_kind_real
     if( abs(obsLat(i))>= 20.0 ) then
         obsErr(i)=-1.321+0.341*H_km-0.005*H_km**2
     else
       if(H_km > 10.0) then
          obsErr(i)=2.013-0.060*H_km+0.0045*H_km**2
       else
          obsErr(i)=-1.18+0.058*H_km+0.025*H_km**2
       endif
     endif
     obsErr(i) = 1.0_kind_real/abs(exp(obsErr(i)))
  end if
end do

end subroutine refractivity_obserr_NBAM

end module gnssro_mod_obserror

