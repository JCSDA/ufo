!==========================================================================
module gnssro_mod_obserror
!==========================================================================

use kinds
use gnssro_mod_constants

contains

subroutine obs_error(obsLat, obsZ, nobs, GlobalModel, ERR_TYPE, obsErr)
integer,                         intent(in)  :: nobs, GlobalModel
real(kind_real), dimension(nobs),intent(in)  :: obsLat,  obsZ
real(kind_real), dimension(nobs),intent(out) :: obsErr
character(len=255),intent(in)  :: ERR_TYPE
real(kind_real)                :: obsZ_km

select case (trim(ERR_TYPE))

case ("RefGSI")
call refractivity_err_gsi(obsLat, obsZ, nobs, GlobalModel, obsErr)

case ("BndGSI")

end select

end subroutine obs_error

subroutine  refractivity_err_gsi(obsLat, obsZ, nobs,GlobalModel, obsErr)
integer,                         intent(in)  :: nobs, GlobalModel
real(kind_real), dimension(nobs),intent(in)  :: obsLat,obsZ
real(kind_real), dimension(nobs),intent(out) :: obsErr
real(kind_real)                :: obsZ_km

do i = 1, nobs
obsZ_km  = obsZ(nobs)/1000.0_kind_real

if( GlobalModel .eq. 1 ) then ! for global

     if( obsLat(i)>= 20.0 .or.obsLat(i)<= -20.0 ) then
         obsErr(i)=-1.321_kind_real+0.341_kind_real*obsZ_km-0.005_kind_real*obsZ_km**2
     else
       if(obsZ_km > 10.0) then
          obsErr(i)=2.013_kind_real-0.060_kind_real*obsZ_km+0.0045_kind_real*obsZ_km**2
       else
          obsErr(i)=-1.18_kind_real+0.058_kind_real*obsZ_km+0.025_kind_real*obsZ_km**2
       endif
     endif
     obsErr(i) = 1.0_kind_real/abs(exp(obsErr(i)))

else ! for regional 
     if( obsLat(i) >= 20.0 .or.obsLat(i) <= -20.0 ) then
         if (obsZ_km > 10.00) then
             obsErr(i) =-1.321_kind_real+0.341_kind_real*obsZ_km-0.005_kind_real*obsZ_km**2
         else
             obsErr(i) =-1.2_kind_real+0.065_kind_real*obsZ_km+0.021_kind_real*obsZ_km**2
         endif
     else
         if(obsZ_km > 10.00) then
            obsErr(i) =2.013_kind_real-0.120_kind_real*obsZ_km+0.0065_kind_real*obsZ_km**2
         else
            obsErr(i) =-1.19_kind_real+0.03_kind_real*obsZ_km+0.023_kind_real*obsZ_km**2
         endif
     endif
     obsErr(i) = 1.0_kind_real/abs(exp(obsErr(i)))
endif
end do
end subroutine refractivity_err_gsi

end module gnssro_mod_obserror
