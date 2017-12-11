subroutine add_rtm_layers(prsitmp,prsltmp,prsitmp_ext,prsltmp_ext,nsig,msig,klevel)
! !USES:

use crtm_module, only: toa_pressure

implicit none

! !INPUT PARAMETERS:
integer, intent(in) :: nsig, msig
integer, dimension(msig)  ,intent(  out) :: klevel

real(r_kind)   ,dimension(nsig+1),intent(in   ) :: prsitmp
real(r_kind)   ,dimension(nsig)  ,intent(in   ) :: prsltmp

real(r_kind)   ,dimension(msig+1),intent(  out) :: prsitmp_ext
real(r_kind)   ,dimension(msig)  ,intent(  out) :: prsltmp_ext

integer ::  k,kk,l
real(r_kind) dprs,toa_prs_kpa

integer, dimension(nsig) :: nlayers

nlayers = 1
nlayers(63) = 3
nlayers(64) = 6


!   Check if model top pressure above rtm top pressure, where prsitmp
!   is in kPa and toa_pressure is in hPa.
if (prsitmp(nsig) < toa_pressure)then
       write(6,*)'ADD_RTM_LAYERS:  model top pressure(hPa)=', &
            prsitmp(nsig),&
            ' above rtm top pressure(hPa)=',toa_pressure
       call stop2(35)
end if

!   Linear in pressure sub-divsions
kk=0
do k = 1,nsig
  if (nlayers(k)<=1) then
    kk = kk + 1
    prsltmp_ext(kk) = prsltmp(k)
    prsitmp_ext(kk) = prsitmp(k)
    klevel(kk) = k
  else
    if (k/=nsig) then
       dprs = (prsitmp(k+1)-prsitmp(k))/nlayers(k)
    else
       dprs = (toa_pressure-prsitmp(k))/nlayers(k)
    end if
    prsitmp_ext(kk+1) = prsitmp(k)
    do l=1,nlayers(k)
       kk=kk + 1
       prsitmp_ext(kk+1) = prsitmp(k) + dprs*l
       prsltmp_ext(kk) = half*(prsitmp_ext(kk+1)+prsitmp_ext(kk))
       klevel(kk) = k
    end do
  endif
end do

!   Set top of atmosphere pressure
prsitmp_ext(msig+1) = toa_pressure

end subroutine add_rtm_layers

