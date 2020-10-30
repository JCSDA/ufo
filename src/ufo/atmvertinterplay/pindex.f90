real function pindex(nx, press, obspressure)
! This routine handles press (pressure) in decendent order (bottom2top)
! press(1): highest pressure value
! press(nx): lowerest pressure value

use kinds
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
end
