real function pindex(nsig, press, obspressure)

use kinds
implicit none

integer :: ix, k, nsig
real(kind_real) :: ozp, obspressure, psi
real(kind_real), dimension(nsig) :: press

psi = 1./press(1)
if(obspressure*psi < 1.) then
  ozp = obspressure
else
  ozp = press(1)
endif
if( ozp >= press(1)) then
  ix = 1
else
  ix = 0
  do k = 1, nsig
    if(ozp >= press(k)) then
      ix = k
      exit
    endif
  enddo
  if(ix == 0) ix = nsig + 1
  if(ix > 1)ix = ix -1
endif
ozp = float(ix) + &
    (ozp-press(ix))/(press(ix+1)-press(ix))
pindex = ozp
return
end
