real function interp_weight(d, x, nx)
implicit none

integer, intent(in) :: nx
real, intent(in)    :: x(nx)
real, intent(in)    :: d

integer ::  ix


! Case in which x is in decreasing order
if(d>=x(1)) then
  ix = 1
else
  ix = 1
  do while (d < x(ix))
    ix = ix + 1
    if (ix == nx) exit
  enddo
  ix = ix - 1
  interp_weight = real(ix) + (d-x(ix)) / (x(ix+1)-x(ix))

end function interp_weight 

