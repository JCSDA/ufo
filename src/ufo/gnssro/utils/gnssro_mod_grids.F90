module gnssro_mod_grids

use kinds, only: kind_real
use ufo_constants_mod, only: one

public :: get_coordinate_value

private

contains

subroutine get_coordinate_value(fin, fout, x, nx, flag)
!
! Get grid coordinates from monotonically increasing or decreasing points
! adapted GSI subprogram:    grdcrd1 
!
  integer,         intent(in)   :: nx    !number of reference grid point
  real(kind_real), intent(in)   :: x(nx) !grid values
  real(kind_real), intent(in)   :: fin   !input point
  character(10),   intent(in)   :: flag  !"increasing" or "decreasing"
  real(kind_real), intent(out)  :: fout  !output point
  integer                       :: ix, isrchf

! Treat "normal" case in which nx>1
  if(nx>1) then
     if (flag == "increasing") then

        if(fin<=x(1)) then
           ix=1
        else
           call searchArray(nx-1,x,fin,flag,isrchf)
           ix=isrchf-1
        end if
        if(ix==nx) ix=ix-1

     else if (flag=="decreasing") then

        if(fin>=x(1)) then
           ix=1
        else
           call searchArray(nx-1,x,fin,flag,isrchf)
           ix=isrchf-1
        end if
     else
        ix = 1
        call abor1_ftn('gnssro get_coordinate_value: flag must be set to "decreasing" or "increasing"')
     end if
     fout=float(ix)+(fin-x(ix))/(x(ix+1)-x(ix))

! Treat special case of nx=1
  elseif (nx==1) then
     fout = one
  endif

  return
end subroutine get_coordinate_value


subroutine searchArray(nx,x,y,flag,isrchf)
  integer,        intent(in)  :: nx     !number of input points
  character(10),  intent(in)  :: flag   !"increasing" or "decreasing"
  real(kind_real),intent(in)  :: y      !target values
  real(kind_real),intent(in)  :: x(nx)  !grid value
  integer,        intent(out) :: isrchf !array index of input grid value near target value
  integer                     :: k

  if(flag=="increasing") then
     do k=1,nx
        if(y<=x(k)) then
           isrchf=k
           return
        end if
     end do
  else
     do k=1,nx
        if(y>=x(k)) then
           isrchf=k
           return
        end if
     end do
  end if

  isrchf=nx+1
  if(nx<=0) isrchf=0

  return
end subroutine searchArray

end module gnssro_mod_grids
