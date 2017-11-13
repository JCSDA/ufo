real function vert_interp(f, nsig, dz)
implicit none

integer :: nsig
real, intent(in)  :: f(nsig)
real, intent(in)  :: dz

integer :: iz, izp
real    :: delz, delzp

iz=int(dz)
iz=max(1,min(iz,nsig))  
izp=min(iz+1,nsig)

delz=dz-float(iz)
delz=max(0.,min(delz,1.))
delzp=1.-delz

vert_interp = f(iz)*delzp + f(izp)*delz

end function vert_interp
