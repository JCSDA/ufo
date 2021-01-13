! (C) Copyright 2020 NOAA NWS NCEP EMC
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module of helper functions for FV3-LAM ESG grid domain configuration
!> These routines are borrowed from regional_esg_grid.fd from the
!> RRFS Regional Workflow / UFS_UTILS repository
!> Credit goes to R. Jim Purser for original source of these subroutines

module esg_grid_mod

  use kinds
  use ufo_constants_mod, only: pi, deg2rad, rad2deg, zero, one, two
  implicit none
  private
  public :: gtoxm_ak_dd, gtoxm_ak_rr

  interface gtoxm_ak_rr
     module procedure gtoxm_ak_rr_m,gtoxm_ak_rr_g;                end interface
  interface gtoxm_ak_dd
     module procedure gtoxm_ak_dd_g;                end interface
  interface grtoc
   module procedure dgrtoc
                                                                  end interface

  logical ,parameter:: T=.true.,F=.false. !<- for pain-relief in logical ops

contains

!=============================================================================
subroutine gtoxm_ak_rr_m(A,K,plat,plon,pazi,lat,lon,xm,ff)!      [gtoxm_ak_rr]
!=============================================================================
! Given the map specification (angles in radians), the grid spacing (in
! map-space units) and the sample lat-lon (in radian), return the the
! image in map space in a 2-vector in grid units. If the transformation
! is invalid, return a .true. failure flag.
!=============================================================================
implicit none
real(kind_real),             intent(in ):: a,k,plat,plon,pazi,lat,lon
real(kind_real),dimension(2),intent(out):: xm
logical,              intent(out):: ff
real(kind_real),dimension(3,3):: prot,azirot
real(kind_real)               :: clat,slat,clon,slon,cazi,sazi
real(kind_real),dimension(3)  :: xc
!=============================================================================
clat=cos(plat); slat=sin(plat)
clon=cos(plon); slon=sin(plon)
cazi=cos(pazi); sazi=sin(pazi)

azirot(:,1)=(/ cazi, sazi, zero/)
azirot(:,2)=(/-sazi, cazi, zero/)
azirot(:,3)=(/   zero,   zero, one/)

prot(:,1)=(/     -slon,       clon,    zero/)
prot(:,2)=(/-slat*clon, -slat*slon,  clat/)
prot(:,3)=(/ clat*clon,  clat*slon,  slat/)
prot=matmul(prot,azirot)

call grtoc(lat,lon,xc)
xc=matmul(transpose(prot),xc)
call xctoxm_ak(a,k,xc,xm,ff)
end subroutine gtoxm_ak_rr_m
!=============================================================================
subroutine gtoxm_ak_rr_g(A,K,plat,plon,pazi,delx,dely,lat,lon,&! [gtoxm_ak_rr]
     xm,ff)
!=============================================================================
! Given the map specification (angles in radians), the grid spacing (in
! map-space units) and the sample lat-lon (in radian), return the the
! image in map space in a 2-vector in grid units. If the transformation
! is invalid, return a .true. failure flag.
!=============================================================================
implicit none
real(kind_real),             intent(in ):: a,k,plat,plon,pazi,delx,dely,lat,lon
real(kind_real),dimension(2),intent(out):: xm
logical,              intent(out):: ff
!=============================================================================
call gtoxm_ak_rr_m(A,K,plat,plon,pazi,lat,lon,xm,ff); if(ff)return
xm(1)=xm(1)/delx; xm(2)=xm(2)/dely
end subroutine gtoxm_ak_rr_g

!=============================================================================
subroutine gtoxm_ak_dd_g(A,K,pdlat,pdlon,pdazi,delx,dely,&!      [gtoxm_ak_dd]
dlat,dlon,     xm,ff)
!=============================================================================
! Like gtoxm_ak_rr_g, except input angles are expressed in degrees
!=============================================================================
implicit none
real(kind_real),             intent(in ):: a,k,pdlat,pdlon,pdazi,delx,dely,dlat,dlon
real(kind_real),dimension(2),intent(out):: xm
logical,              intent(out):: ff
!-----------------------------------------------------------------------------
real(kind_real):: plat,plon,pazi,lat,lon
!=============================================================================
plat=pdlat*deg2rad ! Convert these angles from degrees to radians
plon=pdlon*deg2rad !
pazi=pdazi*deg2rad !
lat=dlat*deg2rad
lon=dlon*deg2rad
call gtoxm_ak_rr_g(A,K,plat,plon,pazi,delx,dely,lat,lon,xm,ff)
end subroutine gtoxm_ak_dd_g
!=============================================================================
subroutine xctoxm_ak(a,k,xc,xm,ff)!                                [xctoxm_ak]
!=============================================================================
! Inverse mapping of xmtoxc_ak. That is, go from given cartesian unit
! 3-vector, xc, to map coordinate 2-vector xm (or return a raised failure
! flag, FF, if the attempt fails).
!=============================================================================
implicit none
real(kind_real),             intent(in ):: a,k
real(kind_real),dimension(3),intent(in ):: xc
real(kind_real),dimension(2),intent(out):: xm
logical,              intent(out):: ff
!-----------------------------------------------------------------------------
real(kind_real),dimension(2):: xs,xt
!=============================================================================
ff=F
call xctoxs(xc,xs)
call xstoxt(k,xs,xt,ff); if(ff)return
call xttoxm(a,xt,xm,ff)
end subroutine xctoxm_ak
!=============================================================================
subroutine xctoxs(xc,xs)!                                             [xctoxs]
!=============================================================================
! Inverse of xstoxc. I.e., cartesians to stereographic
!=============================================================================
implicit none
real(kind_real),dimension(3),intent(in ):: xc
real(kind_real),dimension(2),intent(out):: xs
!-----------------------------------------------------------------------------
real(kind_real):: zp
!=============================================================================
zp=one+xc(3); xs=xc(1:2)/zp
end subroutine xctoxs
!=============================================================================
subroutine xstoxt(k,xs,xt,ff)!                                        [xstoxt]
!=============================================================================
! Inverse of xttoxs.
!=============================================================================
implicit none
real(kind_real),             intent(in ):: k
real(kind_real),dimension(2),intent(in ):: xs
real(kind_real),dimension(2),intent(out):: xt
logical,              intent(out):: ff
!-----------------------------------------------------------------------------
real(kind_real):: s,sc
!=============================================================================
s=k*(xs(1)*xs(1)+xs(2)*xs(2)); sc=one-s
ff=abs(s)>=one; if(ff)return
xt=two*xs/sc
end subroutine xstoxt
!=============================================================================
subroutine xttoxm(a,xt,xm,ff)!                                       [xttoxm]
!=============================================================================
! Inverse of xmtoxt
!=============================================================================
implicit none
real(kind_real),             intent(in ):: a
real(kind_real),dimension(2),intent(in ):: xt
real(kind_real),dimension(2),intent(out):: xm
logical              ,intent(out):: ff
!-----------------------------------------------------------------------------
integer:: i
!=============================================================================
do i=1,2; call zttozm(a,xt(i),xm(i),ff); if(ff)return; enddo
end subroutine xttoxm
!=============================================================================
subroutine zttozm(a,zt,zm,ff)!                                        [zttozm]
!=============================================================================
! Inverse of zmtozt
!=============================================================================
implicit none
real(kind_real),intent(in ):: a,zt
real(kind_real),intent(out):: zm
logical, intent(out):: ff
!-----------------------------------------------------------------------------
real(kind_real):: ra,razt
!=============================================================================
ff=F
if    (a>zero)then; ra=sqrt( a); razt=ra*zt; zm=atan (razt)/ra
elseif(a<zero)then; ra=sqrt(-a); razt=ra*zt; ff=abs(razt)>=one; if(ff)return
                                           zm=atanh(razt)/ra
else                                     ; zm=zt
endif
end subroutine zttozm
!=============================================================================
subroutine dgrtoc(rlat,rlon,xe)!                                       [grtoc]
!=============================================================================
implicit none
real(kind_real),             intent(IN ):: rlat,rlon
real(kind_real),dimension(3),intent(OUT):: xe
!-----------------------------------------------------------------------------
real(kind_real)                         :: sla,cla,slo,clo
!=============================================================================
sla=sin(rlat);  cla=cos(rlat)
slo=sin(rlon);  clo=cos(rlon)
xe(1)=cla*clo; xe(2)=cla*slo; xe(3)=sla
end subroutine dgrtoc
!=============================================================================


end module esg_grid_mod
