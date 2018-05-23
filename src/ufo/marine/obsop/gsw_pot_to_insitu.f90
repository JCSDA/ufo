MODULE gsw_pot_to_insitu
!===============================================================================
!
! MODULE:
!  gsw_pot_to_insitu
!
! PURPOSE:
!  This module provides additional subroutines needed to compute the in situ
!  temperature from a numerical model's potential temperature (ºC) and practical 
!  salinity (psu). It utilizes the TEOS GSW Fortran package, version 3.03,
!  available at http://www.teos-10.org
!
!
!-------------------------------------------------------------------------------
! Author:  Steve Penny
! Institution: University of Maryland, College Park; NWS/NOAA/NCEP
! Contact: Steve.Penny@noaa.gov
!===============================================================================

IMPLICIT NONE

PUBLIC :: t_from_pt, p_from_z    ! Computing in situ from potential temperature
PUBLIC :: sa_from_sp, pt_from_t  ! Computing potential temperature from in situ

PRIVATE

INTEGER, PARAMETER :: r14 = selected_real_kind(14,30)

CONTAINS

!===============================================================================
FUNCTION t_from_pt(pt_in,sp_in,p_in,lon_in,lat_in)
!===============================================================================
  !USE common, ONLY: r_size

  use gsw_mod_toolbox, only : gsw_ct_from_pt  
  use gsw_mod_toolbox, only : gsw_sa_from_sp
  use gsw_mod_toolbox, only : gsw_t_from_ct
  use gsw_mod_kinds

  !USE gsw_oceanographic_toolbox, ONLY: gsw_sa_from_sp, gsw_ct_from_pt, gsw_t_from_ct !, r14
  !STEVE: The toolbox is not coded as a module, so this will not work.
  !       I'll leave it unaltered as a list of functions in order to facilitate forward compatibility.

  IMPLICIT NONE
  REAL(r8) :: t_from_pt
  REAL(r8), INTENT(IN) :: pt_in, sp_in, p_in, lon_in, lat_in
  REAL(r8) :: pt, sp, p, lon, lat
  REAL(r8) :: sa, ct, T
 
  pt = pt_in
  sp = sp_in
  p  = p_in
  lon = lon_in
  lat = lat_in

! Compute in situ temperature from potential temperature by computing
! the conservative temperature from potential temperature then computing
! the in situ temperature from the conservative temperature.
! sa     : Absolute Salinity                   [g/kg]
! sp     : Practical Salinity                  [unitless]
! ct     : Conservative Temperature            [deg C]
! pt     : potential temperature with          [deg C]
!           reference pressure of 0 dbar
! c      : conductivity                        [mS/cm]
! t      : in-situ temperature [ITS-90]        [deg C]
! p      : sea pressure                        [dbar]
! z      : depth                               [m]

  ! Compute the absolute salinity from the practical salinity
  sa = gsw_sa_from_sp(sp,p,lon,lat)

  ! Compute the conservative temperature
  ct = gsw_ct_from_pt(sa,pt)

  ! Compute the in situ temperature
  t_from_pt  = gsw_t_from_ct(sa,ct,p)

END FUNCTION t_from_pt

!===============================================================================
PURE FUNCTION p_from_z(dpth,xlat)
!===============================================================================
! pressure from depth from saunder's formula with eos80.
! reference: saunders,peter m., practical conversion of pressure
!            to depth., j.p.o. , april 1981.
! r millard
! march 9, 1983
! check value: p_from_z=7500.004 dbars;for lat=30 deg., depth=7321.45 meters
! http://sam.ucsd.edu/sio210/propseawater/ppsw_fortran/ppsw.f
  IMPLICIT NONE
  REAL(r14) :: p_from_z
  REAL(r14), INTENT(IN) :: dpth
  REAL(r14), INTENT(IN) :: xlat
  REAL(r14), PARAMETER :: pi=3.141592654
  REAL(r14) :: plat, d, c1

  plat=abs(xlat*pi/180.)
  d=sin(plat)
  c1=5.92e-3+d**2*5.25e-3
  p_from_z=((1-c1)-sqrt(((1-c1)**2)-(8.84e-6*dpth)))/4.42e-6

END FUNCTION p_from_z

!===============================================================================
FUNCTION sa_from_sp(sp,p,long,lat)
!===============================================================================

! Calculates Absolute Salinity, SA, from Practical Salinity, SP
!
! sp     : Practical Salinity                              [unitless]
! p      : sea pressure                                    [dbar]
! long   : longitude                                       [DEG E]     
! lat    : latitude                                        [DEG N]
!
! gsw_sa_from_sp   : Absolute Salinity                     [g/kg]

IMPLICIT NONE
REAL (r14) :: sa_from_sp
REAL (r14), INTENT(IN) :: sp, long, lat, p
REAL (r14) :: gsw_sa_from_sp

sa_from_sp = gsw_sa_from_sp(sp,p,long,lat)

END FUNCTION sa_from_sp

!===============================================================================
FUNCTION pt_from_t(sa,t,p,p_ref)
!===============================================================================

! Calculates potential temperature of seawater from in-situ temperature 
!
! sa     : Absolute Salinity                               [g/kg]
! t      : in-situ temperature                             [deg C]
! p      : sea pressure                                    [dbar]
! p_ref  : reference sea pressure                          [dbar]
!
! gsw_pt_from_t : potential temperature                    [deg C]

IMPLICIT NONE
REAL (r14) :: pt_from_t
REAL (r14), INTENT(IN) :: sa, t, p, p_ref
REAL (r14) :: gsw_pt_from_t

pt_from_t = gsw_pt_from_t(sa,t,p,p_ref)

END FUNCTION pt_from_t

END MODULE gsw_pot_to_insitu
