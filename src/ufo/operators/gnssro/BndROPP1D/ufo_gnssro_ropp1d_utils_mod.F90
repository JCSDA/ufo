! (C) Copyright 2017-2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Fortran module to handle gnssro bending angle observations following 
!> the ROPP (2018 Aug) implementation

module ufo_gnssro_ropp1d_utils_mod

!use iso_c_binding
use fckit_log_module, only: fckit_log
use kinds,            only: kind_real

! ROPP data type and library subroutines
use typesizes,     only: wp => EightByteReal
use datetimetypes, only: dp
use ropp_fm_types, only: State1dFM, obs1dBangle
use geodesy,       only: gravity, R_eff, geometric2geopotential
use arrays,        only: callocate

implicit none
public             :: init_ropp_1d_statevec
public             :: init_ropp_1d_statevec_ad
public             :: init_ropp_1d_obvec
public             :: init_ropp_1d_obvec_tlad
public             :: ropp_tidy_up_1d
public             :: ropp_tidy_up_tlad_1d
private

contains

! ------------------------------------------------------------------------------------
! ------------------------------------------------------------------------------------

subroutine init_ropp_1d_statevec(step_time,rlon,rlat, temp,shum,pres,phi,lm,phi_sfc,x, &
                                 iflip, non_ideal)
!  Description:
!     subroutine to fill a ROPP state vector structure with
!     model background fields: Temperature, pressure, specific
!     humidity at full-levels, and surface geopotential height.
!
!  Inputs:
!     temp   background temperature
!     shum   background specific humidity
!     pres   background pressure
!     phi    geopotential height
!     phi_sfc  terrain geopotential of background field
!     lm       number of vertical levels in  the background
!
!  Outputs:
!     x      Forward model state vector
! ------------------------------------------------------------------------------------
  implicit none

  type(State1dFM),                     intent(out) :: x
  real(kind=dp),                       intent(in)  :: step_time
  real(kind=kind_real),                intent(in)  :: rlat, rlon
  real(kind=kind_real),                intent(in)  :: phi_sfc
  integer,                             intent(in)  :: lm
  real(kind=kind_real), dimension(lm), intent(in)  :: temp,shum,pres,phi

! Local variables
  character(len=250)                 :: record
  real(kind=kind_real)               :: rlon_local
  integer ::  n,i,j,k
  integer, optional, intent(in)  :: iflip
  integer, intent(in)                :: non_ideal

  if (non_ideal .eq. 1) then
    x%non_ideal     = .TRUE.
  else
    x%non_ideal     = .FALSE.
  endif
  x%direct_ion    = .FALSE.
  x%state_ok      = .TRUE.
  x%new_bangle_op = .TRUE.     ! activate ROPP v8 new interpolation scheme

! ROPP Longitude value is -180.0 to 180.0
  x%lat      = real(rlat,kind=wp)
  rlon_local = rlon
  if (rlon_local .gt. 180) rlon_local = rlon_local - 360.
  x%lon  = real(rlon_local,kind=wp)
  x%time = real(step_time, kind=wp)


! Number of levels in background profile.  What about (lm+1) field ?
  x%n_lev = lm

!--------------------------------------------------------------
! allocate arrays for temperature, specific humidity, pressure
! and geopotential height data
!--------------------------------------------------------------
  if (associated(x%temp)) deallocate(x%temp) 
  if (associated(x%shum)) deallocate(x%shum)
  if (associated(x%pres)) deallocate(x%pres)
  if (associated(x%geop)) deallocate(x%geop)

  allocate(x%temp(x%n_lev))
  allocate(x%shum(x%n_lev))
  allocate(x%pres(x%n_lev))
  allocate(x%geop(x%n_lev))

!----------------------------------------------------
! ROPP FM requires vertical height profile to be of the ascending order.
! (see ropp_io_ascend ( ROdata )). So we need to flip the data.
!----------------------------------------------------
  write(record,'(4a9,a11)') 'lvl','temp','shum','pres','geop'

 n = lm
 if ( present(iflip) .and. iflip .eq. 1) then
   do k = 1, lm
      x%temp(n) = real(temp(k),kind=wp)
      x%shum(n) = real(shum(k),kind=wp)
      x%pres(n) = real(pres(k),kind=wp)
      x%geop(n) = real(phi(k),kind=wp)
      n = n - 1
   end do
else
   do k = 1, lm
      x%temp(k) = real(temp(k),kind=wp)
      x%shum(k) = real(shum(k),kind=wp)
      x%pres(k) = real(pres(k),kind=wp)
      x%geop(k) = real(phi(k),kind=wp)
   end do
end if

! sufrace geopotential height value
  x%geop_sfc = real(phi_sfc,kind=wp)
  write(record,'("geop_sfc",f15.2)') x%geop_sfc
!------------------------------------------------
! make sure to deallocate all in memory
  if (associated(x%cov%d)) deallocate(x%cov%d)
  if (associated(x%cov%e)) deallocate(x%cov%e)
  if (associated(x%cov%f)) deallocate(x%cov%f)
  if (associated(x%cov%s)) deallocate(x%cov%s)

end subroutine init_ropp_1d_statevec

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------

subroutine init_ropp_1d_statevec_ad(temp_d,shum_d,pres_d,phi_d,lm,x_ad, iflip)

!  Description:
!     subroutine to fill a ROPP state vector structure with
!     model background fields: Temperature, pressure, specific
!     humidity at full-levels, and surface geopotential height.
!
!  Inputs:
!     temp   background temperature
!     shum   background specific humidity
!     pres   background pressure
!     phi    geopotential height
!     phi_sfc  terrain geopotential of background field
!     lm       number of vertical levels in  the background
!
!  Outputs:
!     x      Forward model state vector
!
! ###############################################################
  type(State1dFM),     intent(inout) :: x_ad
  integer,             intent(in)    :: lm
  real(kind=kind_real), dimension(lm), intent(inout)    :: temp_d,shum_d,pres_d,phi_d

! Local variables
  integer   ::   n,k
  integer, optional, intent(in) :: iflip
!-------------------------------------------------------------------------

  n = lm

  if ( present(iflip) .and. iflip .eq. 1) then
    do k = 1, lm

!!     x_tl%temp(n,:) = real(temp_d(k),kind=wp)
       temp_d(k) = temp_d(k) + real(x_ad%temp(n),kind=kind_real)
       x_ad%temp(n) = 0.0_wp

!!     x_tl%shum(n,:) = real(shum_d(k),kind=wp)
       shum_d(k) = shum_d(k) + real(x_ad%shum(n),kind=kind_real)
       x_ad%shum(n) = 0.0_wp

!!     x_tl%pres(n,:) = real(pres_d(k),kind=wp)
       pres_d(k) = pres_d(k) + real(x_ad%pres(n),kind=kind_real)
       x_ad%pres(n) = 0.0_wp 

!!     x_tl%geop(n,:) = real(phi_d(k),kind=wp)
       phi_d(k) = phi_d(k) + real(x_ad%geop(n),kind=kind_real)
       x_ad%geop(n) = 0.0_wp

       n = n -1
    end do
  else
    do k = 1, lm
       temp_d(k) = temp_d(k) + real(x_ad%temp(k),kind=kind_real)
       x_ad%temp(k) = 0.0_wp
       shum_d(k) = shum_d(k) + real(x_ad%shum(k),kind=kind_real)
       x_ad%shum(k) = 0.0_wp
       pres_d(k) = pres_d(k) + real(x_ad%pres(k),kind=kind_real)
       x_ad%pres(k) = 0.0_wp
       phi_d(k) = phi_d(k) + real(x_ad%geop(k),kind=kind_real)
       x_ad%geop(k) = 0.0_wp
     end do

  end if

end subroutine init_ropp_1d_statevec_ad

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------

subroutine init_ropp_1d_obvec(nvprof,obs_impact,ichk,ob_time,rlat,rlon,roc,undulat,y)
      
!  Description:
!     subroutine to fill a ROPP observation vector structure
!     observation provides the inputs of:
!     impact parameter, lat, lon, time, radius of curvature
!     
!     forward model will provide the simulation of bending angle
!
!  Inputs:
!     obs_impact    impact parameter
!     ob_time       time of the observation
!     rlat          latitude
!     rlon          longitude
!     roc           radius of curvature
!     undulat       undulation correction for radius of curvature
!     
!  Outputs:
!     y:     Partially filled Forward model observation vector
!-----------------------------------------------------------------------------------
  implicit none
! Output state vector
  type(Obs1dBangle),          intent(out)   :: y

  integer,                    intent(in)    :: nvprof
  integer, dimension(nvprof), intent(in)    :: ichk
  real(kind=kind_real), dimension(nvprof), intent(in)    :: obs_impact
  real(kind=kind_real),       intent(in)    :: rlat, rlon
  real(kind=kind_real),       intent(in)    :: roc, undulat
  real(kind=dp),              intent(in)    :: ob_time
! Local variables
  real(kind=wp)                             :: r8lat
  real(kind=kind_real)                      :: rlon_local
  character(len=250)                        :: record
  integer                                   :: i

  y%time     = real(ob_time,kind=wp)
  r8lat      = real(rlat,kind=wp)
  y%lat      = r8lat
  rlon_local = rlon
  if (rlon_local .gt. 180) rlon_local = rlon_local - 360.
  y%lon        = real(rlon_local,kind=wp)
  y%nobs       = nvprof
  y%g_sfc      = gravity(r8lat, 0.0_wp)          ! 2nd argument is height(m) above the sfc
  y%r_curve    = real(roc,kind=wp)             ! ROPP rad of curve (m)
  y%undulation = real(undulat,kind=wp)      ! ROPP undulation corr for rad of curve (m)
  y%r_earth    = r_eff(r8lat)

!----------------------------------------------------
! allocate bending angle, impact parameter & weights
!----------------------------------------------------
  if (associated(y%bangle)) then
      deallocate(y%bangle)
      deallocate(y%impact)
      deallocate(y%weights)
      nullify(y%bangle)
      nullify(y%impact)
      nullify(y%weights)
  end if

  allocate(y%bangle(1:nvprof))                     ! value computed in fwd model
  allocate(y%impact(1:nvprof))
  allocate(y%weights(1:nvprof))                    ! value set in fwd model

  do i=1,nvprof
     y%impact(i) = real(obs_impact(i),kind=wp)  ! ROPP expects impact parameter in meters
     if (ichk(i) .le. 0) then
        y%weights(i) = 1.0_wp                    ! following t_fascod example
     else
        y%weights(i) = 0.0_wp
     end if
  end do

  y%bangle(:)  = 0.0_wp                    ! following t_fascod example

  write(record,'(a9,2a11,3a15)') 'ROPPyvec:','lat', 'lon', 'g_sfc', 'roc', 'r_earth_eff'
  write(record,'(9x,2f11.2,f15.6,2f15.2)')   y%lat, y%lon, y%g_sfc, y%r_curve, y%r_earth

  y%obs_ok = .TRUE.

end subroutine init_ropp_1d_obvec                                                    

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
subroutine init_ropp_1d_obvec_tlad(iloop,nvprof,obs_impact,  &
                         rlat,rlon,roc,undulat,y,y_p)
  implicit none
! Output state vector
  type(Obs1dBangle),                       intent(out)  :: y,y_p

  integer,                                 intent(in)   :: iloop
  integer,                                 intent(in)   :: nvprof
  real(kind=kind_real), dimension(nvprof), intent(in)   :: obs_impact
  real(kind=kind_real),                    intent(in)   :: rlat,rlon
  real(kind=kind_real),                    intent(in)   :: roc, undulat
  real(kind=wp)              :: r8lat
  integer                    :: i
  real(kind=kind_real)       :: rlon_local

!-------------------------------------------------------------------------
  y%time     = real(0.0, kind=wp)!)real(ob_time,kind=wp)
  r8lat      = real(rlat,kind=wp)
  y%lat      = real(rlat,kind=wp)
  rlon_local = rlon
  if (rlon_local .gt. 180) rlon_local = rlon_local - 360.
  y%lon        = real(rlon_local,kind=wp)
  y%nobs       = nvprof
  y%g_sfc      = gravity(r8lat, 0.0_wp)          ! 2nd argument is height(m) above the sfc
  y%r_curve    = real(roc,kind=wp)             ! ROPP rad of curve (m)
  y%undulation = real(undulat,kind=wp)      ! ROPP undulation corr for rad of curve (m)
  y%r_earth    = r_eff(r8lat)

! allocate bending angle, impact parameter & weights
!---------------------------------------------------------
  allocate(y%bangle(1:nvprof))                     ! value computed in fwd model
  allocate(y%impact(1:nvprof))
  allocate(y%weights(1:nvprof))
  allocate(y_p%bangle(1:nvprof))                     ! value computed in fwd model
  allocate(y_p%impact(1:nvprof))
  allocate(y_p%weights(1:nvprof))

! number of points in profile
  y%nobs    = nvprof
  y_p%nobs  = nvprof
 
  do i=1,nvprof
     y_p%impact(i)  = real(obs_impact(i),kind=wp)
     y%impact(i)    = real(obs_impact(i),kind=wp)
     y%weights(i)   = 1.0
     y_p%weights(i) = 1.0
     y%bangle(i)    = 0.0
     y_p%bangle(i)  = 0.0
  end do

  y%obs_ok = .TRUE.

end subroutine init_ropp_1d_obvec_tlad

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
subroutine ropp_tidy_up_1d(x,y)
  implicit none
  type(state1dfm),   intent(inout) :: x
  type(obs1dbangle), intent(inout) :: y
! x
  if (associated(x%temp)) deallocate(x%temp)
  if (associated(x%shum)) deallocate(x%shum)
  if (associated(x%pres)) deallocate(x%pres)
  if (associated(x%geop)) deallocate(x%geop)

! y
  if (associated(y%impact)) deallocate(y%impact)
  if (associated(y%bangle)) deallocate(y%bangle)
  if (associated(y%weights))  deallocate(y%weights)

end subroutine ropp_tidy_up_1d

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
subroutine ropp_tidy_up_tlad_1d(x,x_p,y,y_p)
  implicit none
  type(state1dfm),   intent(inout) :: x,x_p  ! x_p can be either x_tl or x_ad
  type(obs1dbangle), intent(inout) :: y,y_p  ! y_p can be either y_tl or y_ad
! x
  deallocate(x%temp)
  deallocate(x%shum)
  deallocate(x%pres)
  deallocate(x%geop)

  deallocate(x_p%temp)
  deallocate(x_p%shum)
  deallocate(x_p%pres)
  deallocate(x_p%geop)

! y
  deallocate(y%impact)
  deallocate(y%bangle)
  deallocate(y%weights)

  deallocate(y_p%impact)
  deallocate(y_p%bangle)
  deallocate(y_p%weights)

end subroutine ropp_tidy_up_tlad_1d

!-------------------------------------------------------------------------
     
end module ufo_gnssro_ropp1d_utils_mod
