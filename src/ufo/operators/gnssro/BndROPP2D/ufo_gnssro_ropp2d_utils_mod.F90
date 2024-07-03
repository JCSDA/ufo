! (C) Copyright 2017-2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Fortran module to handle gnssro bending angle observations following 
!> the ROPP (2018 Aug) implementation

module ufo_gnssro_ropp2d_utils_mod

!use iso_c_binding
use fckit_log_module, only: fckit_log
use kinds,            only: kind_real

! ROPP data type and library subroutines
use typesizes,     only: wp => EightByteReal
use datetimetypes, only: dp
use ropp_fm_types, only: State2dFM, obs1dBangle
use geodesy,       only: gravity, R_eff, geometric2geopotential
use arrays,        only: callocate

implicit none
public             :: init_ropp_2d_statevec
public             :: init_ropp_2d_statevec_ad
public             :: init_ropp_2d_obvec
public             :: init_ropp_2d_obvec_tlad
public             :: ropp_tidy_up_2d
public             :: ropp_tidy_up_tlad_2d
private

contains

! ------------------------------------------------------------------------------------
! ------------------------------------------------------------------------------------
subroutine init_ropp_2d_statevec(rlon,rlat,temp,shum,pres,phi,lm, x,  &
                                 n_horiz, dtheta, iflip,non_ideal)

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
!
!     phi_sfc  terrain geopotential of background field
!     lm       number of vertical levels in  the background
!
!  Outputs:
!     x      Forward model state vector
!
! ###############################################################
  implicit none
! Output state vector
  type(State2dFM),      intent(out)   :: x
  integer,              intent(in)    :: lm, n_horiz,non_ideal
  real(kind=kind_real), intent(in)    :: dtheta
  real(kind=kind_real), dimension(n_horiz), intent(in)     ::   rlon, rlat
  real(kind=kind_real), dimension(lm,n_horiz), intent(in)  :: temp,shum,pres,phi
! Local variables
  integer :: n,k
  integer, optional, intent(in)  :: iflip
!-------------------------------------------------------------------------
! number of profiles in plane
  x%n_horiz = n_horiz
  x%dtheta  = dtheta
  if (non_ideal .eq. 1) then
    x%non_ideal = .TRUE.
  else
    x%non_ideal = .FALSE.
  endif
! Number of levels in background profile.  What about (lm+1) field ?
  x%n_lev=lm

!-------------------------------------------------------------
! allocate arrays for temperature, specific humidity, pressure
! and geopotential height data and additional 2d fields
!-------------------------------------------------------------
  allocate(x%temp(x%n_lev,x%n_horiz))
  allocate(x%shum(x%n_lev,x%n_horiz))
  allocate(x%pres(x%n_lev,x%n_horiz))
  allocate(x%geop(x%n_lev,x%n_horiz))

! needed in 2D
  allocate(x%refrac(x%n_lev,x%n_horiz))
  allocate(x%nr(x%n_lev,x%n_horiz))
  allocate(x%lat(x%n_horiz))
  allocate(x%lon(x%n_horiz))

  x%lat(:) = real(rlat(:),kind=wp)
  x%lon(:) = real(rlon(:),kind=wp)
  where (x%lon .gt. 180.0) x%lon = x%lon -360.0


! TEMPORARY!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!----------------------------------------------------
! ROPP FM requires vertical height profile to be of the ascending order.
! (see ropp_io_ascend ( ROdata )). So we need to flip the data.
!----------------------------------------------------

  n = lm
  if ( present(iflip) .and. iflip .eq. 1) then

    do k = 1, lm
       x%temp(n,:) = real(temp(k,:),kind=wp)
       x%shum(n,:) = real(shum(k,:),kind=wp)
       x%pres(n,:) = real(pres(k,:),kind=wp)
       x%geop(n,:) = real(phi(k,:),kind=wp)
       n = n - 1
    end do

  else
    do k = 1, lm
       x%temp(k,:) = real(temp(k,:),kind=wp)
       x%shum(k,:) = real(shum(k,:),kind=wp)
       x%pres(k,:) = real(pres(k,:),kind=wp)
       x%geop(k,:) = real(phi(k,:),kind=wp)
    end do

  end if
  return
end subroutine init_ropp_2d_statevec

! ------------------------------------------------------------------------------

subroutine init_ropp_2d_statevec_ad(temp_d,shum_d,pres_d,phi_d,lm,x_ad,n_horiz,iflip)

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
!
!     phi_sfc  terrain geopotential of background field
!     lm       number of vertical levels in  the background
!
!  Outputs:
!     x      Forward model state vector
!
! ###############################################################
  implicit none

! Function arguments

! Output state vector

  type(State2dFM),           intent(inout) :: x_ad
  integer,                   intent(in)    :: lm, n_horiz
  real(kind=kind_real), &
       dimension(lm,n_horiz), intent(inout) :: temp_d,shum_d,pres_d,phi_d

! Local variables
  integer ::  n,j,k
  integer, optional, intent(in)  :: iflip
!-------------------------------------------------------------------------
  n = lm
  x_ad%n_horiz = n_horiz

  if ( present(iflip) .and. iflip .eq. 1) then

    do k = 1, lm
       do j =  1, x_ad%n_horiz
!!!        x_tl%temp(n,:) = real(temp_d(k),kind=wp)
           temp_d(k,j) = temp_d(k,j) + real(x_ad%temp(n,j),kind=kind_real)
           x_ad%temp(n,j) = 0.0_wp
        
!!!        x_tl%shum(n,:) = real(shum_d(k),kind=wp)
           shum_d(k,j) = shum_d(k,j) + real(x_ad%shum(n,j),kind=kind_real)
           x_ad%shum(n,j) = 0.0_wp
        
!!!        x_tl%pres(n,:) = real(pres_d(k),kind=wp)
           pres_d(k,j) = pres_d(k,j) + real(x_ad%pres(n,j),kind=kind_real)
           x_ad%pres(n,j) = 0.0_wp
        
!!!        x_tl%geop(n,:) = real(phi_d(k),kind=wp)
           phi_d(k,j) = phi_d(k,j) + real(x_ad%geop(n,j),kind=kind_real)
           x_ad%geop(n,j) = 0.0_wp
      enddo

      n = n - 1
    end do
  else
    do k = 1, lm
       do j =  1, x_ad%n_horiz
           temp_d(k,j) = temp_d(k,j) + real(x_ad%temp(k,j),kind=kind_real)
           x_ad%temp(k,j) = 0.0_wp
           shum_d(k,j) = shum_d(k,j) + real(x_ad%shum(k,j),kind=kind_real)
           x_ad%shum(k,j) = 0.0_wp
           pres_d(k,j) = pres_d(k,j) + real(x_ad%pres(k,j),kind=kind_real)
           x_ad%pres(k,j) = 0.0_wp
           phi_d(k,j) = phi_d(k,j) + real(x_ad%geop(k,j),kind=kind_real)
           x_ad%geop(k,j) = 0.0_wp
      enddo
    end do
  end if
  return
 end subroutine init_ropp_2d_statevec_ad
!-------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------
subroutine init_ropp_2d_obvec(nvprof,obs_impact,rlat,rlon,roc,undulat,y)

!  Description:
!     subroutine to fill a ROPP observation vector structure
!     observation provides the inputs of:
!     impact parameter, lat, lon, time, radius of curvature
!
!     forward model will provide the simulation of bending angle
!
!  Inputs:
!     obs_impact    impact parameter
!     rlat          latitude
!     rlon          longitude
!     roc           radius of curvature
!     undulat       undulation correction for radius of curvature
!
!  Outputs:
!     y:     Partially filled Forward model observation vector
!
! ###############################################################
  implicit none
! Output state vector
  type(Obs1dBangle),                       intent(out)   :: y

  integer,                                 intent(in)    :: nvprof
  real(kind=kind_real), dimension(nvprof), intent(in)    :: obs_impact
  real(kind=kind_real),                    intent(in)    :: rlat,rlon
  real(kind=kind_real),                    intent(in)    :: roc, undulat

  real(kind=wp)              :: r8lat
  integer                    :: i
  real(kind=kind_real)       :: rlon_local

!-------------------------------------------------------------------------
  r8lat = real(rlat,kind=wp)
  y%lat = r8lat
  rlon_local = rlon
  if (rlon_local .gt. 180.) rlon_local = rlon_local - 360.0 ! ROPP Longitude value -180 to 180
  y%lon        = real(rlon_local,kind=wp)
  y%g_sfc      = gravity(r8lat, 0.0_wp)          ! 2nd argument is height(m) above the sfc
  y%r_curve    = real(roc,kind=wp)             ! ROPP rad of curve (m)
  y%undulation = real(undulat,kind=wp)      ! ROPP undulation corr for rad of curve (m)
  y%r_earth    = r_eff(r8lat)

!---------------------------------------------------
! allocate bending angle, impact parameter & weights
!----------------------------------------------------
  allocate(y%bangle(1:nvprof))                     ! value computed in fwd model
  allocate(y%impact(1:nvprof))
! needed for 2D
  allocate(y%a_path(1:nvprof,2))
  allocate(y%rtan(1:nvprof))

! number of points in profile
  y%nobs = nvprof

  do i=1,nvprof
     y%impact(i) = real(obs_impact(i),kind=wp)  ! ROPP expects impact parameter in meters
  end do

  return

end subroutine init_ropp_2d_obvec

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
subroutine init_ropp_2d_obvec_tlad(iloop,nvprof,obs_impact,  &
                                rlat,rlon,roc,undulat,y,y_p)

!  Description:
!     subroutine to fill a ROPP observation vector structure
!     observation provides the inputs of:
!     impact parameter, lat, lon, time, radius of curvature
!
!     forward model will provide the simulation of bending angle
!
!  Inputs:
!     obs_impact    impact parameter
!     rlat          latitude
!     rlon          longitude
!     roc           radius of curvature
!     undulat       undulation correction for radius of curvature
!
!  Outputs:
!     y:     Partially filled Forward model observation vector
!-------------------------------------------------------------------------
  implicit none
! Output state vector
  type(Obs1dBangle),          intent(out)   :: y,y_p

  integer,                    intent(in)    :: iloop
  integer,                    intent(in)    :: nvprof
  real(kind=kind_real),    dimension(nvprof), intent(in)    :: obs_impact
  real(kind=kind_real),    intent(in)    :: rlat,rlon
  real(kind=kind_real),    intent(in)    :: roc, undulat

  real(kind=wp)              :: r8lat
  integer                    :: i
  real(kind=kind_real)       :: rlon_local

!-------------------------------------------------------------------------
  r8lat = real(rlat,kind=wp)
  y%lat = r8lat
  rlon_local = rlon
  if (rlon_local .gt. 180.) rlon_local = rlon_local - 360.0 ! ROPP Longitude value -180 to 180
  y%lon = real(rlon_local,kind=wp)
  y%g_sfc = gravity(r8lat, 0.0_wp)          ! 2nd argument is height(m) above the sfc
  y%r_curve = real(roc,kind=wp)             ! ROPP rad of curve (m)
  y%undulation = real(undulat,kind=wp)      ! ROPP undulation corr for rad of curve (m)
  y%r_earth = r_eff(r8lat)

!--------------------------------------------------------
! allocate bending angle, impact parameter & weights
!---------------------------------------------------------
!  if (iloop ==1) then

     allocate(y%bangle(1:nvprof))                     ! value computed in fwd model
     allocate(y%impact(1:nvprof))
!    needed for 2D
     allocate(y%a_path(1:nvprof,2))
     allocate(y%rtan(1:nvprof))
!    TL code
     allocate(y_p%bangle(1:nvprof))                     ! value computed in fwd model
     allocate(y_p%impact(1:nvprof))
!    needed for 2D
     allocate(y_p%a_path(1:nvprof,2))
     allocate(y_p%rtan(1:nvprof))

!  endif


! number of points in profile
  y%nobs    = nvprof
  y_p%nobs = nvprof

  do i = 1 ,nvprof
     y%impact(i)    = real(obs_impact(i),kind=wp)  ! ROPP expects impact parameter in meters
     y_p%impact(i)  = real(obs_impact(i),kind=wp)  ! ROPP expects impact parameter in meters
  end do

  return

end subroutine init_ropp_2d_obvec_tlad


!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
subroutine ropp_tidy_up_2d(x,y)
  implicit none
  type(state2dfm),   intent(inout) :: x
  type(obs1dbangle), intent(inout) :: y
! x
  if (associated(x%temp))   deallocate(x%temp)
  if (associated(x%shum))   deallocate(x%shum)
  if (associated(x%pres))   deallocate(x%pres)
  if (associated(x%geop))   deallocate(x%geop)
  if (associated(x%nr))     deallocate(x%nr)
  if (associated(x%refrac)) deallocate(x%refrac)
  if (associated(x%lat))    deallocate(x%lat)
  if (associated(x%lon))    deallocate(x%lon)
! y
  if (associated(y%impact)) deallocate(y%impact)
  if (associated(y%bangle)) deallocate(y%bangle)
  if (associated(y%a_path)) deallocate(y%a_path)
  if (associated(y%rtan))   deallocate(y%rtan)
 
end subroutine ropp_tidy_up_2d

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
subroutine ropp_tidy_up_tlad_2d(x,x_p,y,y_p)
  implicit none
  type(state2dfm),   intent(inout) :: x,x_p  ! x_p can be either x_tl or x_ad
  type(obs1dbangle), intent(inout) :: y,y_p  ! y_p can be either y_tl or y_ad
! x
  deallocate(x%temp)
  deallocate(x%shum)
  deallocate(x%pres)
  deallocate(x%geop)
  deallocate(x%nr)
  deallocate(x%refrac)
  deallocate(x%lat)
  deallocate(x%lon)
  deallocate(x_p%temp)
  deallocate(x_p%shum)
  deallocate(x_p%pres)
  deallocate(x_p%geop)
  deallocate(x_p%nr)
  deallocate(x_p%refrac)
  deallocate(x_p%lat)
  deallocate(x_p%lon)
! y
  deallocate(y%impact)
  deallocate(y%bangle)
  deallocate(y%a_path)
  deallocate(y%rtan)
  deallocate(y_p%impact)
  deallocate(y_p%bangle)
  deallocate(y_p%a_path)
  deallocate(y_p%rtan)

  return

end subroutine ropp_tidy_up_tlad_2d
     
end module ufo_gnssro_ropp2d_utils_mod
