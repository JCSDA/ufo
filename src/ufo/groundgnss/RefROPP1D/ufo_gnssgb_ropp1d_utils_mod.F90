! (C) Copyright 2017-2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Fortran module to handle ground-based gnss  observations following 
!> the ROPP (2018 Aug) implementation

module ufo_gnssgb_ropp1d_utils_mod

!use iso_c_binding
use fckit_log_module, only: fckit_log
use kinds,            only: kind_real

! ROPP data type and library subroutines
use typesizes,     only: wp => EightByteReal
use datetimetypes, only: dp
use ropp_fm_types, only: State1dFM, obs1dRefrac
use geodesy,       only: gravity, R_eff, geometric2geopotential, geopotential2geometric
use arrays,        only: callocate

implicit none
public             :: init_ropp_1d_statevec
public             :: init_ropp_1d_statevec_ad
public             :: init_ropp_1d_obvec
public             :: init_ropp_1d_obvec_tlad
public             :: ropp_tidy_up_1d
public             :: ropp_tidy_up_tlad_1d
public             :: calc_model_z
public             :: calc_station_phi
public             :: gnss_ztd_integral
private

contains

! ------------------------------------------------------------------------------------
! ------------------------------------------------------------------------------------

subroutine init_ropp_1d_statevec(step_time,rlon,rlat, temp,shum,pres,phi,lm,phi_sfc,x, iflip)
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
  character(len=250)                 :: err_msg
  real(kind=kind_real)               :: rlon_local
  integer ::  n,i,j,k
  integer, optional, intent(in)  :: iflip
  x%non_ideal     = .FALSE.
  x%direct_ion    = .FALSE.
  x%state_ok      = .TRUE.

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
  write(err_msg,'(4a9,a11)') 'lvl','temp','shum','pres','geop'

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
  write(err_msg,'("geop_sfc",f15.2)') x%geop_sfc
!------------------------------------------------
! covariance matrix, is this used by ROPP FM?
!------------------------------------------------
  x%cov_ok = .TRUE.

! Allocate memory
! For ECMWF example, Covariance matrix for temperature sigma and
! specific humidity sigma, and surface pressure. There the
! size of covariance matrix is 2 * nlevel + 1.

  n = (2*x%n_lev)+1  ! Number of elements in the state vector

  if (associated(x%cov%d)) deallocate(x%cov%d)
  call callocate(x%cov%d, n*(n+1)/2)    ! From ROPP utility library

  do i = 1, x%n_lev
     x%cov%d(i + i*(i-1)/2) = 1.0_wp
  end do

  do i = 1, x%n_lev
     j = x%n_lev + i
     x%cov%d(j + j*(j-1)/2) = 1.0_wp
  enddo

  x%cov%d(n + n*(n-1)/2) = 1.0_wp

!------------------------------
! Rest of the covariance marix
!------------------------------
  if (associated(x%cov%e)) deallocate(x%cov%e)
  if (associated(x%cov%f)) deallocate(x%cov%f)
  if (associated(x%cov%s)) deallocate(x%cov%s)

  x%cov%fact_chol = .FALSE.
  x%cov%equi_chol = 'N'

  return
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
  return

end subroutine init_ropp_1d_statevec_ad

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------

subroutine calc_station_phi(rlat, station_height, station_phi)
      
!  Description:
!     convert station geometric height to geopotential height
!     provide inputs of:
!     lat, station_height
!
!  Inputs:
!     rlat             latitude
!     station_height   station geometric height
!     
!  Outputs:
!     station_phi:     station geopotential height
!-----------------------------------------------------------------------------------
  implicit none

  real(kind=kind_real),       intent(in)    :: rlat, station_height
  real(kind=kind_real),       intent(out)   :: station_phi

  ! not used in forward operator at this time -- simulate at model levels only
  station_phi = geometric2geopotential(real(rlat,kind=wp), real(station_height,kind=wp))

end subroutine calc_station_phi

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------

subroutine calc_model_z(nlev, rlat, model_phi, model_z)
      
!  Description:
!     convert model geopotential height to geometric height
!     provide inputs of:
!     number_of_model_levels, lat, model_geopotential
!
!  Inputs:
!     nlev             number of model levels
!     rlat             latitude
!     model_phi        model geopotential height
!     
!  Outputs:
!     model_z:         model geometric height
!-----------------------------------------------------------------------------------
  implicit none

  integer,              intent(in)                     :: nlev
  real(kind=kind_real), intent(in)                     :: rlat
  real(kind=kind_real), dimension(nlev),  intent(in)   :: model_phi
  real(kind=kind_real), dimension(nlev),  intent(out)  :: model_z

  integer ilev

  ! not used in forward operator at this time -- simulate at model levels only
  do ilev = 1, nlev
    model_z(ilev) = geopotential2geometric(real(rlat,kind=wp), real(model_phi(ilev),kind=wp))
  end do

end subroutine calc_model_z

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------

subroutine init_ropp_1d_obvec(nlev, ichk, ob_time, rlat, rlon, station_phi, x, y)
      
!  Description:
!     subroutine to fill a ROPP observation vector structure
!     observation provides the inputs of:
!     lat, lon, time
!     
!     forward model will provide the simulation of refractivity
!
!  Inputs:
!     ob_time       time of the observation
!     rlat          latitude
!     rlon          longitude
!     
!  Outputs:
!     y:     Partially filled Forward model observation vector
!-----------------------------------------------------------------------------------
  implicit none
! Input model state vector
  type(State1dFM),            intent(in)    :: x
! Output observation vector
  type(Obs1dRefrac),          intent(out)   :: y

  integer,                    intent(in)    :: nlev
  integer, dimension(nlev), intent(in)      :: ichk
  real(kind=kind_real),       intent(in)    :: rlat, rlon, station_phi
  real(kind=dp),              intent(in)    :: ob_time
! Local variables
  real(kind=wp)                             :: r8lat
  real(kind=kind_real)                      :: rlon_local
  character(len=250)                        :: err_msg
  integer                                   :: i

  y%time     = real(ob_time,kind=wp)
  r8lat      = real(rlat,kind=wp)
  y%lat      = r8lat
  rlon_local = rlon
  if (rlon_local .gt. 180) rlon_local = rlon_local - 360.
  y%lon        = real(rlon_local,kind=wp)

!----------------------------------------------------
! allocate refractivity & weights
!----------------------------------------------------
  if (associated(y%refrac)) then
      deallocate(y%refrac)
      deallocate(y%weights)
      deallocate(y%geop)
      nullify(y%refrac)
      nullify(y%weights)
      nullify(y%geop)
  end if

  allocate(y%refrac(1:nlev))             ! value computed in fwd model
  allocate(y%weights(1:nlev))            ! value set in fwd model
  allocate(y%geop(1:nlev))               ! value set in fwd model

  do i=1,nlev
     if (ichk(i) .le. 0) then
        y%weights(i) = 1.0_wp              ! following t_fascod example
     else
        y%weights(i) = 0.0_wp
     end if
     y%geop(i)  = x%geop(i)
  end do

  y%refrac(:)       = 0.0_wp

  write(err_msg,'(a9,2a11,2a15)') 'ROPPyvec:','lat', 'lon', 'geop(1)', 'station_phi'
  call fckit_log%debug(err_msg)
  write(err_msg,'(9x,2f11.2,2f15.2)')   y%lat, y%lon, y%geop(1), station_phi
  call fckit_log%debug(err_msg)

!---------------------------------------------
! covariance matrix, is this used by ROPP FM?
!--------------------------------------------
  y%obs_ok = .TRUE.
  return

end subroutine init_ropp_1d_obvec                                                    

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
subroutine init_ropp_1d_obvec_tlad(nlev,  &
                         rlat,rlon,x,y,y_p)
  implicit none
! Input state vector
  type(State1dFM),                         intent(in)   :: x
! Output observation vector
  type(Obs1dRefrac),                       intent(out)  :: y,y_p

  integer,                                 intent(in)   :: nlev
  real(kind=kind_real),                    intent(in)   :: rlat,rlon
  real(kind=wp)              :: r8lat
  real(kind=kind_real)       :: rlon_local

!-------------------------------------------------------------------------
  y%time     = real(0.0, kind=wp)!)real(ob_time,kind=wp)
  r8lat      = real(rlat,kind=wp)
  y%lat      = real(rlat,kind=wp)
  rlon_local = rlon
  if (rlon_local .gt. 180) rlon_local = rlon_local - 360.
  y%lon        = real(rlon_local,kind=wp)

! allocate refractivity
!---------------------------------------------------------
  allocate(y%refrac(1:nlev))                     ! value computed in fwd model
  allocate(y%weights(1:nlev))
  allocate(y%geop(1:nlev))
  allocate(y_p%refrac(1:nlev))                   ! value computed in fwd model
  allocate(y_p%weights(1:nlev))
  allocate(y_p%geop(1:nlev))

  y%weights(:)   = 1.0
  y_p%weights(:) = 1.0
  y%refrac(:)    = 0.0
  y_p%refrac(:)  = 0.0
  y%geop(:)      = x%geop(:)
  y_p%geop(:)    = x%geop(:)

  y%obs_ok = .TRUE.
  return

end subroutine init_ropp_1d_obvec_tlad

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
subroutine gnss_ztd_integral(lm, model_refrac, model_z, ob_terr, model_ztd, l_linear)

  implicit none

  integer, intent(in)             :: lm
  real(kind=kind_real), intent(in), dimension(lm) :: model_refrac, model_z
  real(kind=kind_real), intent(in)                :: ob_terr
  real(kind=kind_real), intent(out)               :: model_ztd
  logical, intent(inout)          :: l_linear

  integer, parameter              :: max_string = 800
  character(max_string)           :: err_msg
  real            :: ddzd, c_i
  integer         :: k
  !-------------------------------------------------------------------------------
  ! integral of ROPP refractivity values to compute zenith total delay (ztd)
  !-------------------------------------------------------------------------------
  l_linear = .false.
  model_ztd = 0.0
  ddzd = 0.0
  do k = 1, lm
    if (k==1) then
      if ( model_refrac(k) <= 0 .or. model_refrac(k+1) <= 0 ) then
        write(err_msg,'(a,2es13.3)') ' unphysical refractivity ', &
            model_refrac(k), model_refrac(k+1)
        call fckit_log%info(err_msg)
        cycle
      end if
      ! lower boundary condition -- station below lowest model level
      if ( ob_terr < model_z(k) ) then
        c_i = ( log(model_refrac(k+1)) - log(model_refrac(k)) )      / &
               ( model_z(k) - model_z(k+1) )
        if ( c_i == 0 ) then
          write(err_msg,'(a,2es13.3)') ' model refractivity repeating1 ', &
              model_refrac(k), model_refrac(k+1)
          call fckit_log%info(err_msg)
          l_linear = .true.
          ddzd = model_refrac(k)*(model_z(k+1) - ob_terr)  ! use linear estimate
          cycle
        end if
        ddzd =  -1.0 * model_refrac(k)/c_i * exp(c_i*model_z(k)) *     &
               ( exp(-c_i*model_z(k+1)) - exp(-c_i*ob_terr) )
!       ddzd = model_refrac(k)*(model_z(k+1) - ob_terr) ! linear estimate
!       -- new seems wrong but looks better
!       ddzd = model_refrac(k)*(model_z(k) - ob_terr)   ! linear estimate
!       -- old
      end if
    else if ( ob_terr >= model_z(k)) then
      cycle
    else
      if ( model_refrac(k) <= 0 .or. model_refrac(k-1) <= 0 ) then
        write(err_msg,'(a,2es13.3)') ' unphysical refractivity ', &
            model_refrac(k), model_refrac(k-1)
        call fckit_log%info(err_msg)
        cycle
      end if
      ! alternate lower boundary condition -- station between levels
      if ( ob_terr >= model_z(k-1) .and. ob_terr < model_z(k) ) then
        c_i = ( log(model_refrac(k)) - log(model_refrac(k-1)) )  / &
               ( model_z(k-1) - model_z(k) )
        if ( c_i == 0 ) then
          write(err_msg,'(a,2es13.3)') ' model refractivity repeating2 ', &
              model_refrac(k), model_refrac(k-1)
          call fckit_log%info(err_msg)
          l_linear = .true.
          ddzd = model_refrac(k) * ( ob_terr - model_z(k-1) )  ! use linear estimate
          ddzd = ddzd + 0.5*(model_refrac(k)+model_refrac(k-1))*(model_z(k) - ob_terr)
          cycle
        end if
        ddzd =  -1.0 * model_refrac(k-1)/c_i * exp(c_i*ob_terr) * &
               ( exp(-c_i*model_z(k)) - exp(-c_i*ob_terr) )
!       ddzd = model_refrac(k) * ( ob_terr - model_z(k-1) )  ! linear
!       estimate -- new
!       ddzd = model_refrac(k) * ( model_z(k) - ob_terr )  ! linear
!       estimate -- old
        ddzd = ddzd + 0.5*(model_refrac(k)+model_refrac(k-1))*(model_z(k) - ob_terr)
      else ! normal integral up the column
        ddzd = 0.5*(model_refrac(k)+model_refrac(k-1))*(model_z(k) - model_z(k-1))
      end if
    end if
    model_ztd = model_ztd + ddzd
  enddo

  ! a very bad place to put this suggestions?
  model_ztd = model_ztd * 1.e-6

end subroutine gnss_ztd_integral

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
subroutine ropp_tidy_up_1d(x,y)
  implicit none
  type(state1dfm),   intent(inout) :: x
  type(obs1drefrac), intent(inout) :: y
! x
  if (associated(x%temp)) deallocate(x%temp)
  if (associated(x%shum)) deallocate(x%shum)
  if (associated(x%pres)) deallocate(x%pres)
  if (associated(x%geop)) deallocate(x%geop)
! y
  if (associated(y%refrac)) deallocate(y%refrac)
  if (associated(y%geop)) deallocate(y%geop)
  if (associated(y%weights)) deallocate(y%weights)

end subroutine ropp_tidy_up_1d

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
subroutine ropp_tidy_up_tlad_1d(x,x_p,y,y_p)
  implicit none
  type(state1dfm),   intent(inout) :: x,x_p  ! x_p can be either x_tl or x_ad
  type(obs1drefrac), intent(inout) :: y,y_p  ! y_p can be either y_tl or y_ad
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
  deallocate(y%refrac)
  deallocate(y_p%refrac)

  deallocate(y%weights)
  deallocate(y_p%weights)

  deallocate(y%geop)
  deallocate(y_p%geop)
 return

end subroutine ropp_tidy_up_tlad_1d

!-------------------------------------------------------------------------
     
end module ufo_gnssgb_ropp1d_utils_mod
