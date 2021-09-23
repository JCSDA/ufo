!==========================================================================
module gnssro_mod_obserror
!==========================================================================

use iso_c_binding
use kinds
use gnssro_mod_constants
use ufo_roobserror_utils_mod
use fckit_log_module, only: fckit_log

private
public :: bending_angle_obserr_ECMWF, bending_angle_obserr_NRL
public :: bending_angle_obserr_NBAM, refractivity_obserr_NBAM
public :: gnssro_obserr_avtemp, gnssro_obserr_latitude

contains
subroutine bending_angle_obserr_ECMWF(obsImpH, obsValue, nobs,  obsErr, QCflags, missing)
implicit none
integer,                         intent(in)  :: nobs
real(kind_real), dimension(nobs),intent(in)  :: obsImpH, obsValue
integer(c_int),  dimension(nobs),intent(in)  :: QCflags(:)
real(kind_real), dimension(nobs),intent(out) :: obsErr
real(kind_real)                  :: H_km, missing
integer :: i

obsErr = missing

do i = 1, nobs
   
if (QCflags(i) .eq. 0) then

    H_km  = obsImpH(i)/1000.0_kind_real

    if ( H_km <= 10.0 ) then
       obsErr(i) = (H_km*1.25 + (10-H_km)*20)/10.0
       obsErr(i) = obsErr(i)/100.0*obsValue(i)
    else if ( H_km > 10.0 .and. H_km <= 32.0 ) then
      obsErr(i) = 1.25/100.0*obsValue(i) 
    else
      obsErr(i) = 3.0*1e-6
    end if
end if

end do

end subroutine bending_angle_obserr_ECMWF
!----------------------------------------

subroutine bending_angle_obserr_NRL(obsLat, obsImpH, obsValue, nobs,  obsErr, QCflags, missing)
use ufo_constants_mod, only: deg2rad
implicit none
integer,                         intent(in)  :: nobs
real(kind_real), dimension(nobs),intent(in)  :: obsImpH, obsValue, obsLat
integer(c_int),  dimension(nobs),intent(in)  :: QCflags(:)
real(kind_real), dimension(nobs),intent(out) :: obsErr
real(kind_real):: H_m, missing
real(kind_real):: lat_in_rad,trop_proxy,damping_factor,errfac
real(kind_real):: max_sfc_error, min_ba_error

integer :: i

obsErr = missing 
max_sfc_error = 20.0   ! %
min_ba_error  = 1.25   ! %

do i = 1, nobs

if (QCflags(i) .eq. 0) then
   
   H_m  = obsImpH(i)
   lat_in_rad = deg2rad * obsLat(i)
   trop_proxy = 8666.66 + 3333.33*cos(2.0*lat_in_rad)
   damping_factor = 0.66 + cos(lat_in_rad)/3.0
   errfac = max_sfc_error*damping_factor*(trop_proxy - H_m)/trop_proxy
   errfac = max(min_ba_error, errfac)
   obsErr(i) =  max(obsValue(i)*errfac/100.0, 3.0*1e-6) ! noise floor at top of RO profile

end if

end do

end subroutine bending_angle_obserr_NRL
!--------------------------------------

subroutine  bending_angle_obserr_NBAM(obsLat, obsImpH, obsSaid, nobs, obsErr, QCflags, missing)
implicit none
integer,                         intent(in)  :: nobs
real(kind_real), dimension(nobs),intent(in)  :: obsImpH, obsLat
integer(c_int),  dimension(nobs),intent(in)  :: obsSaid, QCflags(:)
real(kind_real), dimension(nobs),intent(out) ::  obsErr
real(kind_real)                 :: H_km, missing

integer :: i

obsErr = missing

do i = 1, nobs

if (QCflags(i) .eq. 0) then

   H_km  = obsImpH(i)/1000.0_kind_real
   if( (ObsSaid(i)==41).or.(ObsSaid(i)==722).or.(ObsSaid(i)==723).or.   &
       (ObsSaid(i)==4).or.(ObsSaid(i)==42).or.(ObsSaid(i)==3).or.       &
       (ObsSaid(i)==5).or.(ObsSaid(i)==821.or.(ObsSaid(i)==421)).or.    &
       (ObsSaid(i)==440).or.(ObsSaid(i)==43)) then
       if( abs(obsLat(i))>= 40.00 ) then
         if(H_km>12.0) then
           obsErr(i)=0.19032 +0.287535 *H_km-0.00260813*H_km**2
         else
           obsErr(i)=-3.20978 +1.26964 *H_km-0.0622538 *H_km**2
         endif
       else
         if(H_km>18.) then
           obsErr(i)=-1.87788 +0.354718 *H_km-0.00313189 *H_km**2
         else
           obsErr(i)=-2.41024 +0.806594 *H_km-0.027257 *H_km**2
         endif
       endif

   else !!!! CDAAC processing
     if( abs(obsLat(i))>= 40.00 ) then
       if ( H_km>12.00 ) then
          obsErr(i)=-0.685627 +0.377174 *H_km-0.00421934 *H_km**2
       else
          obsErr(i)=-3.27737 +1.20003 *H_km-0.0558024 *H_km**2
       endif
     else
       if( H_km>18.00 ) then
          obsErr(i)=-2.73867 +0.447663 *H_km-0.00475603 *H_km**2
       else
          obsErr(i)=-3.45303 +0.908216 *H_km-0.0293331 *H_km**2
       endif
     endif
     obsErr(i) = 0.001 /abs(exp(obsErr(i)))

  endif

end if

end do

end subroutine bending_angle_obserr_NBAM
!---------------------------------------

subroutine refractivity_obserr_NBAM(obsLat, obsZ, nobs, obsErr, QCflags,missing)
implicit none
integer,                         intent(in)  :: nobs
real(kind_real), dimension(nobs),intent(in)  :: obsLat, obsZ
real(kind_real), dimension(nobs),intent(out) :: obsErr
integer(c_int),  dimension(nobs),intent(in)  :: QCflags(:)
real(kind_real)                   :: H_km, missing

integer :: i

obsErr = missing

do i = 1, nobs

  if (QCflags(i) .eq. 0) then
     H_km  = obsZ(i)/1000.0_kind_real
     if( abs(obsLat(i))>= 20.0 ) then
         obsErr(i)=-1.321+0.341*H_km-0.005*H_km**2
     else
       if(H_km > 10.0) then
          obsErr(i)=2.013-0.060*H_km+0.0045*H_km**2
       else
          obsErr(i)=-1.18+0.058*H_km+0.025*H_km**2
       endif
     endif
     obsErr(i) = 1.0_kind_real/abs(exp(obsErr(i)))
  end if
end do

end subroutine refractivity_obserr_NBAM


subroutine gnssro_obserr_avtemp(nobs, n_horiz, rmatrix_filename, obsSatid, obsOrigC, nlevs, &
                                air_temperature, geopotential_height, obsZ, obsValue, obsErr, &
                                QCflags, missing)

implicit none

! Subroutine arguments
integer, intent(in)          :: nobs                     ! Number of observations
integer, intent(in)          :: n_horiz                  ! Number of geovals per observation
character(len=*), intent(in) :: rmatrix_filename         ! Name of the R-matrix file
integer, intent(in)          :: obsSatid(:)              ! Satellite identifier
integer, intent(in)          :: obsOrigC(:)              ! Originating centre number
integer, intent(in)          :: nlevs                    ! Number of model levels
real, intent(in)             :: air_temperature(:,:)     ! Temperature of the model background
real, intent(in)             :: geopotential_height(:,:) ! Geopotential height of the model levels
real(kind_real), intent(in)  :: obsZ(:)                  ! Height of the observation
real(kind_real), intent(in)  :: obsValue(:)              ! The observed value
real(kind_real), intent(out) :: obsErr(:)                ! The calculated observation error (uncertainty)
integer(c_int),  intent(in)  :: QCflags(:)               ! Quality control flags for the observations
real(kind_real), intent(in)  :: missing                  ! Missing value indicator

! Local parameters
integer, parameter :: Rmax_num = 1000                    ! Max number of R matrices to be read in

! Local variables
type(rmatrix_type), allocatable :: Rmatrix_list(:) ! List of all the R matrices to use
type(rmatrix_type) :: Rmatrix                     ! The chosen R matrix
real(kind_real) :: frac_err                       ! Fractional observation error
real :: av_temp
integer :: npoints
integer :: ilev
integer :: R_num_sats                             ! Actual number of R-matrices read in
character(len=200) :: Message                     ! Message to be output
integer :: iob                                    ! Loop variable, observation number
integer :: igeoval                                ! Loop variable, geoval number
integer :: iheight                                ! Loop variable, height in profile

! Read in R matrix data
CALL ufo_roobserror_getrmatrix(Rmax_num,         &  ! Max number of R matrices to read in
                               rmatrix_filename, &  ! The name of the file to be read in
                               Rmatrix_list,     &  ! List of all R matrices to use
                               R_num_sats)          ! Number of R matrices read in
!--------------------------------------------------------
! Choose the R-matrix values.  We use different code if we are passed
! a matrix which depends on latitude or average temperature
!--------------------------------------------------------

do iob = 1, nobs
  if (QCflags(iob) .eq. 0) then

    IF (RMatrix_list(1) % av_temp > 0) THEN
      !--------------------------------------------------------
      ! Choose R matrix depending on satid, origctr and the average
      ! background temperature between the surface and 20km
      !--------------------------------------------------------

      ! Calculate the average troposphere temperature for this profile

      av_temp = 0
      npoints = 0
      igeoval = (iob-1) * n_horiz + (n_horiz + 1) / 2
      DO ilev = 1, nlevs
        IF (geopotential_height(igeoval, ilev) < RMatrix_list(1) % max_height) THEN
          av_temp = av_temp + air_temperature(igeoval, ilev)
          npoints = npoints + 1
        END IF
      END DO

      IF (npoints > 0) THEN
        av_temp = av_temp / npoints
      ELSE
        av_temp = missing
      END IF

      ! Find the observation error matrix which best matches the average
      ! temperature we found

      CALL ufo_roobserror_interpolate_rmatrix(obsSatid(iob),   &
                                              obsOrigC(iob),   &
                                              av_temp,         &
                                              R_num_sats,      &
                                              RMatrix_list,    &
                                              RMatrix)

    ELSE
      WRITE (Message, '(2A)') "RMatrices must have positive average ", &
                                   "temperature"
      CALL abor1_ftn(Message)
    END IF

    do iheight = 1, Rmatrix % num_heights - 1
      if (obsZ(iob) < Rmatrix % height(iheight + 1)) then
        exit
      end if
    end do

    ! Fractional error
    frac_err = Rmatrix % frac_err(iheight) + &
               (Rmatrix % frac_err(iheight + 1) - Rmatrix % frac_err(iheight)) * &
               (obsZ(iob) - Rmatrix % height(iheight)) / &
               (Rmatrix % height(iheight + 1) - Rmatrix % height(iheight))

    WRITE(Message,'(A,I8,2F16.4,2E26.8)') 'Result', iob, obsZ(iob), frac_err, &
        ObsErr(iob), MAX(frac_err * obsValue(iob), Rmatrix % min_error)
    CALL fckit_log % debug(Message)

    ! Standard deviation
    ObsErr(iob) = MAX(frac_err * obsValue(iob), Rmatrix % min_error)

  else
    obsErr(iob) = missing
  end if
end do

end subroutine gnssro_obserr_avtemp


subroutine gnssro_obserr_latitude(nobs, rmatrix_filename, obsSatid, obsOrigC, obsLat, obsZ, obsValue, obsErr, QCflags, missing)

implicit none

! Subroutine arguments
integer, intent(in)          :: nobs              ! Number of observations
character(len=*), intent(in) :: rmatrix_filename  ! Name of the R-matrix file
integer, intent(in)          :: obsSatid(:)       ! Satellite identifier
integer, intent(in)          :: obsOrigC(:)       ! Originating centre number
real(kind_real), intent(in)  :: obsLat(:)         ! Latitude of the observation
real(kind_real), intent(in)  :: obsZ(:)           ! Height of the observation
real(kind_real), intent(in)  :: obsValue(:)       ! The observed value
real(kind_real), intent(out) :: obsErr(:)         ! The calculated observation error (uncertainty)
integer(c_int),  intent(in)  :: QCflags(:)        ! Quality control flags for the observations
real(kind_real), intent(in)  :: missing           ! Missing value indicator

! Local parameters
integer, parameter :: Rmax_num = 1000             ! Max number of R matrices to be read in

! Local variables
type(rmatrix_type), allocatable :: Rmatrix_list(:) ! List of all the R matrices to use
type(rmatrix_type) :: Rmatrix                     ! The chosen R matrix
real(kind_real) :: frac_err                       ! Fractional observation error
integer :: R_num_sats                             ! Actual number of R-matrices read in
character(len=200) :: Message                     ! Message to be output
integer :: iob                                    ! Loop variable, observation number
integer :: iheight                                ! Loop variable, height in profile

! Read in R matrix data
CALL ufo_roobserror_getrmatrix(Rmax_num,         &  ! Max number of R matrices to read in
                               rmatrix_filename, &  ! The name of the file to be read in
                               Rmatrix_list,     &  ! List of all R matrices to use
                               R_num_sats)          ! Number of R matrices read in

!--------------------------------------------------------
! Choose R matrix depending on satid, origctr and latitude
! No interpolation between matrices to match old code
!--------------------------------------------------------

do iob = 1, nobs
  if (QCflags(iob) .eq. 0) then
    IF (RMatrix_list(1) % latitude > 0) THEN
      CALL ufo_roobserror_findnearest_rmatrix(obsSatid(iob),  &
                                              obsOrigC(iob),  &
                                              obsLat(iob),    &
                                              R_num_sats,     &
                                              RMatrix_list,   &
                                              RMatrix)
    ELSE
      WRITE (Message, '(2A)') "RMatrices must have positive average ", &
                                   "temperature or latitude"
      CALL abor1_ftn(Message)
    END IF

    do iheight = 1, Rmatrix % num_heights
      if (obsZ(iob) < Rmatrix % height(iheight)) then
        exit
      end if
    end do

    ! Calculate fractional error
    if (iheight == 1) then
      frac_err = Rmatrix % frac_err(iheight)
    else if (iheight > Rmatrix % num_heights) then
      frac_err = Rmatrix % frac_err(RMatrix % num_heights)
    else
      frac_err = Rmatrix % frac_err(iheight - 1) + &
                 (Rmatrix % frac_err(iheight) - Rmatrix % frac_err(iheight - 1)) * &
                 (obsZ(iob) - Rmatrix % height(iheight - 1)) / &
                 (Rmatrix % height(iheight) - Rmatrix % height(iheight - 1))
    end if

    WRITE(Message,'(A,I8,2F16.4,2E21.8,F12.4)') 'Result', iob, obsZ(iob), &
        frac_err, ObsErr(iob), MAX(frac_err * obsValue(iob), Rmatrix % min_error), &
        obsLat(iob)
    CALL fckit_log % debug(Message)

    ! Standard deviation
    ObsErr(iob) = MAX (frac_err * obsValue(iob), Rmatrix % min_error)

  else
    WRITE(Message,'(A,I8,2F16.4)') 'Missing', iob, obsZ(iob), obsLat(iob)
    CALL fckit_log % debug(Message)
    obsErr(iob) = missing
  end if
end do

end subroutine gnssro_obserr_latitude

end module gnssro_mod_obserror

