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
public :: bending_angle_obserr_NBAM, refractivity_obserr_NCEP
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
!  processed at EUMETSAT
   if( (ObsSaid(i)==41).or.(ObsSaid(i)==722).or.(ObsSaid(i)==723).or.   &
       (ObsSaid(i)>=3.and.ObsSaid(i)<=5).or.(ObsSaid(i)==42).or.        &
       (ObsSaid(i)==821).or.(ObsSaid(i)==421).or.    &
       (ObsSaid(i)==440).or.(ObsSaid(i)==43) ) then
       if( abs(obsLat(i)) > 40.0 ) then
         if(H_km > 12.0) then
           obsErr(i)=0.19032 +0.287535 *H_km-0.00260813*H_km**2
         else
           obsErr(i)=-3.20978 +1.26964 *H_km-0.0622538 *H_km**2
         endif
       else
         if(H_km > 18.0) then
           obsErr(i)=-1.87788 +0.354718 *H_km-0.00313189 *H_km**2
         else
           obsErr(i)=-2.41024 +0.806594 *H_km-0.027257 *H_km**2
         endif
       endif
!  COSMIC-2 (750-755) and Commercial (265-269)
   else if ( (ObsSaid(i) >= 750 .and. ObsSaid(i) <= 755) .or. (ObsSaid(i) >= 265 .and. ObsSaid(i) <= 269) ) then
       if ( abs(obsLat(i)) > 40.0 ) then
         if (H_km <= 8.0) then
            obsErr(i) = -1.0304261+0.3203316*H_km+0.0141337*H_km**2
         else if (H_km > 8.0.and.H_km <= 12.0) then
            obsErr(i) = 2.1750271+0.0431177*H_km-0.0008567*H_km**2
         else
            obsErr(i) = -0.3447429+0.2829981*H_km-0.0028545*H_km**2
         endif
       else
         if (H_km <= 4.0) then
            obsErr(i) = 0.7285212-1.1138755*H_km+0.2311123*H_km**2
         elseif (H_km <= 18.0 .and. H_km > 4.0) then
            obsErr(i) = -3.3878629+0.8691249*H_km-0.0297196*H_km**2
         else
            obsErr(i) = -2.3875749+0.3667211*H_km-0.0037542*H_km**2
         endif
       endif
   else ! Other
       if( abs(obsLat(i)) > 40.0 ) then
         if ( H_km > 12.00 ) then
            obsErr(i)=-0.685627 +0.377174 *H_km-0.00421934 *H_km**2
         else
            obsErr(i)=-3.27737 +1.20003 *H_km-0.0558024 *H_km**2
         endif
       else
         if( H_km > 18.0 ) then
            obsErr(i)=-2.73867 +0.447663 *H_km-0.00475603 *H_km**2
         else
            obsErr(i)=-3.45303 +0.908216 *H_km-0.0293331 *H_km**2
         endif
       endif
   endif
   obsErr(i) = 0.001 /abs(exp(obsErr(i)))

end if

end do

end subroutine bending_angle_obserr_NBAM
!---------------------------------------

subroutine refractivity_obserr_NCEP(obsLat, obsZ, nobs, obsErr, QCflags,missing)
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

end subroutine refractivity_obserr_NCEP


subroutine gnssro_obserr_avtemp(nobs, n_horiz, rmatrix_filename, obsSatid, obsOrigC, obsZ, &
                                obsValue, averageTemp, obsErr, &
                                QCflags, missing, allow_extrapolation, record_number, &
                                sort_order, unique, verboseOutput)

implicit none

! Subroutine arguments
integer, intent(in)              :: nobs                     ! Number of observations
integer, intent(in)              :: n_horiz                  ! Number of geovals per observation
character(len=*), intent(in)     :: rmatrix_filename         ! Name of the R-matrix file
integer, intent(in)              :: obsSatid(:)              ! Satellite identifier
integer, intent(in)              :: obsOrigC(:)              ! Originating centre number
real(kind_real), intent(in)      :: averageTemp(:)           ! Average temperature below a given height
real(kind_real), intent(in)      :: obsZ(:)                  ! Height of the observation
real(kind_real), intent(in)      :: obsValue(:)              ! The observed value
real(kind_real), intent(out)     :: obsErr(:)                ! The calculated observation error (uncertainty)
integer(c_int),  intent(in)      :: QCflags(:)               ! Quality control flags for the observations
real(kind_real), intent(in)      :: missing                  ! Missing value indicator
logical, intent(in)              :: allow_extrapolation      ! Extrapolate the obs errors outside the range specified
integer(c_size_t), intent(in)    :: record_number(1:nobs)    ! Number used to identify unique profiles in the data
integer, intent(in)              :: sort_order(1:nobs)       ! An index to sort the record numbers
integer, allocatable, intent(in) :: unique(:)                ! Set of unique profile numbers
logical, intent(in)              :: verboseOutput            ! Whether to output extra debugging information

! Local parameters
integer, parameter :: Rmax_num = 1000                    ! Max number of R matrices to be read in

! Local variables
type(rmatrix_type), allocatable :: Rmatrix_list(:) ! List of all the R matrices to use
type(rmatrix_type) :: Rmatrix                      ! The chosen R matrix
real(kind_real) :: frac_err                        ! Fractional observation error
integer :: R_num_sats                              ! Actual number of R-matrices read in
character(len=200) :: Message                      ! Message to be output
integer :: iob                                     ! Loop variable, observation number
integer :: iPoint                                  ! Loop variable, point in the profile
integer :: igeoval                                 ! Loop variable, geoval number
integer :: iheight                                 ! Loop variable, height in profile
integer :: start_point                             ! Starting index of the current profile
integer :: current_point                           ! Ending index of the current profile
integer :: iprofile                                ! Loop variable, profile number

! Read in R matrix data
CALL ufo_roobserror_getrmatrix(Rmax_num,         &  ! Max number of R matrices to read in
                               rmatrix_filename, &  ! The name of the file to be read in
                               Rmatrix_list,     &  ! List of all R matrices to use
                               R_num_sats)          ! Number of R matrices read in
!--------------------------------------------------------
! Choose the R-matrix values.  We use different code if we are passed
! a matrix which depends on latitude or average temperature
!--------------------------------------------------------

! For every profile that we have found, perform a 1DVar minimisation
current_point = 1
do iprofile = 1, size(unique)
  start_point = current_point
  ! Work out which observations belong to the current profile
  do current_point = start_point, nobs
    if (unique(iprofile) /= record_number(sort_order(current_point))) exit
  end do

  do iPoint = start_point, current_point-1
    iob = sort_order(iPoint)
    if (QCflags(iob) .eq. 0) then
      if (RMatrix_list(1) % av_temp > 0) then
        !--------------------------------------------------------
        ! Choose R matrix depending on satid, origctr and the average
        ! background temperature between the surface and 20km
        !--------------------------------------------------------

        CALL ufo_roobserror_interpolate_rmatrix(obsSatid(iob),   &
                                                obsOrigC(iob),   &
                                                averageTemp(iob),&
                                                R_num_sats,      &
                                                RMatrix_list,    &
                                                RMatrix)

      ELSE
        WRITE (Message, '(2A)') "RMatrices must have positive average ", &
                                     "temperature"
        CALL abor1_ftn(Message)
      END IF

      do iheight = 1, Rmatrix % num_heights
        if (obsZ(iob) < Rmatrix % height(iheight)) then
          exit
        end if
      end do

      ! If we are allowing the uncertainties to be extrapolated, then reset the
      ! value of iheight if it is out of range.
      if (allow_extrapolation) then
        if (iheight == 1) then
          iheight = 2
        else if (iheight > Rmatrix % num_heights) then
          iheight = Rmatrix % num_heights
        end if
      end if

      ! Calculate fractional error
      if (iheight == 1) then
        frac_err = Rmatrix % frac_err(1)
      else if (iheight > Rmatrix % num_heights) then
        frac_err = Rmatrix % frac_err(RMatrix % num_heights)
      else
        frac_err = Rmatrix % frac_err(iheight - 1) + &
                   (Rmatrix % frac_err(iheight) - Rmatrix % frac_err(iheight - 1)) * &
                   (obsZ(iob) - Rmatrix % height(iheight - 1)) / &
                   (Rmatrix % height(iheight) - Rmatrix % height(iheight - 1))
      end if

      if (verboseOutput) then
        WRITE(Message,'(A,I8,2F16.4,3E26.8)') 'Result', iob, obsZ(iob), frac_err, &
            ObsErr(iob), MAX(frac_err * obsValue(iob), Rmatrix % min_error), averageTemp(iob)
        CALL fckit_log % info(Message)
      end if

      ! Standard deviation
      ObsErr(iob) = MAX(frac_err * obsValue(iob), Rmatrix % min_error)

    else
      obsErr(iob) = missing
    end if
  end do
end do

end subroutine gnssro_obserr_avtemp


subroutine gnssro_obserr_latitude(nobs, rmatrix_filename, obsSatid, obsOrigC, &
    obsLat, obsZ, obsValue, obsErr, QCflags, missing, allow_extrapolation, &
    record_number, sort_order, unique, verboseOutput)

use missing_values_mod

implicit none

! Subroutine arguments
integer, intent(in)              :: nobs                 ! Number of observations
character(len=*), intent(in)     :: rmatrix_filename     ! Name of the R-matrix file
integer, intent(in)              :: obsSatid(:)          ! Satellite identifier
integer, intent(in)              :: obsOrigC(:)          ! Originating centre number
real(kind_real), intent(in)      :: obsLat(:)            ! Latitude of the observation
real(kind_real), intent(in)      :: obsZ(:)              ! Height of the observation
real(kind_real), intent(in)      :: obsValue(:)          ! The observed value
real(kind_real), intent(out)     :: obsErr(:)            ! The calculated observation error (uncertainty)
integer(c_int),  intent(in)      :: QCflags(:)           ! Quality control flags for the observations
real(kind_real), intent(in)      :: missing              ! Missing value indicator
logical, intent(in)              :: allow_extrapolation  ! Extrapolate errors outside of the given range
integer(c_size_t), intent(in)    :: record_number(1:nobs)! Number used to identify unique profiles in the data
integer, intent(in)              :: sort_order(1:nobs)   ! An index to sort the record numbers
integer, allocatable, intent(in) :: unique(:)            ! Set of unique profile numbers
logical, intent(in)              :: verboseOutput        ! Whether to output extra debugging information

! Local parameters
integer, parameter :: Rmax_num = 1000             ! Max number of R matrices to be read in

! Local variables
type(rmatrix_type), allocatable :: Rmatrix_list(:) ! List of all the R matrices to use
type(rmatrix_type) :: Rmatrix                      ! The chosen R matrix
real(kind_real) :: frac_err                        ! Fractional observation error
integer :: R_num_sats                              ! Actual number of R-matrices read in
character(len=200) :: Message                      ! Message to be output
integer :: iob                                     ! Loop variable, observation number
integer :: iPoint                                  ! Loop variable, point in the profile
integer :: iheight                                 ! Loop variable, height in profile
integer :: start_point                             ! Starting index of the current profile
integer :: current_point                           ! Ending index of the current profile
integer :: iprofile                                ! Loop variable, profile number

! Read in R matrix data
CALL ufo_roobserror_getrmatrix(Rmax_num,         &  ! Max number of R matrices to read in
                               rmatrix_filename, &  ! The name of the file to be read in
                               Rmatrix_list,     &  ! List of all R matrices to use
                               R_num_sats)          ! Number of R matrices read in

!--------------------------------------------------------
! Choose R matrix depending on satid, origctr and latitude
! No interpolation between matrices to match old code
!--------------------------------------------------------

! For every profile that we have found, perform a 1DVar minimisation
current_point = 1
do iprofile = 1, size(unique)
  start_point = current_point
  ! Work out which observations belong to the current profile
  do current_point = start_point, nobs
    if (unique(iprofile) /= record_number(sort_order(current_point))) exit
  end do

  do iPoint = start_point, current_point-1
    iOb = sort_order(iPoint)
    if (QCflags(iob) .eq. 0) then
      IF (RMatrix_list(1) % latitude /= &
          missing_value(RMatrix_list(1) % latitude)) THEN
        ! Use the R-matrix from the first observation in the profile
        call ufo_roobserror_findnearest_rmatrix(obsSatid(iob),       &
                                                obsOrigC(iob),       &
                                                obsLat(sort_order(start_point)), &
                                                R_num_sats,          &
                                                RMatrix_list,        &
                                                RMatrix)
      ELSE
        WRITE (Message, '(2A)') "RMatrices must have a valid latitude set"
        CALL abor1_ftn(Message)
      END IF

      do iheight = 1, Rmatrix % num_heights
        if (obsZ(iob) < Rmatrix % height(iheight)) then
          exit
        end if
      end do

      ! If we are allowing the uncertainties to be extrapolated, then reset the
      ! value of iheight if it is out of range.
      if (allow_extrapolation) then
        if (iheight == 1) then
          iheight = 2
        else if (iheight > Rmatrix % num_heights) then
          iheight = Rmatrix % num_heights
        end if
      end if

      ! Calculate fractional error
      if (iheight == 1) then
        frac_err = Rmatrix % frac_err(1)
      else if (iheight > Rmatrix % num_heights) then
        frac_err = Rmatrix % frac_err(RMatrix % num_heights)
      else
        frac_err = Rmatrix % frac_err(iheight - 1) + &
                   (Rmatrix % frac_err(iheight) - Rmatrix % frac_err(iheight - 1)) * &
                   (obsZ(iob) - Rmatrix % height(iheight - 1)) / &
                   (Rmatrix % height(iheight) - Rmatrix % height(iheight - 1))
      end if

      if (verboseOutput) then
        WRITE(Message,'(A,I8,2F16.4,2E21.8,F12.4)') 'Result', iob, obsZ(iob), &
            frac_err, ObsErr(iob), MAX(frac_err * obsValue(iob), Rmatrix % min_error), &
            obsLat(iob)
        CALL fckit_log % info(Message)
      end if

      ! Standard deviation
      ObsErr(iob) = MAX (frac_err * obsValue(iob), Rmatrix % min_error)

    else
      WRITE(Message,'(A,I8,2F16.4)') 'Missing', iob, obsZ(iob), obsLat(iob)
      obsErr(iob) = missing
    end if
  end do
end do

end subroutine gnssro_obserr_latitude

end module gnssro_mod_obserror

