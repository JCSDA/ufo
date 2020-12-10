! (C) Copyright 2020 Met Office
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module to get the jacobian for the 1D-Var

module ufo_rttovonedvarcheck_minimize_jacobian_mod

use fckit_configuration_module, only: fckit_configuration
use iso_c_binding
use kinds
use ufo_geovals_mod
use ufo_radiancerttov_mod
use ufo_rttovonedvarcheck_constants_mod
use ufo_rttovonedvarcheck_minimize_utils_mod, only: &
        ufo_rttovonedvarcheck_Qsplit
use ufo_rttovonedvarcheck_ob_mod
use ufo_rttovonedvarcheck_profindex_mod
use ufo_vars_mod

implicit none
private

public ufo_rttovonedvarcheck_get_jacobian

contains

!------------------------------------------------------------------------------
!> Get the jacobian used in the 1D-Var.
!!
!! \author Met Office
!!
!! \date 09/06/2020: Created
!!
subroutine ufo_rttovonedvarcheck_get_jacobian(geovals, ob, channels, &
                                              obsdb, conf, profindex, &
                                              prof_x, hofxdiags, rttov_simobs, &
                                              hofx, H_matrix)

implicit none

! subroutine arguments
type(ufo_geovals), intent(in)                     :: geovals       !< model data at obs location
type(ufo_rttovonedvarcheck_ob), intent(inout)     :: ob            !< satellite metadata
integer, intent(in)                               :: channels(:)   !< channels used for this calculation
type(c_ptr), value, intent(in)                    :: obsdb         !< observation database
type(fckit_configuration), intent(in)             :: conf          !< configuration
type(ufo_rttovonedvarcheck_profindex), intent(in) :: profindex     !< index array for x vector
real(kind_real), intent(in)                       :: prof_x(:)     !< x vector
type(ufo_geovals), intent(inout)                  :: hofxdiags     !< model data to pass the jacobian
type(ufo_radiancerttov), intent(inout)            :: rttov_simobs  !< rttov simulate obs object
real(kind_real), intent(out)                      :: hofx(:)       !< BT's
real(kind_real), intent(out)                      :: H_matrix(:,:) !< Jacobian

select case (trim(ob % forward_mod_name))
  case ("RTTOV")
    write(*,*) "RTTOV get H_matrix"
    call ufo_rttovonedvarcheck_GetHmatrixRTTOVsimobs(geovals, ob, obsdb, &
                                              rttov_simobs, channels, &
                                              profindex, prof_x(:), hofxdiags, &
                                              hofx(:), H_matrix) ! out

   case default
      write(*,*) "No suitable forward model exiting"
      stop
end select

end  subroutine ufo_rttovonedvarcheck_get_jacobian

!------------------------------------------------------------------------------
!> Get the jacobian from rttov and if neccessary convert 
!! to variables used in the 1D-Var.
!!
!! \details Heritage: Ops_SatRad_GetHmatrix_RTTOV12.f90
!!
!! \warning mwemiss and emisspc not implemented yet
!!
!! \author Met Office
!!
!! \date 09/06/2020: Created
!!
subroutine ufo_rttovonedvarcheck_GetHmatrixRTTOVsimobs(geovals, ob, obsdb, &
                                       rttov_data, channels, profindex, prof_x, &
                                       hofxdiags, hofx, H_matrix)

implicit none

! subroutine arguments
type(ufo_geovals), intent(in)                     :: geovals       !< model data at obs location
type(ufo_rttovonedvarcheck_ob), intent(inout)     :: ob            !< satellite metadata
type(c_ptr), value, intent(in)                    :: obsdb         !< observation database
type(ufo_radiancerttov), intent(inout)            :: rttov_data    !< structure for running rttov_k
integer, intent(in)                               :: channels(:)   !< channels used for this calculation
type(ufo_rttovonedvarcheck_profindex), intent(in) :: profindex     !< index array for x vector
real(kind_real), intent(in)                       :: prof_x(:)     !< x vector
type(ufo_geovals), intent(inout)                  :: hofxdiags     !< model data to pass the jacobian
real(kind_real), intent(out)                      :: hofx(:)       !< BT's
real(kind_real), intent(out)                      :: H_matrix(:,:) !< Jacobian

! Local arguments
integer :: nchans, nlevels, nq_levels
integer :: i, j
integer :: chan
logical :: RTTOV_GasunitConv = .false.
real(kind_real),allocatable  :: q_kgkg(:)
real(kind_real)              :: s2m_kgkg
type(ufo_geoval), pointer    :: geoval
real(kind_real), allocatable :: temperature(:)
real(kind_real), allocatable :: pressure(:)
real(kind_real), allocatable :: dq_dqt(:)
real(kind_real), allocatable :: dql_dqt(:)
real(kind_real), allocatable :: dqi_dqt(:)
real(kind_real), allocatable :: dBT_dq(:)
real(kind_real), allocatable :: dBT_dql(:)
character(len=max_string)    :: varname
real(c_double)               :: BT(size(ob % channels_all))

nchans = size(channels)

call rttov_data % simobs(geovals, obsdb, size(ob % channels_all), 1, BT, hofxdiags, ob_info=ob)

! --------------------
!Get hofx for just channels used
!--------------------
all_chan_loop: do i = 1, size(ob % channels_all)
  do j = 1, nchans
    if(channels(j) == ob % channels_all(i)) then
      hofx(j) = BT(i)
      cycle all_chan_loop
    end if
  end do
end do all_chan_loop

!----------------------------------------------------------------
! 1.1) Temperature - invert as RTTOV level 1 as top of atmosphere and
!               1Dvar profile  and is from the surface up.
!      var_ts - air_temperature
! Note : RTTOV jacobian is TOA -> surface same as prof_x
!----------------------------------------------------------------
if (profindex % t(1) > 0) then
  nlevels = profindex % t(2) - profindex % t(1) + 1
  do i = 1, nchans
    write(varname,"(3a,i0)") "brightness_temperature_jacobian_",trim(var_ts),"_",channels(i) ! K
    call ufo_geovals_get_var(hofxdiags , varname, geoval)
    H_matrix(i,profindex % t(1):profindex % t(2)) = geoval % vals(:,1)
  end do
end if

!------
! 1.2) Water vapour
!------

! Water Vapour Jacobians must be converted from
! kg/kg to ln(g/kg) - the unit conversion cancels, then:
! dy/d(ln q) = dy/dq * q(kg/kg)
! var_q    = "specific_humidity"     ! kg/kg
! Note : RTTOV jacobian is TOA -> surface same as prof_x
if (profindex % q(1) > 0) then

  nq_levels = profindex % q(2)-profindex % q(1)+1
  allocate(q_kgkg(nq_levels))

  q_kgkg(:) = exp(prof_x(profindex % q(1):profindex % q(2))) / 1000.0_kind_real

  do i = 1, nchans
    write(varname,"(3a,i0)") "brightness_temperature_jacobian_",trim(var_q),"_",channels(i) ! kg/kg
    call ufo_geovals_get_var(hofxdiags, varname, geoval)
    H_matrix(i,profindex % q(1):profindex % q(2)) = geoval % vals(:,1) * q_kgkg(:)
  end do

  deallocate(q_kgkg)

end if

!------
! 1.3) Total water
!------

! For the sake of this first implementation this will not include liquid
! and ice water content just water vapour which is consistent with the
! profile loaded from GeoVaLs.
! var_q    = "specific_humidity"     ! kg/kg
! var_clw  = "mass_content_of_cloud_liquid_water_in_atmosphere_layer"
! var_cli  = "mass_content_of_cloud_ice_in_atmosphere_layer"
! Note : RTTOV jacobian is TOA -> surface same as prof_x
if (profindex % qt(1) > 0) then

  allocate(q_kgkg(nlevels))
  allocate(temperature(nlevels))
  allocate(pressure(nlevels))
  allocate(dq_dqt(nlevels))
  allocate(dql_dqt(nlevels))
  allocate(dqi_dqt(nlevels))
  allocate(dBT_dq(nlevels))
  allocate(dBT_dql(nlevels))

  ! Get qtotal from profile
  q_kgkg(:) = exp(prof_x(profindex % qt(1):profindex % qt(2))) / 1000.0_kind_real

  ! Get temperature and pressure from geovals
  ! var_ts   = "air_temperature" K
  call ufo_geovals_get_var(geovals, trim(var_ts), geoval)
  temperature(:) = geoval%vals(nlevels:1:-1, 1)
  ! var_prs  = "air_pressure" Pa
  call ufo_geovals_get_var(geovals, trim(var_prs), geoval)
  pressure(:) = geoval%vals(nlevels:1:-1, 1)

  ! Calculate the gradients with respect to qtotal
  call ufo_rttovonedvarcheck_Qsplit (0,    & ! in
                          temperature,     & ! in
                          pressure,        & ! in
                          nlevels,         & ! in
                          q_kgkg(:),       & ! in
                          dq_dqt(:),       & ! out
                          dql_dqt(:),      & ! out
                          dqi_dqt(:))        ! out

  ! Calculate jacobian wrt humidity and clw
  do i = 1, nchans

    write(varname,"(3a,i0)") "brightness_temperature_jacobian_", trim(var_q), "_", channels(i) ! kg/kg
    call ufo_geovals_get_var(hofxdiags, varname, geoval)
    dBT_dq(:) = geoval % vals(:,1)

    write(varname,"(3a,i0)") "brightness_temperature_jacobian_", trim(var_clw), "_", channels(i) ! kg/kg
    call ufo_geovals_get_var(hofxdiags, varname, geoval)
    dBT_dql(:) = geoval % vals(:,1)

    H_matrix(i,profindex % qt(1):profindex % qt(2)) = &
            (dBT_dq(:)  * dq_dqt(:) + &
             dBT_dql(:) * dql_dqt(:) ) * q_kgkg(:)
  end do

  ! Clean up
  deallocate(q_kgkg)
  deallocate(temperature)
  deallocate(pressure)
  deallocate(dq_dqt)
  deallocate(dql_dqt)
  deallocate(dqi_dqt)
  deallocate(dBT_dq)
  deallocate(dBT_dql)

end if

!----
! 2.) Single-valued variables
!----

! 2.1) Surface Temperature - var_sfc_t2m = "surface_temperature"

if (profindex % t2 > 0) then
  do i = 1, nchans
    write(varname,"(3a,i0)") "brightness_temperature_jacobian_",trim(var_sfc_t2m),"_",channels(i) ! K
    call ufo_geovals_get_var(hofxdiags, varname, geoval)
    H_matrix(i,profindex % t2) = geoval % vals(1,1)
  end do
end if

! 2.2) Water vapour - var_sfc_q2m = "specific_humidity_at_two_meters_above_surface" ! (kg/kg)

if (profindex % q2 > 0) then
  s2m_kgkg = exp(prof_x(profindex % q2)) / 1000.0_kind_real
  do i = 1, nchans
    write(varname,"(3a,i0)") "brightness_temperature_jacobian_",trim(var_sfc_q2m),"_",channels(i) ! kg/kg
    call ufo_geovals_get_var(hofxdiags, varname, geoval)
    H_matrix(i,profindex % q2) = geoval % vals(1,1) * s2m_kgkg
  end do
end if

! 2.3) Surface pressure - var_sfc_p2m = "air_pressure_at_two_meters_above_surface" ! (Pa)

if (profindex % pstar > 0) then
  do i = 1, nchans
    write(varname,"(3a,i0)") "brightness_temperature_jacobian_",trim(var_sfc_p2m),"_",channels(i)
    call ufo_geovals_get_var(hofxdiags, varname, geoval)
    H_matrix(i,profindex % pstar) = geoval % vals(1,1)
  end do
end if

! This has been left in for future development
! 2.4) Windspeed - var_u = "eastward_wind"
! Remember that all wind has been transferred to u and v is set to zero for
! windspeed retrieval.

!if (profindex % windspeed > 0) then
!  do i = 1, nchans
!    write(varname,"(3a,i0)") "brightness_temperature_jacobian_",trim(var_u),"_",channels(i)
!    call ufo_geovals_get_var(hofxdiags, varname, geoval)
!    H_matrix(i,profindex % windspeed) = geoval % vals(1,1)
!  end do
!end if

! 2.5) Skin temperature - var_sfc_tskin = "skin_temperature"  ! (K)

if (profindex % tstar > 0) then
  do i = 1, nchans
    write(varname,"(3a,i0)") "brightness_temperature_jacobian_",trim(var_sfc_tskin),"_",channels(i)
    call ufo_geovals_get_var(hofxdiags, varname, geoval)
    H_matrix(i,profindex % tstar) = geoval % vals(1,1)
  end do
end if

! This has been left in for future development
! 2.5) Cloud top pressure
! This is not in rttov interface yet
!if (profindex % cloudtopp > 0) then
!  do i = 1, nchans
!    varname = "cloud_top_pressure"
!    write(varname,"(3a,i0)") "brightness_temperature_jacobian_",trim(varname),"_",channels(i)
!    call ufo_geovals_get_var(hofxdiags, varname, geoval)
!    H_matrix(i,profindex % cloudtopp) = geoval % vals(1,1)
!  end do
!end if

! This has been left in for future development
! 2.6) Effective cloud fraction
! This is not in rttov interface yet
!if (profindex % cloudfrac > 0) then
!  do i = 1, nchans
!    varname = "cloud_fraction"
!    call ufo_geovals_get_var(hofxdiags, varname, geoval)
!    H_matrix(i,profindex % cloudfrac) = geoval % vals(1,1)
!  end do
!end if

!----
! 3.) Emissivities
! This has been left in for future development
!----

! 3.1 Microwave Emissivity - var_sfc_emiss = "surface_emissivity"

!if (profindex % mwemiss(1) > 0) then
!  ! The emissivity matrix needs to be "unpacked" as it is only one
!  ! dimensional over channels - implying you want to retrieve a
!  ! single emissivity value. It is unpacked here so that each emissivity
!  ! has a corresponding entry for the relevant channel. Note that this is
!  ! a bit physically dubious as several channels have the same frequency, etc.
!  ! This complexity is dealt with in the B Matrix.
!  ! Check that we want only the diagonal elements to be non-zero
!  do j = 1, size(EmissMap)
!      chan = EmissMap(j)
!    do i = 1, nchans
!      if (channels(i) == chan) then
!        write(varname,"(3a,i0)") "brightness_temperature_jacobian_",trim(var_sfc_emiss),"_",channels(i)
!        call ufo_geovals_get_var(hofxdiags, varname, geoval)
!        H_matrix(i,profindex % mwemiss(1) + j - 1) = geoval % vals(1,1)
!      end if
!    end do
!  end do
!end if

!! 3.2. Infrared Emissivity - work in progress
!
!IF (profindex % emisspc(1) > 0) THEN
!
!  DO i = 1, nchans
!    emissivity_K(i) = rttov_data % profiles_k(i) % emissivity(1)
!  END DO
!
!  CALL Ops_SatRad_EmisKToPC (nchans,                                                           & ! in
!                             Channels,                                                         & ! in
!                             nemisspc,                                                         & ! in
!                             emiss(:),                                                         & ! in
!                             emissivity_K(:),                                                  & ! in
!                             H_matrix(1:nchans,profindex % emisspc(1):profindex % emisspc(2)))   ! out
!END IF
!
! Here for diagnostics
!call ufo_rttovonedvarcheck_PrintHmatrix( &
!  nchans,   &           ! in
!  size(prof_x),  &      ! in
!  ob % channels_used, & ! in
!  H_matrix, &           ! in
!  profindex )           ! in

end subroutine ufo_rttovonedvarcheck_GetHmatrixRTTOVsimobs

!---------------------------------------------------------------------------
!> Routine to print the contents of the jacobian for testing
!!
!! \details Heritage: Ops_SatRad_PrintHMatrix.f90
!!
!! \author Met Office
!!
!! \date 09/06/2020: Created
!!
subroutine ufo_rttovonedvarcheck_PrintHmatrix( &
  nchans,   &       ! in
  nprofelements, &  ! in
  channels, &       ! in
  H_matrix, &       ! in
  profindex )       ! in

implicit none

!Subroutine arguments:
integer, intent(in)              :: nchans
integer, intent(in)              :: nprofelements
integer, intent(in)              :: channels(nchans)
real(kind_real), intent(in)      :: H_matrix(nchans,nprofelements)
type(ufo_rttovonedvarcheck_profindex), intent(in) :: profindex

! Local variables:
integer :: i
character(len=10) :: int_fmt
character(len=12) :: real_fmt
character(len=3) :: txt_nchans
character(len=*), parameter :: RoutineName = "ufo_rttovonedvarcheck_PrintHmatrix"
!-------------------------------------------------------------------------------

write( unit=txt_nchans,fmt='(i3)' )  nchans
write( unit=int_fmt,fmt='(a)' ) '(' // trim(txt_nchans) // 'I13)'
write( unit=real_fmt,fmt='(a)' ) '(' // trim(txt_nchans) // 'E13.5)'

write(*,*)

write(*, int_fmt) channels(:)

  if ( profindex % t(1) > 0 ) THEN
    write(*, '(a)') 'Temperature Profile'
    do i = profindex%t(1),profindex%t(2)
      write(*, real_fmt)  H_matrix(:,i)
    end do
  end if

  if ( profindex % q(1) > 0 ) THEN
    write(*, '(a)') 'q Profile'
    do i = profindex%q(1),profindex%q(2)
      write(*, real_fmt)  H_matrix(:,i)
    end do
  end if

  if ( profindex % qt(1) > 0 ) THEN
    write(*, '(a)') 'qt Profile /1000'
    do i = profindex%qt(1),profindex%qt(2)
      write(*, real_fmt)  H_matrix(:,i)/1000
    end do
  end if

  if ( profindex % o3profile(1) > 0 ) THEN
    write(*, '(a)') 'Ozone Profile'
    do i = profindex%o3profile(1),profindex%o3profile(2)
      write(*, real_fmt)  H_matrix(:,i)
    end do
  end if

  if ( profindex % o3total > 0 ) THEN
    write(*, '(a)') 'Total Column Ozone'
    write(*, real_fmt)  H_matrix(:,profindex%o3total)
  end if


  if ( profindex % lwp > 0 ) THEN
    write(*, '(a)') 'LWP'
    write(*, real_fmt)  H_matrix(i,profindex % lwp)
  end if

  if ( profindex % t2 > 0 ) THEN
    write(*, '(a)') '2m T'
    write(*, real_fmt)  H_matrix(:,profindex % t2)
  end if

  if ( profindex % q2 > 0 ) THEN
    write(*, '(a)') '2m q'
    write(*, real_fmt)  H_matrix(:,profindex % q2)
  end if

  if ( profindex % pstar > 0 ) THEN
    write(*, '(a)') 'P Star'
    write(*, real_fmt)  H_matrix(:,profindex % pstar)
  end if

  if ( profindex % windspeed > 0 ) THEN
    write(*, '(a)') 'Windspeed'
    write(*, real_fmt)  H_matrix(:,profindex % windspeed)
  end if

  if ( profindex % tstar > 0 ) THEN
    write(*, '(a)') 'Skin Temperature'
    write(*, real_fmt)  H_matrix(:,profindex % tstar)
  end if

  if ( profindex % mwemiss(1) > 0) THEN
    write(*, '(a)') 'Microwave emissivity retrieval'
    do i = profindex%mwemiss(1),profindex%mwemiss(2)
      write(*, real_fmt)  H_matrix(:,i)
    end do
  end if

  if ( profindex % cloudtopp > 0 ) THEN
    write(*, '(a)') 'Cloud top pressure'
    write(*, real_fmt)  H_matrix(:,profindex % cloudtopp)
  end if

  if ( profindex % cloudfrac > 0 ) THEN
    write(*, '(a)') 'Effective cloud fraction'
    write(*, real_fmt)  H_matrix(:,profindex % cloudfrac)
  end if

write(*,*)
write(*, '(a)') 'End H-Matrix'
write(*, '(a)') '------------------------'

end  subroutine ufo_rttovonedvarcheck_PrintHmatrix

! ------------------------------------------------------------

end module ufo_rttovonedvarcheck_minimize_jacobian_mod
