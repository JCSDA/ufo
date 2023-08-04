! (C) Copyright 2020 Met Office
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module to get the jacobian for the 1D-Var

module ufo_rttovonedvarcheck_minimize_jacobian_mod

use iso_c_binding
use kinds
use ufo_constants_mod, only: zero
use ufo_geovals_mod
use ufo_radiancerttov_mod
use ufo_rttovonedvarcheck_constants_mod
use ufo_rttovonedvarcheck_ob_mod
use ufo_rttovonedvarcheck_profindex_mod
use ufo_rttovonedvarcheck_setup_mod, only: ufo_rttovonedvarcheck
use ufo_rttovonedvarcheck_utils_mod, only: ufo_rttovonedvarcheck_all_to_subset_by_channels, &
                                           ufo_rttovonedvarcheck_geovals_index_by_channels
use ufo_vars_mod
use ufo_utils_mod, only: Ops_SatRad_Qsplit

implicit none
private

public ufo_rttovonedvarcheck_get_jacobian
public ufo_rttovonedvarcheck_get_bts

contains

!------------------------------------------------------------------------------
!> Get the jacobian used in the 1D-Var.
!!
!! \author Met Office
!!
!! \date 09/06/2020: Created
!!
subroutine ufo_rttovonedvarcheck_get_jacobian(config, geovals, ob, channels, &
                                              profindex, &
                                              prof_x, hofxdiags, rttov_simobs, &
                                              hofx, H_matrix)

implicit none

! subroutine arguments
type(ufo_rttovonedvarcheck), intent(in)           :: config        !< configuration from main 1D-Var object
type(ufo_geovals), intent(in)                     :: geovals       !< model data at obs location
type(ufo_rttovonedvarcheck_ob), intent(inout)     :: ob            !< satellite metadata
integer, intent(in)                               :: channels(:)   !< channels used for this calculation
type(ufo_rttovonedvarcheck_profindex), intent(in) :: profindex     !< index array for x vector
real(kind_real), intent(in)                       :: prof_x(:)     !< x vector
type(ufo_geovals), intent(inout)                  :: hofxdiags     !< model data to pass the jacobian
type(ufo_radiancerttov), intent(inout)            :: rttov_simobs  !< rttov simulate obs object
real(kind_real), intent(out)                      :: hofx(:)       !< BT's
real(kind_real), intent(out)                      :: H_matrix(:,:) !< Jacobian

select case (trim(ob % forward_mod_name))
  case ("RTTOV")
    call ufo_rttovonedvarcheck_GetHmatrixRTTOVsimobs(geovals, config, ob, &
                                              rttov_simobs, channels, &
                                              profindex, hofxdiags, &
                                              hofx(:), H_matrix) ! out

  case default
    call abor1_ftn("rttovonedvarcheck get jacobian: no suitable forward model => exiting")
end select

end  subroutine ufo_rttovonedvarcheck_get_jacobian

!------------------------------------------------------------------------------
!> Get the BTs only.  This is much faster than running the k code
!!
!! \author Met Office
!!
!! \date 21/01/2021: Created
!!
subroutine ufo_rttovonedvarcheck_get_bts(config, geovals, ob, channels, &
                                         rttov_simobs, hofx)

implicit none

! subroutine arguments
type(ufo_rttovonedvarcheck), intent(in)           :: config        !< configuration from main 1D-Var object
type(ufo_geovals), intent(in)                     :: geovals       !< model data at obs location
type(ufo_rttovonedvarcheck_ob), intent(inout)     :: ob            !< satellite metadata
integer, intent(in)                               :: channels(:)   !< channels used for this calculation
type(ufo_radiancerttov), intent(inout)            :: rttov_simobs  !< rttov simulate obs object
real(kind_real), intent(out)                      :: hofx(:)       !< BTs


integer           :: i, j !< counters
type(ufo_geovals) :: empty_hofxdiags  !< model data to pass the jacobian
real(c_double)    :: BT(size(ob % channels_all)) !< BTs produced for all channels

select case (trim(ob % forward_mod_name))
  case ("RTTOV")
    call rttov_simobs % simobs(geovals, config % obsdb, size(ob % channels_all), 1, BT, empty_hofxdiags, ob_info=ob)
    if (ob % rterror) return
    call ufo_rttovonedvarcheck_all_to_subset_by_channels(ob % channels_all, BT, channels, hofx)
  case default
    call abor1_ftn("rttovonedvarcheck get jacobian: no suitable forward model => exiting")
end select

end  subroutine ufo_rttovonedvarcheck_get_bts

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
subroutine ufo_rttovonedvarcheck_GetHmatrixRTTOVsimobs(geovals, config, ob, &
                                       rttov_data, channels, profindex, &
                                       hofxdiags, hofx, H_matrix)

implicit none

! subroutine arguments
type(ufo_geovals), intent(in)                     :: geovals        !< model data at obs location
type(ufo_rttovonedvarcheck), intent(in)           :: config         !< configuration information
type(ufo_rttovonedvarcheck_ob), intent(inout)     :: ob             !< satellite metadata
type(ufo_radiancerttov), intent(inout)            :: rttov_data     !< structure for running rttov_k
integer, intent(in)                               :: channels(:)    !< channels used for this calculation
type(ufo_rttovonedvarcheck_profindex), intent(in) :: profindex      !< index array for x vector
type(ufo_geovals), intent(inout)                  :: hofxdiags      !< model data to pass the jacobian
real(kind_real), intent(out)                      :: hofx(:)        !< BT's
real(kind_real), intent(out)                      :: H_matrix(:,:)  !< Jacobian - (chans x profindex)

! Local arguments
integer :: nchans, nlevels, nq_levels
integer :: i, j
integer :: chan
integer :: nemisspc
integer, allocatable :: hofxdiag_index_1(:), hofxdiag_index_2(:), hofxdiag_index_3(:)
real(kind_real),allocatable  :: q_kgkg(:)
real(kind_real)              :: s2m_kgkg
type(ufo_geoval), pointer    :: geoval
real(kind_real), allocatable :: pressure(:)
real(kind_real), allocatable :: dq_dqt(:)
real(kind_real), allocatable :: dql_dqt(:)
real(kind_real), allocatable :: dqi_dqt(:)
real(kind_real), allocatable :: dBT_dq(:)
real(kind_real), allocatable :: dBT_dql(:)
real(kind_real), allocatable :: dBT_dqi(:)
real(kind_real), allocatable :: emissivity_k(:)
real(kind_real), allocatable :: emissivity(:)
character(len=MAXVARLEN)     :: varname, basename
real(c_double), allocatable  :: BT(:)
real(kind_real)              :: u, v, dBT_du, dBT_dv, windsp

! Setup varibales
nchans = size(channels)
H_matrix(:,:) = zero
allocate(BT(size(ob % channels_all)))
allocate(hofxdiag_index_1(nchans))

! Run rttov and check it works any failures return here
call rttov_data % simobs(geovals, config % obsdb, size(ob % channels_all), 1, BT, hofxdiags, ob_info=ob)
if (ob % rterror) return

! -------------------------------
!Get hofx for just channels used
!--------------------------------
call ufo_rttovonedvarcheck_all_to_subset_by_channels(ob % channels_all, BT, channels, hofx)

!----------------------------------------------------------------
! 1.1) Temperature - invert as RTTOV level 1 as top of atmosphere and
!               1Dvar profile  and is from the surface up.
!      var_ts - air_temperature
! Note : RTTOV jacobian is TOA -> surface same as prof_x
!----------------------------------------------------------------
if (profindex % t(1) > 0) then
  nlevels = profindex % t(2) - profindex % t(1) + 1

  ! Get list of indexes in hofxdiags
  write(basename,"(3a)") "brightness_temperature_jacobian_", trim(var_ts), "_"
  call ufo_rttovonedvarcheck_geovals_index_by_channels(channels, trim(basename), &
                                            hofxdiags % variables, hofxdiag_index_1)

  do i = 1, nchans
    call ufo_geovals_get_var_by_index(hofxdiags, hofxdiag_index_1(i), geoval)
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

  ! Get humidity data from geovals
  q_kgkg(:) = zero
  call ufo_geovals_get_var(geovals, var_q, geoval)
  q_kgkg(:) = geoval%vals(:, 1)

  ! Get list of indexes in hofxdiags
  write(basename,"(3a)") "brightness_temperature_jacobian_", trim(var_q), "_"
  call ufo_rttovonedvarcheck_geovals_index_by_channels(channels, trim(basename), &
                                          hofxdiags % variables, hofxdiag_index_1)

  do i = 1, nchans
    call ufo_geovals_get_var_by_index(hofxdiags, hofxdiag_index_1(i), geoval)
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
  allocate(pressure(nlevels))
  allocate(dq_dqt(nlevels))
  allocate(dql_dqt(nlevels))
  allocate(dqi_dqt(nlevels))
  allocate(dBT_dq(nlevels))
  allocate(dBT_dql(nlevels))
  allocate(dBT_dqi(nlevels))

  ! Get humidity data from geovals
  q_kgkg(:) = zero
  call ufo_geovals_get_var(geovals, var_q, geoval)
  q_kgkg(:) = q_kgkg(:) + geoval%vals(:, 1)
  call ufo_geovals_get_var(geovals, var_clw, geoval)
  q_kgkg(:) = q_kgkg(:) + geoval%vals(:, 1)
  call ufo_geovals_get_var(geovals, var_cli, geoval)
  q_kgkg(:) = q_kgkg(:) + geoval%vals(:, 1)

  ! var_prs  = "air_pressure" Pa
  call ufo_geovals_get_var(geovals, trim(var_prs), geoval)
  pressure(:) = geoval%vals(:, 1)

  ! Calculate the gradients with respect to qtotal
  call Ops_SatRad_Qsplit ( 0,      &
                    pressure(:),   &
                    ob % background_T(:),   &
                    q_kgkg(:),     & ! in
                    dq_dqt(:),     & ! out
                    dql_dqt(:),    & ! out
                    dqi_dqt(:),    & ! out
                    config % UseQtsplitRain)

  ! Get geovals indexes
  write(basename,"(3a)") "brightness_temperature_jacobian_", trim(var_q), "_"  ! kg/kg
  call ufo_rttovonedvarcheck_geovals_index_by_channels(channels, trim(basename), &
                                      hofxdiags % variables, hofxdiag_index_1)

  allocate(hofxdiag_index_2(nchans))
  write(basename,"(3a)") "brightness_temperature_jacobian_", trim(var_clw), "_" ! kg/kg
  call ufo_rttovonedvarcheck_geovals_index_by_channels(channels, trim(basename), &
                                      hofxdiags % variables, hofxdiag_index_2)

  if (config % RTTOV_mwscattSwitch) then
    allocate(hofxdiag_index_3(nchans))
    write(basename,"(3a)") "brightness_temperature_jacobian_", trim(var_cli), "_" ! kg/kg
    call ufo_rttovonedvarcheck_geovals_index_by_channels(channels, trim(basename), &
                                      hofxdiags % variables, hofxdiag_index_3)
  end if

  ! Calculate jacobian wrt humidity and clw
  do i = 1, nchans
    ! var_q
    call ufo_geovals_get_var_by_index(hofxdiags, hofxdiag_index_1(i), geoval)
    dBT_dq(:) = zero
    dBT_dq(:) = geoval % vals(:,1)

    ! var_clw
    call ufo_geovals_get_var_by_index(hofxdiags, hofxdiag_index_2(i), geoval)
    dBT_dql(:) = zero
    dBT_dql(:) = geoval % vals(:,1)

    if (config % RTTOV_mwscattSwitch) then
      ! Get liquid ice jacobian - var_cli
      call ufo_geovals_get_var_by_index(hofxdiags, hofxdiag_index_3(i), geoval)
      dBT_dqi(:) = zero
      dBT_dqi(:) = geoval % vals(:,1)

      H_matrix(i,profindex % qt(1):profindex % qt(2)) = &
              (dBT_dq(:)  * dq_dqt(:) + &
               dBT_dql(:) * dql_dqt(:) + &
               dBT_dqi(:) * dqi_dqt(:)) * q_kgkg(:)
    else
      H_matrix(i,profindex % qt(1):profindex % qt(2)) = &
              (dBT_dq(:)  * dq_dqt(:) + &
               dBT_dql(:) * dql_dqt(:) ) * q_kgkg(:)
    end if

  end do

  ! Clean up
  deallocate(q_kgkg)
  deallocate(pressure)
  deallocate(dq_dqt)
  deallocate(dql_dqt)
  deallocate(dqi_dqt)
  deallocate(dBT_dq)
  deallocate(dBT_dql)
  deallocate(hofxdiag_index_2)
  if (allocated(hofxdiag_index_3)) deallocate(hofxdiag_index_3)

end if

!----
! 2.) Single-valued variables
!----

! 2.1) Surface Temperature - var_sfc_t2m = "surface_temperature"
if (profindex % t2 > 0) then
  write(basename,"(3a)") "brightness_temperature_jacobian_", trim(var_sfc_t2m), "_"
  call ufo_rttovonedvarcheck_geovals_index_by_channels(channels, trim(basename), &
                                          hofxdiags % variables, hofxdiag_index_1)
  do i = 1, nchans
    call ufo_geovals_get_var_by_index(hofxdiags, hofxdiag_index_1(i), geoval)
    H_matrix(i,profindex % t2) = geoval % vals(1,1)
  end do
end if

! 2.2) Water vapour - var_sfc_q2m = "specific_humidity_at_two_meters_above_surface" ! (kg/kg)
if (profindex % q2 > 0) then
  s2m_kgkg = zero
  call ufo_geovals_get_var(geovals, var_sfc_q2m, geoval)
  s2m_kgkg = geoval%vals(1, 1)
  write(basename,"(3a)") "brightness_temperature_jacobian_", trim(var_sfc_q2m), "_"  ! kg/kg
  call ufo_rttovonedvarcheck_geovals_index_by_channels(channels, trim(basename), &
                                          hofxdiags % variables, hofxdiag_index_1)
  do i = 1, nchans
    call ufo_geovals_get_var_by_index(hofxdiags, hofxdiag_index_1(i), geoval)
    H_matrix(i,profindex % q2) = geoval % vals(1,1) * s2m_kgkg
  end do
end if

! 2.3) Surface pressure - var_ps = "surface_pressure" ! (Pa)
if (profindex % pstar > 0) then
  write(basename,"(3a)") "brightness_temperature_jacobian_", trim(var_ps), "_"
  call ufo_rttovonedvarcheck_geovals_index_by_channels(channels, trim(basename), &
                                          hofxdiags % variables, hofxdiag_index_1)
  do i = 1, nchans
    call ufo_geovals_get_var_by_index(hofxdiags, hofxdiag_index_1(i), geoval)
    H_matrix(i,profindex % pstar) = geoval % vals(1,1)
  end do
end if

! 2.4) Windspeed - var_sfc_u10 = "uwind_at_10m"
!                - var_sfc_v10 = "vwind_at_10m"
!                - windsp = sqrt (u*u + v*v)
if (profindex % windspeed > 0) then
  call ufo_geovals_get_var(geovals, trim(var_sfc_u10), geoval)
  u = geoval % vals(1, 1)
  call ufo_geovals_get_var(geovals, trim(var_sfc_v10), geoval)
  v = geoval % vals(1, 1)
  windsp = sqrt (u ** 2 + v ** 2)

  write(basename,"(3a)") "brightness_temperature_jacobian_", trim(var_sfc_u10), "_"
  call ufo_rttovonedvarcheck_geovals_index_by_channels(channels, trim(basename), &
                                          hofxdiags % variables, hofxdiag_index_1)

  allocate(hofxdiag_index_2(nchans))
  write(basename,"(3a)") "brightness_temperature_jacobian_", trim(var_sfc_v10), "_"
  call ufo_rttovonedvarcheck_geovals_index_by_channels(channels, trim(basename), &
                                          hofxdiags % variables, hofxdiag_index_2)

  do i = 1, nchans
    ! var_sfc_u10
    call ufo_geovals_get_var_by_index(hofxdiags, hofxdiag_index_1(i), geoval)
    dBT_du = geoval % vals(1,1)
    ! var_sfc_v10
    call ufo_geovals_get_var_by_index(hofxdiags, hofxdiag_index_2(i), geoval)
    dBT_dv = geoval % vals(1,1)
    if (windsp > zero) then
      ! directional derivation of the Jacobian
      H_matrix(i,profindex % windspeed) = (dBT_du * u + dBT_dv * v) / windsp
    else
      H_matrix(i,profindex % windspeed) = zero
    end if
  end do
  deallocate(hofxdiag_index_2)
end if

! 2.5) Skin temperature - var_sfc_tskin = "skin_temperature"  ! (K)
if (profindex % tstar > 0) then
  write(basename,"(3a)") "brightness_temperature_jacobian_", trim(var_sfc_tskin), "_"
  call ufo_rttovonedvarcheck_geovals_index_by_channels(channels, trim(basename), &
                                          hofxdiags % variables, hofxdiag_index_1)
  do i = 1, nchans
    call ufo_geovals_get_var_by_index(hofxdiags, hofxdiag_index_1(i), geoval)
    H_matrix(i,profindex % tstar) = geoval % vals(1,1)
  end do
end if

! 2.6) Cloud top pressure
if (profindex % cloudtopp > 0) then
  write(basename,"(a)") "brightness_temperature_jacobian_cloud_top_pressure_"
  call ufo_rttovonedvarcheck_geovals_index_by_channels(channels, trim(basename), &
                                          hofxdiags % variables, hofxdiag_index_1)
  do i = 1, nchans
    call ufo_geovals_get_var_by_index(hofxdiags, hofxdiag_index_1(i), geoval)
    H_matrix(i,profindex % cloudtopp) = geoval % vals(1,1)
  end do
end if

! 2.7) Cloud fraction
if (profindex % cloudfrac > 0) then
  write(basename,"(a)") "brightness_temperature_jacobian_cloud_fraction_"
  call ufo_rttovonedvarcheck_geovals_index_by_channels(channels, trim(basename), &
                                          hofxdiags % variables, hofxdiag_index_1)
  do i = 1, nchans
    call ufo_geovals_get_var_by_index(hofxdiags, hofxdiag_index_1(i), geoval)
    H_matrix(i,profindex % cloudfrac) = geoval % vals(1,1)
  end do
end if

!----
! 3.) Emissivities
! This has been left in for future development
!----

! 3.1 Microwave Emissivity - var_sfc_emiss = "surface_emissivity"
if (profindex % mwemiss(1) > 0) then
  ! The emissivity matrix needs to be "unpacked" as it is only one
  ! dimensional over channels - implying you want to retrieve a
  ! single emissivity value. It is unpacked here so that each emissivity
  ! has a corresponding entry for the relevant channel. Note that this is
  ! a bit physically dubious as several channels have the same frequency, etc.
  ! This complexity is dealt with in the B Matrix.
  ! Check that we want only the diagonal elements to be non-zero
  emissloop: do j = 1, size(config % EmissToChannelMap)
    chan = config % EmissToChannelMap(j)
    do i = 1, nchans
      if (channels(i) == chan) then
        write(varname,"(3a,i0)") "brightness_temperature_jacobian_",trim(var_sfc_emiss),"_",channels(i)
        call ufo_geovals_get_var(hofxdiags, varname, geoval)
        H_matrix(i,profindex % mwemiss(1) + j - 1) = geoval % vals(1,1)
        cycle emissloop
      end if
    end do
  end do emissloop
end if

! 3.2. Infrared Emissivity - var_sfc_emiss = "surface_emissivity"
if (profindex % emisspc(1) > 0) THEN
  allocate(emissivity_k(nchans))
  allocate(emissivity(nchans))
  write(basename,"(3a)") "brightness_temperature_jacobian_", trim(var_sfc_emiss), "_"
  call ufo_rttovonedvarcheck_geovals_index_by_channels(channels, trim(basename), &
                                          hofxdiags % variables, hofxdiag_index_1)
  do i = 1, nchans ! loop over channels used
    call ufo_geovals_get_var_by_index(hofxdiags, hofxdiag_index_1(i), geoval)
    emissivity_k(i) = geoval % vals(1,1)
  end do

  ! create emissivity subset from all channels
  call ufo_rttovonedvarcheck_all_to_subset_by_channels(ob % channels_all, &
                             ob % emiss, channels, emissivity)

  nemisspc = profindex % emisspc(2) - profindex % emisspc(1) + 1
  call ob % pcemiss_object % emisktopc (nchans,                                                           & ! in
                                        channels,                                                         & ! in
                                        nemisspc,                                                         & ! in
                                        emissivity,                                                       & ! in
                                        emissivity_K(:),                                                  & ! in
                                        H_matrix(1:nchans,profindex % emisspc(1):profindex % emisspc(2)))   ! out
  deallocate(emissivity)
  deallocate(emissivity_k)
end if

! Here for diagnostics
if (config % FullDiagnostics) then
  call ufo_rttovonedvarcheck_PrintHmatrix( &
    nchans,   &                  ! in
    profindex % nprofelements, & ! in
    channels, &                  ! in
    H_matrix, &                  ! in
    profindex )                  ! in
end if

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
write( unit=int_fmt,fmt='(a)' ) '(' // trim(txt_nchans) // 'I30)'
write( unit=real_fmt,fmt='(a)' ) '(' // trim(txt_nchans) // 'E30.15)'

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

  if ( profindex % emisspc(1) > 0) THEN
    write(*, '(a)') 'PC emissivity retrieval'
    do i = profindex%emisspc(1),profindex%emisspc(2)
      write(*, real_fmt)  H_matrix(:,i)
    end do
  end if

  if ( profindex % cloudtopp > 0 ) THEN
    write(*, '(a)') 'Cloud top pressure'
    write(*, real_fmt)  H_matrix(:,profindex % cloudtopp)
  end if

  if ( profindex % cloudfrac > 0 ) THEN
    write(*, '(a)') 'Cloud fraction'
    write(*, real_fmt)  H_matrix(:,profindex % cloudfrac)
  end if

write(*,*)
write(*, '(a)') 'End H-Matrix'
write(*, '(a)') '------------------------'

end  subroutine ufo_rttovonedvarcheck_PrintHmatrix

! ------------------------------------------------------------

end module ufo_rttovonedvarcheck_minimize_jacobian_mod
