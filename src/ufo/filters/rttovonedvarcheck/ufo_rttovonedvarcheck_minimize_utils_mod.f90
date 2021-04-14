! (C) Copyright 2020 Met Office
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module containing subroutines used by both minimizers.

module ufo_rttovonedvarcheck_minimize_utils_mod

use kinds
use ufo_constants_mod, only: grav, zero, t0c, half, one, two, min_q
use fckit_log_module, only : fckit_log
use ufo_geovals_mod
use ufo_rttovonedvarcheck_constants_mod
use ufo_rttovonedvarcheck_ob_mod
use ufo_rttovonedvarcheck_profindex_mod
use ufo_rttovonedvarcheck_rsubmatrix_mod
use ufo_rttovonedvarcheck_utils_mod, only: ufo_rttovonedvarcheck
use ufo_vars_mod
use ufo_utils_mod, only: Ops_SatRad_Qsplit, Ops_QSat, Ops_QSatWat, cmp_strings

implicit none
private

! subroutines - public
public ufo_rttovonedvarcheck_GeoVaLs2ProfVec
public ufo_rttovonedvarcheck_ProfVec2GeoVaLs
public ufo_rttovonedvarcheck_check_geovals
public ufo_rttovonedvarcheck_CostFunction
public ufo_rttovonedvarcheck_CheckIteration
public ufo_rttovonedvarcheck_CheckCloudyIteration
public ufo_rttovonedvarcheck_PrintIterInfo

character(len=max_string) :: message

contains

!-------------------------------------------------------------------------------
!> Copy geovals data (and ob) to profile.
!!
!! \details Heritage: Ops_SatRad_RTprof2Vec_RTTOV12.f90
!!
!! Convert profile data from the GeoVaLs format (and ob) into the minimisation
!! format vector profile. We only copy fields that are being retrieved, as indicated by
!! the profindex structure.
!!
!! \author Met Office
!!
!! \date 09/06/2020: Created
!!
subroutine ufo_rttovonedvarcheck_GeoVaLs2ProfVec( geovals, & ! in
                                             profindex,    & ! in
                                             ob,           & ! in
                                             prof_x )        ! out

implicit none

! subroutine arguments:
type(ufo_geovals), intent(in)    :: geovals   !< model data at obs location
type(ufo_rttovonedvarcheck_profindex), intent(in) :: profindex !< index array for x vector
type(ufo_rttovonedvarcheck_ob), intent(in) :: ob   !< satellite metadata
real(kind_real), intent(out)     :: prof_x(:) !< x vector

! Local arguments:
character(len=*), parameter :: RoutineName = "ufo_rttovonedvarcheck_GeoVaLs2ProfVec"
character(len=max_string)   :: varname

type(ufo_geoval), pointer    :: geoval
integer                      :: nlevels
integer                      :: ii
real(kind_real), allocatable :: humidity_total(:)
real(kind_real), allocatable :: emiss_pc(:)
real(kind_real)              :: u, v               ! components for windspeed calculation

!-------------------------------------------------------------------------------

prof_x(:) = zero
nlevels = profindex % nlevels

!------------------------------
! 2. Set multi-level variables
!------------------------------

! Note: GeoVaLs are surface -> TOA; prof_x is TOA -> surface

! Note that if number of profile levels is less than number of pressure levels
! we assume the levels are from the surface upwards (remember that RTTOV levels
! are upside-down)

! var_ts - air_temperature - K
if (profindex % t(1) > 0) then
  call ufo_geovals_get_var(geovals, var_ts, geoval)
  prof_x(profindex % t(1):profindex % t(2)) = geoval%vals(nlevels:1:-1, 1) ! K
end if

! var_q - specific_humidity - kg/kg
! for retrieval is ln(g/kg)
if (profindex % q(1) > 0) then
  call ufo_geovals_get_var(geovals, var_q, geoval)
  prof_x(profindex % q(1):profindex % q(2)) = &
         log (geoval%vals(nlevels:1:-1, 1) * 1000.0_kind_real) ! ln(g/kg)
end if

! var_q - specific_humidity - kg/kg
! var_clw  = "mass_content_of_cloud_liquid_water_in_atmosphere_layer" - kg/kg
! var_cli  = "mass_content_of_cloud_ice_in_atmosphere_layer" - kg/kg
! for retrieval is ln(g/kg)
if (profindex % qt(1) > 0) then
  allocate(humidity_total(nlevels))
  humidity_total(:) = zero
  
  ! Get humidity data from geovals
  call ufo_geovals_get_var(geovals, var_q, geoval)
  humidity_total(:) = humidity_total(:) + geoval%vals(nlevels:1:-1, 1)
  call ufo_geovals_get_var(geovals, var_clw, geoval)
  humidity_total(:) = humidity_total(:) + geoval%vals(nlevels:1:-1, 1)
  call ufo_geovals_get_var(geovals, var_cli, geoval)
  humidity_total(:) = humidity_total(:) + geoval%vals(nlevels:1:-1, 1)
  
  ! Convert from kg/kg to ln(g/kg)
  prof_x(profindex % qt(1):profindex % qt(2)) = log (humidity_total * 1000.0_kind_real) ! ln(g/kg)

  deallocate(humidity_total)
end if

!----------------------------
! 3. Single-valued variables
!----------------------------

! var_sfc_t2m = "surface_temperature"
if (profindex % t2 > 0) then
  call ufo_geovals_get_var(geovals, var_sfc_t2m, geoval)
  prof_x(profindex % t2) = geoval%vals(1, 1)
end if

! var_sfc_q2m = "specific_humidity_at_two_meters_above_surface" (kg/kg)
! for retrieval is ln(g/kg)
if (profindex % q2 > 0) then
  call ufo_geovals_get_var(geovals, var_sfc_q2m, geoval)
  prof_x(profindex % q2) = log (geoval%vals(1, 1) * 1000.0_kind_real) ! ln(g/kg)
end if

! var_sfc_p2m = "air_pressure_at_two_meters_above_surface" ! (Pa)
if (profindex % pstar > 0) then
  call ufo_geovals_get_var(geovals, var_sfc_p2m, geoval)
  prof_x(profindex % pstar) = geoval%vals(1, 1) / 100.0_kind_real  ! Pa to hPa
end if

! var_sfc_tskin = "skin_temperature"  ! (K)
if (profindex % tstar > 0) then
  call ufo_geovals_get_var(geovals, var_sfc_tskin, geoval)
  prof_x(profindex % tstar) = geoval%vals(1, 1)
end if

! This has been left in for future development
! cloud top pressure
!if (profindex % cloudtopp > 0) then
!  prof_x(profindex % cloudtopp) = ob % cloudtopp ! carried around as hPa
!end if

! This has been left in for future development
! cloud fraction
!if (profindex % cloudfrac > 0) then
!  prof_x(profindex % cloudfrac) = ob % cloudfrac
!end if

! Windspeed - var_u = "eastward_wind"
!           - var_v = "northward_wind"
!           - windsp = sqrt (u*u + v*v)
if (profindex % windspeed > 0) then
  call ufo_geovals_get_var(geovals, trim(var_u), geoval)
  u = geoval % vals(1, 1)
  call ufo_geovals_get_var(geovals, trim(var_v), geoval)
  v = geoval % vals(1, 1)
  prof_x(profindex % windspeed) = sqrt(u ** 2 + v ** 2)
end if

!----------------------------
! 4. Emissivities
! This has been left in for future development
!----------------------------

! Microwave Emissivity
!if (profindex % mwemiss(1) > 0) then
!  ! Check that emissivity map is the correct size for the profile
!  if ((profindex % mwemiss(2) - profindex % mwemiss(1) + 1) /= size(EmissMap)) then
!    call abor1_ftn("mwemiss size differs from emissivity map")
!  end if
!  ! Copy microwave emissivity to profile
!    prof_x(profindex % mwemiss(1):profindex % mwemiss(2)) = ob % emiss(EmissMap)
!end if

! Retrieval of emissivity principal components
!if (profindex % emisspc(1) > 0) THEN
!  ! convert ob % emiss to emiss pc
!  allocate(emiss_pc(profindex % emisspc(2)-profindex % emisspc(1)+1))
!  call ob % pcemis % emistoPC(ob % channels_used(:), ob % emiss(:), emiss_pc(:))
!  prof_x(profindex % emisspc(1):profindex % emisspc(2)) = emiss_pc
!  deallocate(emiss_pc)
!end if

end subroutine ufo_rttovonedvarcheck_GeoVaLs2ProfVec

!-------------------------------------------------------------------------------
!> Copy profile data to geovals (and ob).
!!
!! \details Heritage: Ops_SatRad_Vec2RTprof_RTTOV12.f90
!!
!! Convert profile data to the GeoVaLs (and ob) format.  We only copy fields 
!! that are being retrieved, as indicated by the profindex structure.
!!
!! \author Met Office
!!
!! \date 09/06/2020: Created
!!
subroutine ufo_rttovonedvarcheck_ProfVec2GeoVaLs(geovals, & ! inout
                                            profindex,    & ! in
                                            ob,           & ! inout
                                            prof_x,       & ! in
                                            UseQtsplitRain) ! in

implicit none

! subroutine arguments:
type(ufo_geovals), intent(inout) :: geovals   !< model data at obs location
type(ufo_rttovonedvarcheck_profindex), intent(in) :: profindex !< index array for x vector
type(ufo_rttovonedvarcheck_ob), intent(inout) :: ob   !< satellite metadata
real(kind_real), intent(in)      :: prof_x(:) !< x vector
logical, intent(in) :: UseQtsplitRain !< flag to choose whether to split rain in qt

! Local arguments:
character(len=*), parameter  :: RoutineName = "ufo_rttovonedvarcheck_ProfVec2GeoVaLs"
character(len=max_string)    :: varname
integer                      :: gv_index, i, ii
integer                      :: nlevels
integer                      :: EmissElement
type(ufo_geoval), pointer    :: geoval
real(kind_real), allocatable :: pressure(:)
real(kind_real), allocatable :: humidity_total(:)
real(kind_real), allocatable :: q(:)
real(kind_real), allocatable :: ql(:)
real(kind_real), allocatable :: qi(:)
real(kind_real), allocatable :: emiss_pc(:)
real(kind_real)              :: u, v, windsp  ! variable needed for the windspeed calculation

!-------------------------------------------------------------------------------

nlevels = profindex % nlevels

!------------------------------
! 2. Set multi-level variables
!------------------------------

! Note GeoVaLs are surface -> TOA ; profile is TOA -> surface

! Note that if number of profile levels is less than number of pressure levels
! we assume the levels are from the surface upwards (remember that RTTOV levels
! are upside-down)

! var_ts - air_temperature - K
if (profindex % t(1) > 0) then
  gv_index = 0
  do i=1,geovals%nvar
    if (cmp_strings(var_ts, geovals%variables(i))) gv_index = i
  end do
  geovals%geovals(gv_index)%vals(nlevels:1:-1,1) = prof_x(profindex % t(1):profindex % t(2)) ! K
end if

! var_q = "specific_humidity" ! kg/kg
! for retrieval is ln(g/kg)
if (profindex % q(1) > 0) then
  gv_index = 0
  do i=1,geovals%nvar
    if (cmp_strings(var_q, geovals%variables(i))) gv_index = i
  end do
  geovals%geovals(gv_index)%vals(nlevels:1:-1,1) = EXP (prof_x(profindex % q(1):profindex % q(2))) / &
                                                  1000.0_kind_real ! ln(g/kg) => kg/kg
end if

! var_q = "specific_humidity" ! kg/kg
! var_clw  = "mass_content_of_cloud_liquid_water_in_atmosphere_layer" - kg/kg
! var_cli  = "mass_content_of_cloud_ice_in_atmosphere_layer" - kg/kg
! for retrieval is ln(g/kg)
if (profindex % qt(1) > 0) then
  nlevels = profindex % nlevels
  allocate(pressure(nlevels))
  allocate(humidity_total(nlevels))
  allocate(q(nlevels))
  allocate(ql(nlevels))
  allocate(qi(nlevels))
  
  ! Convert from ln(g/kg) to kg/kg
  humidity_total(nlevels:1:-1) = EXP (prof_x(profindex % qt(1):profindex % qt(2))) / &
                                1000.0_kind_real ! ln(g/kg) => kg/kg

  ! var_prs  = "air_pressure" Pa
  call ufo_geovals_get_var(geovals, var_prs, geoval)
  pressure(:) = geoval%vals(:, 1)    ! Pa

  ! Split qtotal to q(water_vapour), q(liquid), q(ice)
  call Ops_SatRad_Qsplit ( 1,             &
                    pressure(:),          &
                    ob % background_T(:), &
                    humidity_total,       &
                    q(:),                 &
                    ql(:),                &
                    qi(:),                &
                    UseQtsplitRain)

  ! Assign values to geovals
  gv_index = 0
  do i=1,geovals%nvar
    if (cmp_strings(var_q, geovals%variables(i))) gv_index = i
  end do
  geovals%geovals(gv_index)%vals(:,1) = q(:)

  gv_index = 0
  do i=1,geovals%nvar
    if (cmp_strings(var_clw, geovals%variables(i))) gv_index = i
  end do
  geovals%geovals(gv_index)%vals(:,1) = ql(:)

  gv_index = 0
  do i=1,geovals%nvar
    if (cmp_strings(var_cli, geovals%variables(i))) gv_index = i
  end do
  geovals%geovals(gv_index)%vals(:,1) = qi(:)

  deallocate(pressure)
  deallocate(humidity_total)
  deallocate(q)
  deallocate(ql)
  deallocate(qi)

end if

!----------------------------
! 3. Single-valued variables
!----------------------------

! var_sfc_t2m = "surface_temperature"
if (profindex % t2 > 0) then
  gv_index = 0
  do i=1,geovals%nvar
    if (cmp_strings(var_sfc_t2m, geovals%variables(i))) gv_index = i
  end do
  geovals%geovals(gv_index)%vals(1,1) = prof_x(profindex % t2) ! K
end if

! var_sfc_q2m = "specific_humidity_at_two_meters_above_surface" ! (kg/kg)
! for retrieval is ln(g/kg)
if (profindex % q2 > 0) then
  gv_index = 0
  do i=1,geovals%nvar
    if (cmp_strings(var_sfc_q2m, geovals%variables(i))) gv_index = i
  end do
  geovals%geovals(gv_index)%vals(1,1) = EXP (prof_x(profindex % q2)) / 1000.0_kind_real ! ln(g/kg) => kg/kg
end if

! var_sfc_p2m = "air_pressure_at_two_meters_above_surface" ! (Pa)
if (profindex % pstar > 0) then
  gv_index = 0
  do i=1,geovals%nvar
    if (cmp_strings(var_sfc_p2m, geovals%variables(i))) gv_index = i
  end do
  geovals%geovals(gv_index)%vals(1,1) = prof_x(profindex % pstar) * 100.0_kind_real
end if

! var_sfc_tskin = "skin_temperature"  ! (K)
if (profindex % tstar > 0) then
  gv_index = 0
  do i=1,geovals%nvar
    if (cmp_strings(var_sfc_tskin, geovals%variables(i))) gv_index = i
  end do
  geovals%geovals(gv_index)%vals(1,1) = prof_x(profindex % tstar)
end if

! This has been left in for future development
! cloud top pressure - passed through via the ob
!if (profindex % cloudtopp > 0) then
!  ob % cloudtopp = prof_x(profindex % cloudtopp) ! stored in ob as hPa
!end if

! This has been left in for future development
! cloud fraction - passed through via the ob
!if (profindex % cloudfrac > 0) then
!  ob % cloudfrac = prof_x(profindex % cloudfrac)
!end if

! windspeed
if (profindex % windspeed > 0) then
  call ufo_geovals_get_var(geovals, trim(var_u), geoval)
  u = geoval % vals(1, 1)
  call ufo_geovals_get_var(geovals, trim(var_v), geoval)
  v = geoval % vals(1, 1)
  windsp = sqrt (u ** 2 + v ** 2)

  ! The ratio of new windsp to old windsp gives the fractional change.
  ! This is then applied to each component of the wind.
  if (windsp > zero) then
    u = u * (prof_x(profindex % windspeed) / windsp)
    v = v * (prof_x(profindex % windspeed) / windsp)
  else
    u = zero
    v = zero
  end if

  ! Write back updated u component
  gv_index = 0
  do i=1,geovals%nvar
    if (trim(var_u) == trim(geovals%variables(i))) gv_index = i
  end do
  geovals%geovals(gv_index)%vals(1,1) = u

  ! Write back updated v component
  gv_index = 0
  do i=1,geovals%nvar
    if (trim(var_v) == trim(geovals%variables(i))) gv_index = i
  end do
  geovals%geovals(gv_index)%vals(1,1) = v
end if

!----------------------------
! 4. Emissivities
! This has been left in for future development
!------------------------

! Retrieval of microwave emissivity directly
!if (profindex % mwemiss(1) > 0) THEN
!  do ii = 1, size(ob % channels_used)
!    EmissElement = EmissElements(ob % channels_used(ii))
!    ob % emiss(ii) = prof_x(profindex % mwemiss(1) + EmissElement - 1)
!  end do
!end if

! Retrieval of emissivity principal components
!if (profindex % emisspc(1) > 0) THEN
!  allocate(emiss_pc(profindex % emisspc(2)-profindex % emisspc(1)+1))
!  emiss_pc = prof_x(profindex % emisspc(1):profindex % emisspc(2))
!  ! convert emiss_pc to ob % emissivity using
!  call ob % pcemis % pctoemis(size(ob % channels_used), ob % channels_used, &
!                              size(emiss_pc), emiss_pc(:), ob % emiss(:))
!  deallocate(emiss_pc)
!end if

end subroutine ufo_rttovonedvarcheck_ProfVec2GeoVaLs

!-------------------------------------------------------------------------------
!> Check the geovals are ready for the first iteration
!!
!! \details Heritage: Ops_SatRad_SetUpRTprofBg_RTTOV12.f90
!!
!! Check the geovals profile is ready for the first iteration.  The
!! only check included at the moment if the first calculation for 
!! q total.
!!
!! \author Met Office
!!
!! \date 09/06/2020: Created
!!
subroutine ufo_rttovonedvarcheck_check_geovals(self, geovals, profindex, surface_type)

implicit none

! subroutine arguments:
type(ufo_rttovonedvarcheck), intent(in) :: self
type(ufo_geovals), intent(inout) :: geovals   !< model data at obs location
type(ufo_rttovonedvarcheck_profindex), intent(in) :: profindex !< index array for x vector
integer, intent(in) :: surface_type  !< surface type for cold surface check

character(len=*), parameter  :: routinename = "ufo_rttovonedvarcheck_check_geovals"
character(len=max_string)    :: varname
type(ufo_geoval), pointer    :: geoval
integer                      :: gv_index, i     ! counters
integer                      :: nlevels
real(kind_real), allocatable :: temperature(:)  ! temperature (K)
real(kind_real), allocatable :: pressure(:)     ! pressure (Pa)
real(kind_real), allocatable :: qsaturated(:)
real(kind_real), allocatable :: humidity_total(:)
real(kind_real), allocatable :: q(:)            ! specific humidity (kg/kg)
real(kind_real), allocatable :: ql(:)
real(kind_real), allocatable :: qi(:)
real(kind_real)              :: skin_t, pressure_2m, temperature_2m, NewT
integer                      :: level_1000hpa, level_950hpa

write(message, *) routinename, " : started"
call fckit_log % debug(message)

! -------------------------------------------
! Load variables needed by multiple routines
! -------------------------------------------

nlevels = profindex % nlevels
allocate(temperature(nlevels))
allocate(pressure(nlevels))
call ufo_geovals_get_var(geovals, trim(var_ts), geoval)
temperature(:) = geoval%vals(:, 1) ! K
call ufo_geovals_get_var(geovals, trim(var_prs), geoval)
pressure(:) = geoval%vals(:, 1)    ! Pa

!---------------------------------------------------
! 1. Make sure q and q2m does not exceed saturation
!---------------------------------------------------
if (profindex % q(1) > 0 .or. profindex % qt(1) > 0) then
  allocate(q(nlevels))
  allocate(qsaturated(nlevels))

  ! Get humidity - kg/kg
  varname = trim(var_q)
  call ufo_geovals_get_var(geovals, varname, geoval)
  q(:) = geoval%vals(:, 1)

  ! Calculated saturation humidity
  if (self % UseRHwaterForQC) then
    call Ops_QsatWat (qsaturated(:),   & ! out
                      temperature(:),  & ! in
                      pressure(:),     & ! in
                      nlevels)           ! in
  else
    call Ops_Qsat (qsaturated(:),   & ! out
                   temperature(:),  & ! in
                   pressure(:),     & ! in
                   nlevels)           ! in  
  end if

  ! Limit q
  where (q > qsaturated)
    q = qsaturated
  end where
  where (q < min_q)
    q = min_q
  end where

  ! Assign values to geovals q
  gv_index = 0
  do i=1,geovals%nvar
    if (cmp_strings(varname, geovals%variables(i))) gv_index = i
  end do
  geovals%geovals(gv_index) % vals(:,1) = q(:)

  deallocate(q)
  deallocate(qsaturated)
end if

if (profindex % q2 > 0) then
  allocate(q(1))
  allocate(qsaturated(1))

  ! Get humidity - kg/kg
  varname = trim(var_sfc_q2m)
  call ufo_geovals_get_var(geovals, varname, geoval)
  q(1) = geoval%vals(1, 1)

  ! Calculated saturation humidity
  if (self % UseRHwaterForQC) then
    call Ops_QsatWat (qsaturated(1:1),   & ! out
                      temperature(1:1),  & ! in
                      pressure(1:1),     & ! in
                      1)                   ! in
  else
    call Ops_Qsat (qsaturated(1:1),   & ! out
                   temperature(1:1),  & ! in
                   pressure(1:1),     & ! in
                   1)                   ! in
  end if

  ! Limit q
  if (q(1) > qsaturated(1)) q(1) = qsaturated(1)
  if (q(1) < min_q) q(1) = min_q

  ! Assign values to geovals q
  gv_index = 0
  do i=1,geovals%nvar
    if (cmp_strings(varname, geovals%variables(i))) gv_index = i
  end do
  geovals%geovals(gv_index) % vals(1,1) = q(1)

  deallocate(q)
  deallocate(qsaturated)
end if

!-------------------------
! 2. Specific humidity total
!-------------------------

if (profindex % qt(1) > 0) then

  allocate(humidity_total(nlevels))
  allocate(q(nlevels))
  allocate(ql(nlevels))
  allocate(qi(nlevels))

  humidity_total(:) = 0.0
  ! Water vapour
  call ufo_geovals_get_var(geovals, trim(var_q), geoval)
  humidity_total(:) = humidity_total(:) + geoval%vals(:, 1)

  ! Cloud liquid water
  call ufo_geovals_get_var(geovals, trim(var_clw), geoval)
  where (geoval%vals(:, 1) < 0.0) geoval%vals(:, 1) = 0.0
  humidity_total(:) = humidity_total(:) + geoval%vals(:, 1)

  ! Split qtotal to q(water_vapour), q(liquid), q(ice)
  call Ops_SatRad_Qsplit ( 1,       &
                    pressure(:),    &
                    temperature(:), &
                    humidity_total, &
                    q(:),           &
                    ql(:),          &
                    qi(:),          &
                    self % UseQtsplitRain)

  ! Assign values to geovals q
  varname = trim(var_q)  ! kg/kg
  gv_index = 0
  do i=1,geovals%nvar
    if (cmp_strings(varname, geovals%variables(i))) gv_index = i
  end do
  geovals%geovals(gv_index) % vals(:,1) = q(:)

  ! Assign values to geovals q clw
  varname = trim(var_clw)  ! kg/kg
  gv_index = 0
  do i=1,geovals%nvar
    if (cmp_strings(varname, geovals%variables(i))) gv_index = i
  end do
  geovals%geovals(gv_index) % vals(:,1) = ql(:)

  ! Assign values to geovals ciw
  varname = trim(var_cli)  ! kg/kg
  gv_index = 0
  do i=1,geovals%nvar
    if (cmp_strings(varname, geovals%variables(i))) gv_index = i
  end do
  geovals%geovals(gv_index)%vals(:,1) = qi(:)

  deallocate(humidity_total)
  deallocate(q)
  deallocate(ql)
  deallocate(qi)

end if

!----
! Legacy Ops_SatRad_SetUpRTprofBg_RTTOV12.f90 - done here to make sure geovals are updated
! Reset low level temperatures over seaice and cold, low land as per Ops_SatRad_SetUpRTprofBg.F90
! N.B. I think this should be flagged so it's clear that the background has been modified
!----
if(surface_type /= RTsea .and. self % UseColdSurfaceCheck) then

  ! Get skin temperature
  call ufo_geovals_get_var(geovals, var_sfc_tskin, geoval)
  skin_t = geoval%vals(1, 1)

  ! Get 2m pressure
  call ufo_geovals_get_var(geovals, var_sfc_p2m, geoval)
  pressure_2m = geoval%vals(1, 1)

  ! Get 2m temperature
  call ufo_geovals_get_var(geovals, var_sfc_t2m, geoval)
  temperature_2m = geoval%vals(1, 1)

  if(skin_t < 271.4_kind_real .and. &
     pressure_2m  > 95000.0_kind_real) then

     level_1000hpa = minloc(abs(pressure - 100000.0_kind_real),DIM=1)
     level_950hpa = minloc(abs(pressure - 95000.0_kind_real),DIM=1)

     NewT = temperature(level_950hpa)
     if(pressure_2m > 100000.0_kind_real) then
       NewT = max(NewT, temperature(level_1000hPa))
     end if
     NewT = min(NewT, 271.4_kind_real)

     temperature(level_1000hPa) = max(temperature(level_1000hPa), NewT)
     temperature_2m = max(temperature_2m, NewT)
     skin_t = max(skin_t, NewT)

     ! Put updated values back into geovals
     ! Temperature
     gv_index = 0
     varname = trim(var_ts)
     do i=1,geovals%nvar
       if (cmp_strings(varname, geovals%variables(i))) gv_index = i
     end do
     geovals%geovals(gv_index) % vals(level_1000hPa,1) = temperature(level_1000hPa)

     ! 2m Temperature
     gv_index = 0
     varname = trim(var_sfc_t2m)
     do i=1,geovals%nvar
       if (cmp_strings(varname, geovals%variables(i))) gv_index = i
     end do
     geovals%geovals(gv_index) % vals(1,1) = temperature_2m

     ! Skin Temperature
     gv_index = 0
     varname = trim(var_sfc_tskin)
     do i=1,geovals%nvar
       if (cmp_strings(varname, geovals%variables(i))) gv_index = i
     end do
     geovals%geovals(gv_index) % vals(1,1) = skin_t

   endif

endif

! Tidy up
if (allocated(temperature))    deallocate(temperature)
if (allocated(pressure))       deallocate(pressure)
if (allocated(qsaturated))     deallocate(qsaturated)
if (allocated(humidity_total)) deallocate(humidity_total)
if (allocated(q))              deallocate(q)
if (allocated(ql))             deallocate(ql)
if (allocated(qi))             deallocate(qi)

write(message, *) routinename, " : ended"
call fckit_log % debug(message)

end subroutine ufo_rttovonedvarcheck_check_geovals

!-------------------------------------------------------------------------------
!> Calculate the cost function.
!!
!! \details Heritage: Ops_SatRad_CostFunction.f90
!!
!! Calculate the cost function from the input delta's and error
!! covariances.
!!
!! \author Met Office
!!
!! \date 09/06/2020: Created
!!
subroutine ufo_rttovonedvarcheck_CostFunction(DeltaProf, b_inv, &
                                              DeltaObs, r_matrix, &
                                              Jcost)

implicit none

! subroutine arguments:
real(kind_real), intent(in)       :: DeltaProf(:)
real(kind_real), intent(in)       :: b_inv(:,:)
real(kind_real), intent(in)       :: DeltaObs(:)
type(ufo_rttovonedvarcheck_rsubmatrix), intent(in) :: r_matrix
real(kind_real), intent(out)      :: Jcost(3)

! Local arguments:
character(len=*), parameter  :: RoutineName = "ufo_rttovonedvarcheck_CostFunction"
integer                      :: y_size
real(kind_real)              :: Jb, Jo, Jcurrent
real(kind_real), allocatable :: RinvDeltaY(:)

!-------------------------------------------------------------------------------

allocate(RinvDeltaY, source=DeltaObs)

y_size = size(DeltaObs)
Jcost(:) = zero
call r_matrix % multiply_inverse_vector(DeltaObs, RinvDeltaY)

Jo = half * dot_product(DeltaObs, RinvDeltaY)
Jb = half * dot_product(DeltaProf, (matmul(b_inv, DeltaProf)))
Jcurrent = Jb + Jo

Jcost(1) = (Jo + Jb) * two / real (y_size, kind_real)   ! Normalize cost by nchans
Jcost(2) = Jb * two / real (y_size, kind_real)          ! Normalize cost by nchans
Jcost(3) = Jo * two / real (y_size, kind_real)          ! Normalize cost by nchans

write(message,*) "Jo, Jb, Jcurrent = ", Jo, Jb, Jcost(1)
call fckit_log % debug(message)

deallocate(RinvDeltaY)

end subroutine ufo_rttovonedvarcheck_CostFunction

!-------------------------------------------------------------------------------
!> Constrain profile humidity and check temperature values are okay
!!
!! \details Heritage: Ops_SatRad_CheckIteration.f90
!!
!! Check humidity and temperature profiles.
!!
!! \author Met Office
!!
!! \date 09/06/2020: Created
!!
subroutine ufo_rttovonedvarcheck_CheckIteration (self, &
                                      geovals,    &
                                      profindex,  &
                                      profile,    &
                                      OutOfRange)

implicit none

! subroutine arguments:
type(ufo_rttovonedvarcheck), intent(in) :: self
type(ufo_geovals), intent(in)           :: geovals
type(ufo_rttovonedvarcheck_profindex), intent(in) :: profindex
real(kind_real), intent(inout)   :: profile(:)
logical, intent(out)             :: OutOfRange

! Local declarations:
real(kind_real), allocatable :: qsaturated(:)
real(kind_real), allocatable :: scaled_qsaturated(:)
real(kind_real)              :: q2_sat(1)
real(kind_real), allocatable :: Plevels_1DVar(:)
real(kind_real)              :: Pstar_Pa(1)
real(kind_real), allocatable :: Temp(:)
real(kind_real)              :: Temp2(1)
real(kind_real)              :: rtbase
integer                      :: nlevels_q
integer                      :: toplevel_q
character(len=*), parameter  :: RoutineName = "ufo_rttovonedvarcheck_CheckIteration"
type(ufo_geoval), pointer    :: geoval
character(len=max_string)    :: varname
integer                      :: ii
integer                      :: nlevels_1dvar

! Setup
OutOfRange = .false.
nlevels_1dvar = self % nlevels
allocate(Temp(nlevels_1dvar))

!---------------------
! 1. Check Temperatures
!---------------------

!----
! 1.1) levels
!----

if (profindex % t(1) > 0) then
  Temp = profile(profindex % t(1):profindex % t(2))
  if (any (Temp < minTemperature) .or. any (Temp > MaxTemperature)) then
    OutOfRange = .true.
  end if
else
  varname = var_ts
  call ufo_geovals_get_var(geovals, varname, geoval)
  ! Note: profile is TOA -> surface
  Temp = geoval%vals(nlevels_1dvar:1:-1, 1) ! K
end if

!----
! 1.2) surface air
!----

if (profindex % t2 > 0) then
  Temp2 = profile(profindex % t2)
  if (Temp2(1) < minTemperature .or. Temp2(1) > MaxTemperature) then
    OutOfRange = .true.
  end if
else
  varname = var_sfc_t2m
  call ufo_geovals_get_var(geovals, varname, geoval)
  Temp2 = geoval%vals(1, 1) ! K
end if

!----
! 1.3) surface skin
!----

if (profindex % tstar > 0) then
  if (profile(profindex % tstar) < minTemperature .or. &
      profile(profindex % tstar) > MaxTemperature) then
  OutOfRange = .true.
  end if
end if

!---------------------
! 2. Constrain humidity
!---------------------

Constrain: if (.not. OutOfRange) then

  ! Work out the number of humidity levels and allocate size to arrays

  if (profindex % q(1) > 0) then
    nlevels_q = profindex % q(2) - profindex % q(1) + 1
  else if (profindex % qt(1) > 0) then
    nlevels_q = profindex % qt(2) - profindex % qt(1) + 1
    allocate (scaled_qsaturated(nlevels_q))
  else
    nlevels_q = nlevels_1dvar
  end if
  toplevel_q = nlevels_1dvar - nlevels_q + 1
  allocate (qsaturated(nlevels_q))


  ! Get pressure
  varname = var_prs
  call ufo_geovals_get_var(geovals, varname, geoval)
  if (.not. allocated(Plevels_1DVar)) allocate(Plevels_1DVar(nlevels_q))
  ! Note: profile is TOA -> surface
  Plevels_1DVar(:) = geoval%vals(nlevels_1dvar:1:-1, 1) ! K

  !----
  ! 2.1) Levels
  !----

  ! Note that if number of profile levels is less than number of pressure levels
  ! we assume the levels are from the surface upwards (remember that RTTOV levels
  ! are upside-down)

  if (profindex % q(1) > 0) then

    if (self % useRHwaterForQC) then
      call Ops_QsatWat(qsaturated(1:nlevels_q),    & ! out (qsat levels)
                       Temp(1:nlevels_q),          & ! in  (t levels)
                       Plevels_1DVar(1:nlevels_q), & ! in  (p levels)
                       nlevels_q)                    ! in
    else
      call Ops_Qsat(qsaturated(1:nlevels_q),    & ! out
                    Temp(1:nlevels_q),          & ! in
                    Plevels_1DVar(1:nlevels_q), & ! in
                    nlevels_q)                    ! in
    end if

    qsaturated(1:nlevels_q) = log (qsaturated(1:nlevels_q) * 1000.0_kind_real)
    where (profile(profindex % q(1):profindex % q(2)) > qsaturated(1:nlevels_q))
      profile(profindex % q(1):profindex % q(2)) = qsaturated(1:nlevels_q)
    end where

  end if

  if (profindex % qt(1) > 0) then

    ! Qtotal is not included in the relative humidity quality control changes
    call Ops_Qsat (qsaturated(1:nlevels_q),    & ! out
                   Temp(1:nlevels_q),          & ! in
                   Plevels_1DVar(1:nlevels_q), & ! in
                   nlevels_q)                    ! in

    ! scaled_qsaturated is generated as an upper limit on the value of qtot , and
    ! is set to 2*qsat.
    scaled_qsaturated(1:nlevels_q) = log (two * qsaturated(1:nlevels_q) * 1000.0_kind_real)
    where (profile(profindex % qt(1):profindex % qt(2)) > scaled_qsaturated(1:nlevels_q))
      profile(profindex % qt(1):profindex % qt(2)) = scaled_qsaturated(1:nlevels_q)
    end where

  end if

  !----
  ! 2.2) Surface
  !----

  if (profindex % q2 > 0) then
    varname = var_sfc_p2m
    call ufo_geovals_get_var(geovals, varname, geoval)
    Pstar_Pa(1) = geoval%vals(1, 1)
    if (self % useRHwaterForQC) then
      call Ops_QsatWat (q2_sat(1:1),   & ! out
                        Temp2,         & ! in
                        Pstar_Pa(1:1), & ! in
                        1)               ! in
    else
      call Ops_Qsat (q2_sat(1:1),   & ! out
                     Temp2,         & ! in
                     Pstar_Pa(1:1), & ! in
                     1)               ! in
    end if
    q2_sat(1) = log (q2_sat(1) * 1000.0_kind_real)
    if (profile(profindex % q2) > q2_sat(1)) then
      profile(profindex % q2) = q2_sat(1)
    end if
  end if

! This has been left in for future development
!  !----
!  ! 2.3) Grey cloud
!  !----
!
!  if (profindex % cloudtopp > 0) then
!    profile(profindex % CloudFrac) = min (profile(profindex % CloudFrac), 1.0)
!    profile(profindex % CloudFrac) = max (profile(profindex % CloudFrac), 0.0)
!    profile(profindex % CloudTopP) = max (profile(profindex % CloudTopP), 100.0)
!    if (LimitCTPtorTBase) then
!      rtbase = maxval (Plevels_RTModel(:))
!      profile(profindex % CloudTopP) = min (profile(profindex % CloudTopP), Pstar_mb, rtbase)
!    else
!      profile(profindex % CloudTopP) = min (profile(profindex % CloudTopP), Pstar_mb)
!    end if
!  end if
!
! This has been left in for future development
!  !----
!  ! 2.4) Surface emissivity PCs
!  !----
!
!  if (profindex % emisspc(1) > 0) then
!    where (profile(profindex % emisspc(1):profindex % emisspc(2)) > EmisEigenVec % PCmax(1:nemisspc))
!      profile(profindex % emisspc(1):profindex % emisspc(2)) = EmisEigenVec % PCmax(1:nemisspc)
!    end where
!    where (profile(profindex % emisspc(1):profindex % emisspc(2)) < EmisEigenVec % PCmin(1:nemisspc))
!      profile(profindex % emisspc(1):profindex % emisspc(2)) = EmisEigenVec % PCmin(1:nemisspc)
!    end where
!  end if
!
! This has been left in for future development
!  !--------
!  ! 2.5) Cloud profiles
!  !--------
!
!  if (profindex % cf(1) > 0) then
!    where (profile(profindex % cf(1):profindex % cf(2)) <= 0.001)
!      profile(profindex % cf(1):profindex % cf(2)) = 0.001
!    end where
!    where (profile(profindex % cf(1):profindex % cf(2)) >= 1.0)
!      profile(profindex % cf(1):profindex % cf(2)) = 0.999
!    end where
!  end if
!
!  if (profindex % ql(1) > 0 .AND. .not. useLogForCloud) then
!    where (profile(profindex % ql(1):profindex % ql(2)) < 0)
!      profile(profindex % ql(1):profindex % ql(2)) = 0
!    end where
!  end if
!
!  if (profindex % qi(1) > 0 .AND. .not. useLogForCloud) then
!    where (profile(profindex % qi(1):profindex % qi(2)) < 0)
!      profile(profindex % qi(1):profindex % qi(2)) = 0
!    end where
!  end if

end if Constrain

! --------
! Tidy up
! --------
if (allocated(qsaturated))        deallocate(qsaturated)
if (allocated(scaled_qsaturated)) deallocate(scaled_qsaturated)
if (allocated(Plevels_1DVar))     deallocate(Plevels_1DVar)
if (allocated(Temp))              deallocate(Temp)

end subroutine ufo_rttovonedvarcheck_CheckIteration

!-------------------------------------------------------------------------------
!> Check cloud during iteration.
!!
!! \details Heritage: Ops_SatRad_CheckCloudyIteration.f90
!!
!! For a 1dvar profile using cloud variables, check that the liquid water path
!! and the ice water path are sensible. if they are beyond sensible values stop
!! 1dvar and reject profile
!!
!! \author Met Office
!!
!! \date 09/06/2020: Created
!!
subroutine ufo_rttovonedvarcheck_CheckCloudyIteration( &
  geovals,       & ! in
  profindex,     & ! in
  nlevels_1dvar, & ! in
  OutOfRange,    & ! out
  OutLWP         ) ! out

implicit none

type(ufo_geovals), intent(in)          :: geovals
type(ufo_rttovonedvarcheck_profindex), intent(in) :: profindex
integer, intent(in)                    :: nlevels_1dvar
logical, intent(out)                   :: OutOfRange
real(kind_real), optional, intent(out) :: OutLWP

! Local variables:
real(kind_real) :: LWP
real(kind_real) :: IWP
real(kind_real) :: dp
real(kind_real) :: meanql
real(kind_real) :: meanqi
integer         :: i
integer         :: nlevels_q, toplevel_q
character(len=*), parameter :: RoutineName = "ufo_rttovonedvarcheck_CheckCloudyIteration"

real(kind_real), parameter   :: MaxLWP = two
real(kind_real), parameter   :: MaxIWP = 3.0_kind_real
real(kind_real)              :: Plevels_1DVar(nlevels_1dvar)
type(ufo_geoval), pointer    :: geoval
real(kind_real)              :: clw(nlevels_1dvar)
real(kind_real)              :: ciw(nlevels_1dvar)
character(len=max_string)    :: varname

!-------------------------------------------------------------------------------

!initialise
OutOfRange = .false.
IWP = zero
LWP = zero

! Get pressure from geovals
varname = var_prs
call ufo_geovals_get_var(geovals, varname, geoval)
Plevels_1DVar(:) = geoval%vals(:, 1) ! K

! Get clw from geovals
varname = var_clw
call ufo_geovals_get_var(geovals, varname, geoval)
clw = geoval%vals(:, 1)

! Get ciw from geovals
varname = var_cli
call ufo_geovals_get_var(geovals, varname, geoval)
ciw = geoval%vals(:, 1)

! Work out the number of humidity levels
if ( profindex % q(1) > 0 ) then
  nlevels_q = profindex % q(2) - profindex % q(1) + 1
else if ( profindex % qt(1) > 0 ) then
  nlevels_q = profindex % qt(2) - profindex % qt(1) + 1
else if ( profindex % ql(1) > 0 ) then
  nlevels_q = profindex % ql(2) - profindex % ql(1) + 1
else if ( profindex % qi(1) > 0 ) then
  nlevels_q = profindex % qi(2) - profindex % qi(1) + 1
else
  nlevels_q = nlevels_1dvar
end if
toplevel_q = nlevels_1dvar - nlevels_q + 1

!is the profile cloudy?
!if it is do the test

if (any(ciw(:) > zero) .or. &
    any(clw(:) > zero)) then

!1.1 compute iwp, lwp

  do i=1, nlevels_1dvar-1

    dp =  Plevels_1DVar(i) - Plevels_1DVar(i+1)

    ! Calculate layer mean from CloudIce on levels
    meanqi = half * &
      (ciw(i) + ciw(i+1))
    if (meanqi > zero) then
      IWP = IWP + dp * meanqi
    end if

    ! Calculate layer mean from CLW on levels
    meanql = half * (clw(i) + clw(i+1))
    if (meanql > zero) then
      LWP = LWP + dp * meanql
    end if

  end do

  IWP = IWP / grav
  LWP = LWP / grav


!2.1 test if lwp iwp exceeds thresholds

  if ((IWP > MaxIWP) .or. (LWP > MaxLWP)) then
    call fckit_log % debug("lwp or iwp exceeds thresholds")
    OutOfRange = .true.
    write(message,*) "lwp and iwp = ",LWP,IWP
    call fckit_log % debug(message)
  else
    call fckit_log % debug("lwp and iwp less than thresholds")
    write(message,*) "lwp and iwp = ",LWP,IWP
    call fckit_log % debug(message)
  end if

end if

if (present(OutLWP)) OutLWP = LWP

end subroutine ufo_rttovonedvarcheck_CheckCloudyIteration

!-----------------------------------------------------------
!> Print detailed information for each iteration for diagnostics
!!
!! \author Met Office
!!
!! \date 29/03/2021: Created
!!
subroutine ufo_rttovonedvarcheck_PrintIterInfo(yob, hofx, channels, &
                                         guessprofile, backprofile, &
                                         diffprofile, binv, hmatrix)

implicit none

real(kind_real), intent(in)       :: yob(:)
real(kind_real), intent(in)       :: hofx(:)
integer, intent(in)               :: channels(:)
real(kind_real), intent(in)       :: guessprofile(:)
real(kind_real), intent(in)       :: backprofile(:)
real(kind_real), intent(in)       :: diffprofile(:)
real(kind_real), intent(in)       :: binv(:,:) ! (nprofelements,nprofelements)
real(kind_real), intent(in)       :: hmatrix(:,:) ! (nchans,nprofelements)

integer :: obs_size, profile_size, ii, jj
character(len=12) :: chans_fmt, prof_fmt
character(len=3) :: txt_nchans, txt_nprof
character(len=10) :: int_fmt

obs_size = size(yob)
write( unit=txt_nchans,fmt='(i3)' ) obs_size
write( unit=chans_fmt,fmt='(a)' ) '(' // trim(txt_nchans) // 'E30.16)'
write( unit=int_fmt,fmt='(a)' ) '(' // trim(txt_nchans) // 'I30)'

profile_size = size(guessprofile)
write( unit=txt_nprof,fmt='(i3)' ) profile_size
write( unit=prof_fmt,fmt='(a)' ) '(' // trim(txt_nprof) // 'E30.16)'

write(*,*) "Start print iter info"

! Print obs info
write(*,"(2A30)") "yob", "hofx"
do ii = 1, obs_size
  write(*,"(2E30.16)") yob(ii),hofx(ii)
end do

! Print profile info
write(*,"(3A30)") "guessprofile", "backprofile", "diffprofile"
do ii = 1, profile_size
  write(*,"(3E30.16)") guessprofile(ii),backprofile(ii), diffprofile(ii)
end do

! Print b inv
write(*,"(2A30)") "B inverse"
do ii = 1, profile_size
    write(*,prof_fmt) binv(:,ii)
end do

! Print hmatrix
write(*,"(2A30)") "hmatrix"
write(*, int_fmt) channels(:)
do ii = 1, profile_size
    write(*,chans_fmt) hmatrix(:,ii)
end do

write(*,*) "Finished print iter info"

end subroutine ufo_rttovonedvarcheck_PrintIterInfo

! ----------------------------------------------------------

end module ufo_rttovonedvarcheck_minimize_utils_mod
