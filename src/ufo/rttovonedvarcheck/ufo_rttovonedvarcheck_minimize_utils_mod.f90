! (C) Copyright 2020 Met Office
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module containing subroutines used by both minimizers.

module ufo_rttovonedvarcheck_minimize_utils_mod

use kinds
use ufo_constants_mod, only: grav, zero, t0c, half, one, two
use ufo_geovals_mod
use ufo_radiancerttov_tlad_mod
use ufo_rttovonedvarcheck_constants_mod
use ufo_rttovonedvarcheck_ob_mod
use ufo_rttovonedvarcheck_profindex_mod
use ufo_rttovonedvarcheck_rsubmatrix_mod

implicit none
private

! subroutines - public
public ufo_rttovonedvarcheck_GeoVaLs2ProfVec
public ufo_rttovonedvarcheck_ProfVec2GeoVaLs
public ufo_rttovonedvarcheck_check_geovals
public ufo_rttovonedvarcheck_CostFunction
public ufo_rttovonedvarcheck_Qsplit
public ufo_rttovonedvarcheck_CheckIteration
public ufo_rttovonedvarcheck_CheckCloudyIteration
public ufo_rttovonedvarcheck_Cholesky

! subroutines - private to the module
private ufo_rttovonedvarcheck_Qsat
private ufo_rttovonedvarcheck_QsatWat

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

!-------------------------------------------------------------------------------

write(*,*) trim(RoutineName)," start"

prof_x(:) = zero

!------------------------------
! 2. Set multi-level variables
!------------------------------

! Note that if number of profile levels is less than number of pressure levels
! we assume the levels are from the surface upwards (remember that RTTOV levels
! are upside-down)

! temperature - K
if (profindex % t(1) > 0) then
  varname = "air_temperature" ! K
  call ufo_geovals_get_var(geovals, varname, geoval)
  prof_x(profindex % t(1):profindex % t(2)) = geoval%vals(:, 1) ! K
end if

! specific_humidity - kg/kg - for retrieval is ln(g/kg)
if (profindex % q(1) > 0) then
  varname = "specific_humidity"  ! kg/kg
  call ufo_geovals_get_var(geovals, varname, geoval)
  prof_x(profindex % q(1):profindex % q(2)) = log (geoval%vals(:, 1) * 1000.0_kind_real) ! ln(g/kg)
end if

! specific_humidity - kg/kg - for retrieval is ln(g/kg)
if (profindex % qt(1) > 0) then
  nlevels = profindex%qt(2) - profindex%qt(1) + 1
  allocate(humidity_total(nlevels))
  humidity_total(:) = zero
  
  ! Get humidity data from geovals
  call ufo_geovals_get_var(geovals, "specific_humidity", geoval)
  humidity_total(:) = humidity_total(:) + geoval%vals(:, 1)
  call ufo_geovals_get_var(geovals, "mass_content_of_cloud_liquid_water_in_atmosphere_layer", geoval)
  humidity_total(:) = humidity_total(:) + geoval%vals(:, 1)
  call ufo_geovals_get_var(geovals, "mass_content_of_cloud_ice_in_atmosphere_layer", geoval)
  humidity_total(:) = humidity_total(:) + geoval%vals(:, 1)
  
  ! Convert from kg/kg to ln(g/kg)
  prof_x(profindex % qt(1):profindex % qt(2)) = log (humidity_total * 1000.0_kind_real) ! ln(g/kg)

  deallocate(humidity_total)
end if

!----------------------------
! 3. Single-valued variables
!----------------------------

! surface air temperature - K
if (profindex % t2 > 0) then
  varname = "air_temperature_at_two_meters_above_surface"
  call ufo_geovals_get_var(geovals, varname, geoval)
  prof_x(profindex % t2) = geoval%vals(1, 1)
end if

! surface specific humidity - kg/kg - for retrieval is ln(g/kg)
if (profindex % q2 > 0) then
  varname = "specific_humidity_at_two_meters_above_surface"
  call ufo_geovals_get_var(geovals, varname, geoval)
  prof_x(profindex % q2) = log (geoval%vals(1, 1) * 1000.0_kind_real) ! ln(g/kg)
end if

! surface pressure
if (profindex % pstar > 0) then
  varname = "air_pressure_at_two_meters_above_surface"
  call ufo_geovals_get_var(geovals, varname, geoval)
  prof_x(profindex % pstar) = geoval%vals(1, 1) / 100.0_kind_real  ! Pa to hPa
end if

! surface skin temperature - K
if (profindex % tstar > 0) then
  varname = "skin_temperature"
  call ufo_geovals_get_var(geovals, varname, geoval)
  prof_x(profindex % tstar) = geoval%vals(1, 1)
end if

! cloud top pressure
if (profindex % cloudtopp > 0) then
  prof_x(profindex % cloudtopp) = ob % cloudtopp ! carried around as hPa
end if

! cloud fraction
if (profindex % cloudfrac > 0) then
  prof_x(profindex % cloudfrac) = ob % cloudfrac
end if

! windspeed. Remember that all wind have been transferred to u and v is set to zero
! for windspeed retrievals
if (profindex % windspeed > 0) then
  varname = "eastward_wind"
  call ufo_geovals_get_var(geovals, varname, geoval)
  prof_x(profindex % windspeed) = geoval%vals(1, 1)
end if

!----------------------------
! 4. Emissivities
!----------------------------

! Microwave Emissivity
if (profindex % mwemiss(1) > 0) then
  ! Check that emissivity map is the correct size for the profile
  if ((profindex % mwemiss(2) - profindex % mwemiss(1) + 1) /= size(EmissMap)) then
    call abor1_ftn("mwemiss size differs from emissivity map")
  end if
  ! Copy microwave emissivity to profile
    prof_x(profindex % mwemiss(1):profindex % mwemiss(2)) = ob % emiss(EmissMap)
end if

! Retrieval of emissivity principal components
IF (profindex % emisspc(1) > 0) THEN
  ! convert ob % emiss to emiss pc using
  ! Ops_SatRad_PCToEmis
  ! Prof(profindex % emisspc(1):profindex % emisspc(2)) = emiss_pc
END IF

write(*,*) trim(RoutineName)," end"

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
                                            prof_x )        ! in

implicit none

! subroutine arguments:
type(ufo_geovals), intent(inout) :: geovals   !< model data at obs location
type(ufo_rttovonedvarcheck_profindex), intent(in) :: profindex !< index array for x vector
type(ufo_rttovonedvarcheck_ob), intent(inout) :: ob   !< satellite metadata
real(kind_real), intent(in)      :: prof_x(:) !< x vector

! Local arguments:
character(len=*), parameter  :: RoutineName = "ufo_rttovonedvarcheck_ProfVec2GeoVaLs"
character(len=max_string)    :: varname
integer                      :: gv_index, i, ii
integer                      :: nlevels
integer                      :: EmissElement
type(ufo_geoval), pointer    :: geoval
real(kind_real), allocatable :: temperature(:)
real(kind_real), allocatable :: pressure(:)
real(kind_real), allocatable :: humidity_total(:)
real(kind_real), allocatable :: q(:)
real(kind_real), allocatable :: ql(:)
real(kind_real), allocatable :: qi(:)

!-------------------------------------------------------------------------------

write(*,*) trim(RoutineName)," start"

!------------------------------
! 2. Set multi-level variables
!------------------------------

! Note that if number of profile levels is less than number of pressure levels
! we assume the levels are from the surface upwards (remember that RTTOV levels
! are upside-down)

! temperature - K
if (profindex % t(1) > 0) then
  varname = "air_temperature" ! K
  gv_index = 0
  do i=1,geovals%nvar
    if (varname == trim(geovals%variables(i))) gv_index = i
  end do
  geovals%geovals(gv_index)%vals(:,1) = prof_x(profindex % t(1):profindex % t(2)) ! K
end if

! specific_humidity - kg/kg - for retrieval is ln(g/kg)
if (profindex % q(1) > 0) then
  varname = "specific_humidity"  ! kg/kg
  gv_index = 0
  do i=1,geovals%nvar
    if (varname == trim(geovals%variables(i))) gv_index = i
  end do
  geovals%geovals(gv_index)%vals(:,1) = EXP (prof_x(profindex % q(1):profindex % q(2))) / &
                                                  1000.0_kind_real ! ln(g/kg) => kg/kg
end if

! specific_humidity - kg/kg - for retrieval is ln(g/kg)
if (profindex % qt(1) > 0) then
  nlevels = profindex%qt(2) - profindex%qt(1) + 1
  allocate(temperature(nlevels))
  allocate(pressure(nlevels))
  allocate(humidity_total(nlevels))
  allocate(q(nlevels))
  allocate(ql(nlevels))
  allocate(qi(nlevels))
  
  ! Convert from ln(g/kg) to kg/kg
  humidity_total(:) = EXP (prof_x(profindex % qt(1):profindex % qt(2))) / &
                                1000.0_kind_real ! ln(g/kg) => kg/kg

  ! Get temperature and pressure from geovals
  call ufo_geovals_get_var(geovals, "air_temperature", geoval)
  temperature(:) = geoval%vals(:, 1) ! K
  call ufo_geovals_get_var(geovals, "air_pressure", geoval)
  pressure(:) = geoval%vals(:, 1)    ! Pa

  ! Split qtotal to q(water_vapour), q(liquid), q(ice)
  call ufo_rttovonedvarcheck_Qsplit (1,      & ! in
                          temperature(:),    & ! in
                          pressure(:),       & ! in
                          nlevels,           & ! in
                          humidity_total(:), & ! in
                          q(:),              & ! out
                          ql(:),             & ! out
                          qi(:))               ! out

  ! Assign values to geovals
  varname = "specific_humidity"  ! kg/kg
  gv_index = 0
  do i=1,geovals%nvar
    if (varname == trim(geovals%variables(i))) gv_index = i
  end do
  geovals%geovals(gv_index)%vals(:,1) = q(:)

  varname = "mass_content_of_cloud_liquid_water_in_atmosphere_layer"  ! kg/kg
  gv_index = 0
  do i=1,geovals%nvar
    if (varname == trim(geovals%variables(i))) gv_index = i
  end do
  geovals%geovals(gv_index)%vals(:,1) = ql(:)

  varname = "mass_content_of_cloud_ice_in_atmosphere_layer"  ! kg/kg
  gv_index = 0
  do i=1,geovals%nvar
    if (varname == trim(geovals%variables(i))) gv_index = i
  end do
  geovals%geovals(gv_index)%vals(:,1) = qi(:)

  deallocate(temperature)
  deallocate(pressure)
  deallocate(humidity_total)
  deallocate(q)
  deallocate(ql)
  deallocate(qi)

end if

!----------------------------
! 3. Single-valued variables
!----------------------------

! surface air temperature - K
if (profindex % t2 > 0) then
  varname = "air_temperature_at_two_meters_above_surface"
  gv_index = 0
  do i=1,geovals%nvar
    if (varname == trim(geovals%variables(i))) gv_index = i
  end do
  geovals%geovals(gv_index)%vals(1,1) = prof_x(profindex % t2) ! K
end if

! surface specific humidity - kg/kg - for retrieval is ln(g/kg)
if (profindex % q2 > 0) then
  varname = "specific_humidity_at_two_meters_above_surface"
  gv_index = 0
  do i=1,geovals%nvar
    if (varname == trim(geovals%variables(i))) gv_index = i
  end do
  geovals%geovals(gv_index)%vals(1,1) = EXP (prof_x(profindex % q2)) / 1000.0_kind_real ! ln(g/kg) => kg/kg
end if

! surface pressure
if (profindex % pstar > 0) then
  varname = "air_pressure_at_two_meters_above_surface"
  gv_index = 0
  do i=1,geovals%nvar
    if (varname == trim(geovals%variables(i))) gv_index = i
  end do
  geovals%geovals(gv_index)%vals(1,1) = prof_x(profindex % pstar) * 100.0_kind_real
end if

! surface skin temperature - K
if (profindex % tstar > 0) then
  varname = "skin_temperature"
  gv_index = 0
  do i=1,geovals%nvar
    if (varname == trim(geovals%variables(i))) gv_index = i
  end do
  geovals%geovals(gv_index)%vals(1,1) = prof_x(profindex % tstar)
end if

! cloud top pressure - passed through via the ob
if (profindex % cloudtopp > 0) then
  ob % cloudtopp = prof_x(profindex % cloudtopp) ! stored in ob as hPa
end if

! cloud fraction - passed through via the ob
if (profindex % cloudfrac > 0) then
  ob % cloudfrac = prof_x(profindex % cloudfrac)
end if

! windspeed
IF (profindex % windspeed > 0) THEN
  ! Remember that we transfer all wind to u and set v to zero for
  ! windspeed retrieval.
  varname = "eastward_wind"
  gv_index = 0
  do i=1,geovals%nvar
    if (varname == trim(geovals%variables(i))) gv_index = i
  end do
  geovals%geovals(gv_index)%vals(1,1) = prof_x(profindex % windspeed)
END IF

!----------------------------
! 4. Emissivities
!------------------------

! Retrieval of microwave emissivity directly
if (profindex % mwemiss(1) > 0) THEN
  do ii = 1, size(ob % channels_used)
    EmissElement = EmissElements(ob % channels_used(ii))
    ob % emiss(ii) = prof_x(profindex % mwemiss(1) + EmissElement - 1)
  end do
end if

! Retrieval of emissivity principal components
IF (profindex % emisspc(1) > 0) THEN
  ! emiss_pc = prof_x(profindex % emisspc(1):profindex % emisspc(2)) 
  ! convert emiss_pc to ob % emissivity using
  ! Ops_SatRad_EmisToPC
END IF

write(*,*) trim(RoutineName)," end"

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
subroutine ufo_rttovonedvarcheck_check_geovals(geovals, profindex)

implicit none

! subroutine arguments:
type(ufo_geovals), intent(inout) :: geovals   !< model data at obs location
type(ufo_rttovonedvarcheck_profindex), intent(in) :: profindex !< index array for x vector

character(len=*), parameter  :: routinename = "ufo_rttovonedvarcheck_check_geovals"
character(len=max_string)    :: varname
type(ufo_geoval), pointer    :: geoval
integer                      :: gv_index, i     ! counters
integer                      :: nlevels
real(kind_real), allocatable :: temperature(:)  ! temperature (K)
real(kind_real), allocatable :: pressure(:)     ! pressure (Pa)
real(kind_real), allocatable :: humidity_total(:)
real(kind_real), allocatable :: q(:)            ! specific humidity (kg/kg)
real(kind_real), allocatable :: ql(:)
real(kind_real), allocatable :: qi(:)
real(kind_real)              :: u_wind
real(kind_real)              :: v_wind
real(kind_real)              :: new_u_wind
real(kind_real)              :: new_v_wind

write(*,*) routinename, " : started"

!-------------------------
! 1. Specific humidity total
!-------------------------

if (profindex % qt(1) > 0) then

  nlevels = profindex % qt(2) - profindex % qt(1) + 1
  allocate(temperature(nlevels))
  allocate(pressure(nlevels))
  allocate(humidity_total(nlevels))
  allocate(q(nlevels))
  allocate(ql(nlevels))
  allocate(qi(nlevels))

  ! Get temperature and pressure from geovals
  call ufo_geovals_get_var(geovals, "air_temperature", geoval)
  temperature(:) = geoval%vals(:, 1) ! K
  call ufo_geovals_get_var(geovals, "air_pressure", geoval)
  pressure(:) = geoval%vals(:, 1)    ! Pa

  ! Get humidity data from geovals
  humidity_total(:) = 0.0
  call ufo_geovals_get_var(geovals, "specific_humidity", geoval)
  humidity_total(:) = humidity_total(:) + geoval%vals(:, 1)
  call ufo_geovals_get_var(geovals, "mass_content_of_cloud_liquid_water_in_atmosphere_layer", geoval)
  humidity_total(:) = humidity_total(:) + geoval%vals(:, 1)

  ! Split qtotal to q(water_vapour), q(liquid), q(ice)
  call ufo_rttovonedvarcheck_Qsplit (1,      & ! in
                          temperature(:),    & ! in
                          pressure(:),       & ! in
                          nlevels,           & ! in
                          humidity_total(:), & ! in
                          q(:),              & ! out
                          ql(:),             & ! out
                          qi(:))               ! out

  ! Assign values to geovals
  varname = "specific_humidity"  ! kg/kg
  gv_index = 0
  do i=1,geovals%nvar
    if (varname == trim(geovals%variables(i))) gv_index = i
  end do
  geovals%geovals(gv_index) % vals(:,1) = q(:)

  varname = "mass_content_of_cloud_liquid_water_in_atmosphere_layer"  ! kg/kg
  gv_index = 0
  do i=1,geovals%nvar
    if (varname == trim(geovals%variables(i))) gv_index = i
  end do
  geovals%geovals(gv_index) % vals(:,1) = ql(:)

  varname = "mass_content_of_cloud_ice_in_atmosphere_layer"  ! kg/kg
  gv_index = 0
  do i=1,geovals%nvar
    if (varname == trim(geovals%variables(i))) gv_index = i
  end do
  geovals%geovals(gv_index)%vals(:,1) = qi(:)

  deallocate(temperature)
  deallocate(pressure)
  deallocate(humidity_total)
  deallocate(q)
  deallocate(ql)
  deallocate(qi)

end if

!-------
! 2. Wind
!-------

! RTTOV is isotropic, therefore if we only want to retrieve a "total" windspeed,
! with no directional information, we can put all the wind into u and set v to
! zero. If we are not retrieving windspeed, we just leave u and v separate to
! avoid confusion.

if (profindex % windspeed > 0) THEN
  ! Get winds from geovals
  varname = "eastward_wind"  ! m/s
  call ufo_geovals_get_var(geovals, varname, geoval)
  u_wind = geoval%vals(1,1)

  varname = "northward_wind"  ! m/s
  call ufo_geovals_get_var(geovals, varname, geoval)
  v_wind = geoval%vals(1,1)

  ! Convert to "total" windspeed
  new_u_wind = sqrt(u_wind * u_wind + v_wind * v_wind)
  new_v_wind = zero

  ! Write back to geovals
  varname = "eastward_wind"  ! m/s
  gv_index = 0
  do i=1,geovals%nvar
    if (varname == trim(geovals%variables(i))) gv_index = i
  end do
  geovals%geovals(gv_index)%vals(1,1) = new_u_wind

  varname = "northward_wind"  ! m/s
  gv_index = 0
  do i=1,geovals%nvar
    if (varname == trim(geovals%variables(i))) gv_index = i
  end do
  geovals%geovals(gv_index)%vals(1,1) = new_v_wind

end if

! Tidy up
if (allocated(temperature))    deallocate(temperature)
if (allocated(pressure))       deallocate(pressure)
if (allocated(humidity_total)) deallocate(humidity_total)
if (allocated(q))              deallocate(q)
if (allocated(ql))             deallocate(ql)
if (allocated(qi))             deallocate(qi)

write(*,*) routinename, " : ended"

end subroutine ufo_rttovonedvarcheck_check_geovals

!-------------------------------------------------------------------------------
!> Calculate the cost function.
!!
!! \details Heritage: Ops_SatRad_CostFunction.f90
!!
!! Caculate the cost function from the input delta's and error
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

write(*,*) trim(RoutineName)," start"

allocate(RinvDeltaY, source=DeltaObs)

y_size = size(DeltaObs)
Jcost(:) = zero
call r_matrix % multiply_inverse_vector(DeltaObs, RinvDeltaY)

Jo = half * dot_product(DeltaObs, RinvDeltaY)
Jb = half * dot_product(DeltaProf, (matmul(b_inv, DeltaProf)))
Jcurrent = Jb + Jo

Jcost(1) = (Jo + Jb) * two / real (y_size)     ! Normalize cost by nchans
Jcost(2) = Jb * two / real (y_size)            ! Normalize cost by nchans
Jcost(3) = Jo * two / real (y_size)            ! Normalize cost by nchans

write(*,*) "Jo, Jb, Jcurrent = ", Jcost(3), Jcost(2), Jcost(1)

deallocate(RinvDeltaY)

end subroutine ufo_rttovonedvarcheck_CostFunction

!-------------------------------------------------------------------------------
!> Split the humidity into water vapour, liquid water and ice.
!!
!! \details Heritage: Ops_SatRad_Qsplit.f90
!!
!! if output_type=1 : Split total water content (qtotal) into <br>
!!   water vapor content (q) and <br>
!!   cloud liquid water content (ql) and <br>
!!   cloud ice water content (qi) <br>
!!
!! if output_type ne 1 : Compute derivatives: (q) =dq/dqtotal <br>
!!                                            (ql)=dql/dqtotal <br>
!!                                            (qi)=dqi/dqtotal <br>
!!
!! \warning The derivatives are not valid if LtemperatureVar=.true. since
!!     qsaturated depends on temperature.
!!
!! The partitioning of the excess moisture between ice and clw uses a temperature
!! based parametrization based on aircraft data Ref: Jones DC Reading phdthesis
!! p126
!!
!! \author Met Office
!!
!! \date 09/06/2020: Created
!!
subroutine ufo_rttovonedvarcheck_Qsplit (output_type, &
                              t,           &
                              p,           &
                              nlevels_q,   &
                              qtotal,      &
                              q,           &
                              ql,          &
                              qi)

implicit none

! subroutine arguments:
integer, intent(in)               :: output_type       !< output profiles or gradients
real(kind=kind_real), intent(in)  :: t(:)              !< air temperature
real(kind=kind_real), intent(in)  :: p(:)              !< air pressure
integer, intent(in)               :: nlevels_q         !< no. of levels
real(kind=kind_real), intent(in)  :: qtotal(nlevels_q) !< humidity total (kg/kg)
real(kind=kind_real), intent(out) :: q(nlevels_q)      !< water vapour component q
real(kind=kind_real), intent(out) :: ql(nlevels_q)     !< liquid component ql
real(kind=kind_real), intent(out) :: qi(nlevels_q)     !< ice component qi

! Local declarations:
integer                     :: nlevels_diff
integer                     :: i
integer                     :: toplevel_q
integer                     :: nlevels_mwclw
integer                     :: nlevels_1dvar
integer                     :: Qsplit_MixPhaseParam = 1
real(kind=kind_real), parameter :: lower_rh = 0.95
real(kind=kind_real), parameter :: upper_rh = 1.05
real(kind=kind_real), parameter :: Split_Factor = 0.5
real(kind=kind_real), parameter :: minTempQl = 233.15     ! temperature (K) below which all cloud is ice
real(kind=kind_real), parameter :: min_q = 3.0E-6         ! ( kg / kg )
real(kind=kind_real) :: qsaturated(nlevels_q)
real(kind=kind_real) :: RH_qtotal(nlevels_q)
real(kind=kind_real) :: qnv(nlevels_q)         ! non vapour component
real(kind=kind_real) :: qc(nlevels_q)          ! cloud component
real(kind=kind_real) :: V1(nlevels_q)
real(kind=kind_real) :: V2(nlevels_q)
real(kind=kind_real) :: V1zero
real(kind=kind_real) :: V2zero
real(kind=kind_real) :: W(nlevels_q)
real(kind=kind_real) :: Y1,Y2,Y3,Y4
real(kind=kind_real) :: intConst
real(kind=kind_real) :: Aconst
real(kind=kind_real) :: Bconst
real(kind=kind_real) :: Cconst
real(kind=kind_real) :: Dconst
real(kind=kind_real) :: Denom
real(kind=kind_real) :: SmallValue
real(kind=kind_real) :: LF(nlevels_q)          ! fraction of ql to ql+qi
real(kind=kind_real) :: QsplitRainParamA
real(kind=kind_real) :: QsplitRainParamB
real(kind=kind_real) :: QsplitRainParamC
character(len=*), parameter :: RoutineName = "ufo_rttovonedvarcheck_Qsplit"
character(len=max_string)   :: message
logical                     :: useQtsplitRain

! Addition to make it work
nlevels_1dvar = nlevels_q
nlevels_mwclw = nlevels_q
useQtsplitRain = .true.
QsplitRainParamA = 0.15_kind_real
QsplitRainParamB = 0.09_kind_real
QsplitRainParamC = 50.0_kind_real

! Calculate the highest wet level

!toplevel_q = nlevels_1dvar - nlevels_q + 1 ! assuming they are the same size now change from OPS

! Compute saturated water vapor profile for nlevels_q only

call ufo_rttovonedvarcheck_Qsat (qsaturated(1:nlevels_q), & ! out
                           t(1:nlevels_q),                & ! in
                           p(1:nlevels_q),                & ! in
                           nlevels_q)                       ! in

SmallValue = one / 8.5_kind_real
Denom = SmallValue * (upper_rh - lower_rh)

! don't let rh exceed 2.0 to avoid cosh function blowing up
RH_qtotal(:) = min (qtotal(:) / qsaturated(:), two)

V1(:) = (RH_qtotal(:) - lower_rh) / Denom
V2(:) = (RH_qtotal(:) - upper_rh) / Denom

Y1 = one
Y2 = Split_Factor
Y3 = zero
Y4 = Split_Factor

Aconst = (Y2 - Y1) / two
Bconst = -(Y4 - Y3) / two
Cconst = (Y2 + Y1) / two
Dconst = -(Y4 + Y3) / two

! Deal with mixed phase clouds
! Compute fraction of ql to ql+qi based on temperature profile

! Either
!Qsplit_MixPhaseParam=1 Jones Method

if (Qsplit_MixPhaseParam == 1) then

  where (t(:) - t0c >= -0.01_kind_real) ! -0.01degc and above

    ! all ql
    LF(:) = one

  end where

  where (t(:) <= minTempql)

    ! all qi
    LF(:) = zero

  end where

  where (t(:) > minTempql .and. t(:) - t0c < -0.01_kind_real)

    ! Jones' parametrization
    LF(:) = sqrt (-1.0_kind_real * log (-0.025_kind_real * (t(:) - t0c)) / 70.0_kind_real)

  end where

!Or
!Bower et al 1996 parameterisation for ls cloud. Follows UM partitioning

else if (Qsplit_MixPhaseParam == 2) then

  LF(:) = (one/9.0_kind_real)*(T(:)-t0c)+one
  
  where( LF(:) > one )
    LF(:) = one
  end where

  where( LF(:) < zero )
    LF(:) = zero
  end where

!incorrect value supplied bale out
else

  write(message,'(A,I10)') 'Invalid option for Qsplit_MixPhaseParam. Value=', Qsplit_MixPhaseParam
  call abor1_ftn(message)

end if

! finally set LF to 0.0 for the rttov levels on which clw jacobians are not
! calculated since nlevels_mwclw < nlevels_q

nlevels_diff = nlevels_q - nlevels_mwclw

if (nlevels_diff > 0) then

  LF(1:nlevels_diff) = zero

end if

V1zero = -1.0_kind_real * lower_rh / Denom
V2zero = -1.0_kind_real * upper_rh / Denom
intConst = -(Aconst * Denom * log (cosh (V1zero)) + Bconst * Denom * log (cosh (V2zero)))
W(:) = Aconst * Denom * log (cosh (V1(:))) + Bconst * Denom * log (cosh (V2(:))) + &
       (Cconst + Dconst) * RH_qtotal(:) + intConst

! store the components of qtotal
! ensuring that they are above lower limits

if (useQtsplitRain) then

  ! Split qtotal into q and qnv (non-vapour part - includes
  ! ql, qi, qr)

  ! Split qtotal into q, ql ,qi

  do i = 1, nlevels_q

    q(i) = max (W(i) * qsaturated(i), min_q)
    qnv(i) = max (qtotal(i) - q(i), zero)

    ! Split qnv into a cloud and precipitation part

    qc(i) = max (QsplitRainParamA * (QsplitRainParamB - (QsplitRainParamB / &
                                    ((QsplitRainParamC * qnv(i)) + one))), zero)

    ! Finally split non-precip part into liquid and ice

    ql(i) = max (LF(i) * qc(i), zero)
    qi(i) = max ((one - LF(i)) * (qc(i)), zero)

  end do

else
  do i = 1, nlevels_q

    q(i) = max (W(i) * qsaturated(i), min_q)
    ql(i) = max (LF(i) * (qtotal(i) - q(i)), zero)
    qi(i) = max ((one - LF(i)) * (qtotal(i) - q(i)), zero)

  end do

end if

! Values of q, ql and qi are overwritten if output_type /= 1
! and replaced with the derivatives

if (output_type /= 1) then

  ! Compute derivates
  ! q = dq/dqtotal, ql = dql/dqtotal, qi=dqi/dqtotal

  q(:) = Aconst * tanh (V1(:)) + Cconst + Bconst * tanh (V2(:)) + Dconst

  if (useQtsplitRain) then

    ql(:) = LF(:) * QsplitRainParamA * QsplitRainParamB * QsplitRainParamC * (one - q(:)) /  &
                         ((QsplitRainParamC * qnv(:)) + one) ** 2
    qi(:) = (1.0 - LF(:)) * QsplitRainParamA * QsplitRainParamB * QsplitRainParamC * (one - q(:)) / &
                         ((QsplitRainParamC * qnv(:)) + one) ** 2

  else

    ql(:) = LF(:) * (one - q(:))
    qi(:) = (one - LF(:)) * (one - q(:))

  end if

end if

end subroutine ufo_rttovonedvarcheck_Qsplit

!-------------------------------------------------------------------------------
!> Calculate the Saturation Specific Humidity Scheme (Qsat): Vapour to Liquid/Ice.
!!
!! \details Heritage: Ops_Qsat.inc
!!
!! Returns a saturation mixing ratio given a temperature and pressure
!! using saturation vapour pressures caluclated using the Goff-Gratch
!! formulae, adopted by the WMO as taken from Landolt-Bornstein, 1987
!! Numerical data and Functional relationships in Science and
!! Technology.  Group V/Vol 4B Meteorology.  Physical and Chemical
!! properties of Air, P35.
!!
!! Value in the lookup table are over water above 0 degrees C and over
!! ice below this temperatures.
!!
!! Method: <br>
!!   uses lookup tables to find eSAT, calculates qSAT directly from that.
!!
!! \author Met Office
!!
!! \date 09/06/2020: Created
!!

subroutine ufo_rttovonedvarcheck_Qsat (QS, &
                                       T,  &
                                       P,  &
                                       npnts)

implicit none

! subroutine arguments:
integer, intent(in)                 :: npnts     !< Points being processed by qSAT scheme.
real(kind=kind_real), intent(in)    :: T(npnts)  !< Temperature (K)
real(kind=kind_real), intent(in)    :: P(npnts)  !< Pressure (Pa).
real(kind=kind_real), intent(inout) :: QS(npnts) !< Saturation mixing ratio (KG/KG)

! Local declarations:
real(kind=kind_real), parameter :: one_minus_epsilon = one - epsilon
real(kind=kind_real), parameter :: T_low = 183.15_kind_real  ! Lowest temperature for which look-up table is valid
real(kind=kind_real), parameter :: T_high = 338.15_kind_real  ! Highest temperature for which look-up table is valid
real(kind=kind_real), parameter :: delta_T = 0.1_kind_real    ! Temperature increment of look-up table
integer, parameter   :: N = ((T_high - T_low + (delta_T * 0.5_kind_real)) / delta_T) + one ! Size of lookup-table (gives 1551)
integer              :: ITABLE
real(kind=kind_real) :: ATABLE
real(kind=kind_real) :: FSUBW      ! Converts from sat vapour pressure in pure water to pressure in air
real(kind=kind_real) :: TT
integer              :: I
integer              :: IES
real(kind=kind_real) :: ES(0:N + 1)    ! Table of saturation water vapour pressure (PA)
character(len=*), parameter :: RoutineName = "ufo_rttovonedvarcheck_Qsat"

! Note: 0 element is a repeat of 1st element to cater for special case
!       of low temperatures (.LE.T_low) for which the array index is
!       rounded down due to machine precision.

data (ES(IES), IES= 0, 95) / 0.966483E-02, &
0.966483E-02,0.984279E-02,0.100240E-01,0.102082E-01,0.103957E-01, &
0.105865E-01,0.107803E-01,0.109777E-01,0.111784E-01,0.113825E-01, &
0.115902E-01,0.118016E-01,0.120164E-01,0.122348E-01,0.124572E-01, &
0.126831E-01,0.129132E-01,0.131470E-01,0.133846E-01,0.136264E-01, &
0.138724E-01,0.141225E-01,0.143771E-01,0.146356E-01,0.148985E-01, &
0.151661E-01,0.154379E-01,0.157145E-01,0.159958E-01,0.162817E-01, &
0.165725E-01,0.168680E-01,0.171684E-01,0.174742E-01,0.177847E-01, &
0.181008E-01,0.184216E-01,0.187481E-01,0.190801E-01,0.194175E-01, &
0.197608E-01,0.201094E-01,0.204637E-01,0.208242E-01,0.211906E-01, &
0.215631E-01,0.219416E-01,0.223263E-01,0.227172E-01,0.231146E-01, &
0.235188E-01,0.239296E-01,0.243465E-01,0.247708E-01,0.252019E-01, &
0.256405E-01,0.260857E-01,0.265385E-01,0.269979E-01,0.274656E-01, &
0.279405E-01,0.284232E-01,0.289142E-01,0.294124E-01,0.299192E-01, &
0.304341E-01,0.309571E-01,0.314886E-01,0.320285E-01,0.325769E-01, &
0.331348E-01,0.337014E-01,0.342771E-01,0.348618E-01,0.354557E-01, &
0.360598E-01,0.366727E-01,0.372958E-01,0.379289E-01,0.385717E-01, &
0.392248E-01,0.398889E-01,0.405633E-01,0.412474E-01,0.419430E-01, &
0.426505E-01,0.433678E-01,0.440974E-01,0.448374E-01,0.455896E-01, &
0.463545E-01,0.471303E-01,0.479191E-01,0.487190E-01,0.495322E-01/
      data (ES(IES),IES= 96,190) / &
0.503591E-01,0.511977E-01,0.520490E-01,0.529145E-01,0.537931E-01, &
0.546854E-01,0.555924E-01,0.565119E-01,0.574467E-01,0.583959E-01, &
0.593592E-01,0.603387E-01,0.613316E-01,0.623409E-01,0.633655E-01, &
0.644053E-01,0.654624E-01,0.665358E-01,0.676233E-01,0.687302E-01, &
0.698524E-01,0.709929E-01,0.721490E-01,0.733238E-01,0.745180E-01, &
0.757281E-01,0.769578E-01,0.782061E-01,0.794728E-01,0.807583E-01, &
0.820647E-01,0.833905E-01,0.847358E-01,0.861028E-01,0.874882E-01, &
0.888957E-01,0.903243E-01,0.917736E-01,0.932464E-01,0.947407E-01, &
0.962571E-01,0.977955E-01,0.993584E-01,0.100942E+00,0.102551E+00, &
0.104186E+00,0.105842E+00,0.107524E+00,0.109231E+00,0.110963E+00, &
0.112722E+00,0.114506E+00,0.116317E+00,0.118153E+00,0.120019E+00, &
0.121911E+00,0.123831E+00,0.125778E+00,0.127755E+00,0.129761E+00, &
0.131796E+00,0.133863E+00,0.135956E+00,0.138082E+00,0.140241E+00, &
0.142428E+00,0.144649E+00,0.146902E+00,0.149190E+00,0.151506E+00, &
0.153859E+00,0.156245E+00,0.158669E+00,0.161126E+00,0.163618E+00, &
0.166145E+00,0.168711E+00,0.171313E+00,0.173951E+00,0.176626E+00, &
0.179342E+00,0.182096E+00,0.184893E+00,0.187724E+00,0.190600E+00, &
0.193518E+00,0.196473E+00,0.199474E+00,0.202516E+00,0.205604E+00, &
0.208730E+00,0.211905E+00,0.215127E+00,0.218389E+00,0.221701E+00/
      data (ES(IES),IES=191,285) / &
0.225063E+00,0.228466E+00,0.231920E+00,0.235421E+00,0.238976E+00, &
0.242580E+00,0.246232E+00,0.249933E+00,0.253691E+00,0.257499E+00, &
0.261359E+00,0.265278E+00,0.269249E+00,0.273274E+00,0.277358E+00, &
0.281498E+00,0.285694E+00,0.289952E+00,0.294268E+00,0.298641E+00, &
0.303078E+00,0.307577E+00,0.312135E+00,0.316753E+00,0.321440E+00, &
0.326196E+00,0.331009E+00,0.335893E+00,0.340842E+00,0.345863E+00, &
0.350951E+00,0.356106E+00,0.361337E+00,0.366636E+00,0.372006E+00, &
0.377447E+00,0.382966E+00,0.388567E+00,0.394233E+00,0.399981E+00, &
0.405806E+00,0.411714E+00,0.417699E+00,0.423772E+00,0.429914E+00, &
0.436145E+00,0.442468E+00,0.448862E+00,0.455359E+00,0.461930E+00, &
0.468596E+00,0.475348E+00,0.482186E+00,0.489124E+00,0.496160E+00, &
0.503278E+00,0.510497E+00,0.517808E+00,0.525224E+00,0.532737E+00, &
0.540355E+00,0.548059E+00,0.555886E+00,0.563797E+00,0.571825E+00, &
0.579952E+00,0.588198E+00,0.596545E+00,0.605000E+00,0.613572E+00, &
0.622255E+00,0.631059E+00,0.639962E+00,0.649003E+00,0.658144E+00, &
0.667414E+00,0.676815E+00,0.686317E+00,0.695956E+00,0.705728E+00, &
0.715622E+00,0.725641E+00,0.735799E+00,0.746082E+00,0.756495E+00, &
0.767052E+00,0.777741E+00,0.788576E+00,0.799549E+00,0.810656E+00, &
0.821914E+00,0.833314E+00,0.844854E+00,0.856555E+00,0.868415E+00/
      data (ES(IES),IES=286,380) / &
0.880404E+00,0.892575E+00,0.904877E+00,0.917350E+00,0.929974E+00, &
0.942771E+00,0.955724E+00,0.968837E+00,0.982127E+00,0.995600E+00, &
0.100921E+01,0.102304E+01,0.103700E+01,0.105116E+01,0.106549E+01, &
0.108002E+01,0.109471E+01,0.110962E+01,0.112469E+01,0.113995E+01, &
0.115542E+01,0.117107E+01,0.118693E+01,0.120298E+01,0.121923E+01, &
0.123569E+01,0.125234E+01,0.126923E+01,0.128631E+01,0.130362E+01, &
0.132114E+01,0.133887E+01,0.135683E+01,0.137500E+01,0.139342E+01, &
0.141205E+01,0.143091E+01,0.145000E+01,0.146933E+01,0.148892E+01, &
0.150874E+01,0.152881E+01,0.154912E+01,0.156970E+01,0.159049E+01, &
0.161159E+01,0.163293E+01,0.165452E+01,0.167640E+01,0.169852E+01, &
0.172091E+01,0.174359E+01,0.176653E+01,0.178977E+01,0.181332E+01, &
0.183709E+01,0.186119E+01,0.188559E+01,0.191028E+01,0.193524E+01, &
0.196054E+01,0.198616E+01,0.201208E+01,0.203829E+01,0.206485E+01, &
0.209170E+01,0.211885E+01,0.214637E+01,0.217424E+01,0.220242E+01, &
0.223092E+01,0.225979E+01,0.228899E+01,0.231855E+01,0.234845E+01, &
0.237874E+01,0.240937E+01,0.244040E+01,0.247176E+01,0.250349E+01, &
0.253560E+01,0.256814E+01,0.260099E+01,0.263431E+01,0.266800E+01, &
0.270207E+01,0.273656E+01,0.277145E+01,0.280671E+01,0.284248E+01, &
0.287859E+01,0.291516E+01,0.295219E+01,0.298962E+01,0.302746E+01/
      data (ES(IES),IES=381,475) / &
0.306579E+01,0.310454E+01,0.314377E+01,0.318351E+01,0.322360E+01, &
0.326427E+01,0.330538E+01,0.334694E+01,0.338894E+01,0.343155E+01, &
0.347456E+01,0.351809E+01,0.356216E+01,0.360673E+01,0.365184E+01, &
0.369744E+01,0.374352E+01,0.379018E+01,0.383743E+01,0.388518E+01, &
0.393344E+01,0.398230E+01,0.403177E+01,0.408175E+01,0.413229E+01, &
0.418343E+01,0.423514E+01,0.428746E+01,0.434034E+01,0.439389E+01, &
0.444808E+01,0.450276E+01,0.455820E+01,0.461423E+01,0.467084E+01, &
0.472816E+01,0.478607E+01,0.484468E+01,0.490393E+01,0.496389E+01, &
0.502446E+01,0.508580E+01,0.514776E+01,0.521047E+01,0.527385E+01, &
0.533798E+01,0.540279E+01,0.546838E+01,0.553466E+01,0.560173E+01, &
0.566949E+01,0.573807E+01,0.580750E+01,0.587749E+01,0.594846E+01, &
0.602017E+01,0.609260E+01,0.616591E+01,0.623995E+01,0.631490E+01, &
0.639061E+01,0.646723E+01,0.654477E+01,0.662293E+01,0.670220E+01, &
0.678227E+01,0.686313E+01,0.694495E+01,0.702777E+01,0.711142E+01, &
0.719592E+01,0.728140E+01,0.736790E+01,0.745527E+01,0.754352E+01, &
0.763298E+01,0.772316E+01,0.781442E+01,0.790676E+01,0.800001E+01, &
0.809435E+01,0.818967E+01,0.828606E+01,0.838343E+01,0.848194E+01, &
0.858144E+01,0.868207E+01,0.878392E+01,0.888673E+01,0.899060E+01, &
0.909567E+01,0.920172E+01,0.930909E+01,0.941765E+01,0.952730E+01/
      data (ES(IES),IES=476,570) / &
0.963821E+01,0.975022E+01,0.986352E+01,0.997793E+01,0.100937E+02, &
0.102105E+02,0.103287E+02,0.104481E+02,0.105688E+02,0.106909E+02, &
0.108143E+02,0.109387E+02,0.110647E+02,0.111921E+02,0.113207E+02, &
0.114508E+02,0.115821E+02,0.117149E+02,0.118490E+02,0.119847E+02, &
0.121216E+02,0.122601E+02,0.124002E+02,0.125416E+02,0.126846E+02, &
0.128290E+02,0.129747E+02,0.131224E+02,0.132712E+02,0.134220E+02, &
0.135742E+02,0.137278E+02,0.138831E+02,0.140403E+02,0.141989E+02, &
0.143589E+02,0.145211E+02,0.146845E+02,0.148501E+02,0.150172E+02, &
0.151858E+02,0.153564E+02,0.155288E+02,0.157029E+02,0.158786E+02, &
0.160562E+02,0.162358E+02,0.164174E+02,0.166004E+02,0.167858E+02, &
0.169728E+02,0.171620E+02,0.173528E+02,0.175455E+02,0.177406E+02, &
0.179372E+02,0.181363E+02,0.183372E+02,0.185400E+02,0.187453E+02, &
0.189523E+02,0.191613E+02,0.193728E+02,0.195866E+02,0.198024E+02, &
0.200200E+02,0.202401E+02,0.204626E+02,0.206871E+02,0.209140E+02, &
0.211430E+02,0.213744E+02,0.216085E+02,0.218446E+02,0.220828E+02, &
0.223241E+02,0.225671E+02,0.228132E+02,0.230615E+02,0.233120E+02, &
0.235651E+02,0.238211E+02,0.240794E+02,0.243404E+02,0.246042E+02, &
0.248704E+02,0.251390E+02,0.254109E+02,0.256847E+02,0.259620E+02, &
0.262418E+02,0.265240E+02,0.268092E+02,0.270975E+02,0.273883E+02/
      data (ES(IES),IES=571,665) / &
0.276822E+02,0.279792E+02,0.282789E+02,0.285812E+02,0.288867E+02, &
0.291954E+02,0.295075E+02,0.298222E+02,0.301398E+02,0.304606E+02, &
0.307848E+02,0.311119E+02,0.314424E+02,0.317763E+02,0.321133E+02, &
0.324536E+02,0.327971E+02,0.331440E+02,0.334940E+02,0.338475E+02, &
0.342050E+02,0.345654E+02,0.349295E+02,0.352975E+02,0.356687E+02, &
0.360430E+02,0.364221E+02,0.368042E+02,0.371896E+02,0.375790E+02, &
0.379725E+02,0.383692E+02,0.387702E+02,0.391744E+02,0.395839E+02, &
0.399958E+02,0.404118E+02,0.408325E+02,0.412574E+02,0.416858E+02, &
0.421188E+02,0.425551E+02,0.429962E+02,0.434407E+02,0.438910E+02, &
0.443439E+02,0.448024E+02,0.452648E+02,0.457308E+02,0.462018E+02, &
0.466775E+02,0.471582E+02,0.476428E+02,0.481313E+02,0.486249E+02, &
0.491235E+02,0.496272E+02,0.501349E+02,0.506479E+02,0.511652E+02, &
0.516876E+02,0.522142E+02,0.527474E+02,0.532836E+02,0.538266E+02, &
0.543737E+02,0.549254E+02,0.554839E+02,0.560456E+02,0.566142E+02, &
0.571872E+02,0.577662E+02,0.583498E+02,0.589392E+02,0.595347E+02, &
0.601346E+02,0.607410E+02,0.613519E+02,0.619689E+02,0.625922E+02, &
0.632204E+02,0.638550E+02,0.644959E+02,0.651418E+02,0.657942E+02, &
0.664516E+02,0.671158E+02,0.677864E+02,0.684624E+02,0.691451E+02, &
0.698345E+02,0.705293E+02,0.712312E+02,0.719398E+02,0.726542E+02/
      data (ES(IES),IES=666,760) / &
0.733754E+02,0.741022E+02,0.748363E+02,0.755777E+02,0.763247E+02, &
0.770791E+02,0.778394E+02,0.786088E+02,0.793824E+02,0.801653E+02, &
0.809542E+02,0.817509E+02,0.825536E+02,0.833643E+02,0.841828E+02, &
0.850076E+02,0.858405E+02,0.866797E+02,0.875289E+02,0.883827E+02, &
0.892467E+02,0.901172E+02,0.909962E+02,0.918818E+02,0.927760E+02, &
0.936790E+02,0.945887E+02,0.955071E+02,0.964346E+02,0.973689E+02, &
0.983123E+02,0.992648E+02,0.100224E+03,0.101193E+03,0.102169E+03, &
0.103155E+03,0.104150E+03,0.105152E+03,0.106164E+03,0.107186E+03, &
0.108217E+03,0.109256E+03,0.110303E+03,0.111362E+03,0.112429E+03, &
0.113503E+03,0.114588E+03,0.115684E+03,0.116789E+03,0.117903E+03, &
0.119028E+03,0.120160E+03,0.121306E+03,0.122460E+03,0.123623E+03, &
0.124796E+03,0.125981E+03,0.127174E+03,0.128381E+03,0.129594E+03, &
0.130822E+03,0.132058E+03,0.133306E+03,0.134563E+03,0.135828E+03, &
0.137109E+03,0.138402E+03,0.139700E+03,0.141017E+03,0.142338E+03, &
0.143676E+03,0.145025E+03,0.146382E+03,0.147753E+03,0.149133E+03, &
0.150529E+03,0.151935E+03,0.153351E+03,0.154783E+03,0.156222E+03, &
0.157678E+03,0.159148E+03,0.160624E+03,0.162117E+03,0.163621E+03, &
0.165142E+03,0.166674E+03,0.168212E+03,0.169772E+03,0.171340E+03, &
0.172921E+03,0.174522E+03,0.176129E+03,0.177755E+03,0.179388E+03/
      data (ES(IES),IES=761,855) / &
0.181040E+03,0.182707E+03,0.184382E+03,0.186076E+03,0.187782E+03, &
0.189503E+03,0.191240E+03,0.192989E+03,0.194758E+03,0.196535E+03, &
0.198332E+03,0.200141E+03,0.201963E+03,0.203805E+03,0.205656E+03, &
0.207532E+03,0.209416E+03,0.211317E+03,0.213236E+03,0.215167E+03, &
0.217121E+03,0.219087E+03,0.221067E+03,0.223064E+03,0.225080E+03, &
0.227113E+03,0.229160E+03,0.231221E+03,0.233305E+03,0.235403E+03, &
0.237520E+03,0.239655E+03,0.241805E+03,0.243979E+03,0.246163E+03, &
0.248365E+03,0.250593E+03,0.252830E+03,0.255093E+03,0.257364E+03, &
0.259667E+03,0.261979E+03,0.264312E+03,0.266666E+03,0.269034E+03, &
0.271430E+03,0.273841E+03,0.276268E+03,0.278722E+03,0.281185E+03, &
0.283677E+03,0.286190E+03,0.288714E+03,0.291266E+03,0.293834E+03, &
0.296431E+03,0.299045E+03,0.301676E+03,0.304329E+03,0.307006E+03, &
0.309706E+03,0.312423E+03,0.315165E+03,0.317930E+03,0.320705E+03, &
0.323519E+03,0.326350E+03,0.329199E+03,0.332073E+03,0.334973E+03, &
0.337897E+03,0.340839E+03,0.343800E+03,0.346794E+03,0.349806E+03, &
0.352845E+03,0.355918E+03,0.358994E+03,0.362112E+03,0.365242E+03, &
0.368407E+03,0.371599E+03,0.374802E+03,0.378042E+03,0.381293E+03, &
0.384588E+03,0.387904E+03,0.391239E+03,0.394604E+03,0.397988E+03, &
0.401411E+03,0.404862E+03,0.408326E+03,0.411829E+03,0.415352E+03/
      data (ES(IES),IES=856,950) / &
0.418906E+03,0.422490E+03,0.426095E+03,0.429740E+03,0.433398E+03, &
0.437097E+03,0.440827E+03,0.444570E+03,0.448354E+03,0.452160E+03, &
0.455999E+03,0.459870E+03,0.463765E+03,0.467702E+03,0.471652E+03, &
0.475646E+03,0.479674E+03,0.483715E+03,0.487811E+03,0.491911E+03, &
0.496065E+03,0.500244E+03,0.504448E+03,0.508698E+03,0.512961E+03, &
0.517282E+03,0.521617E+03,0.525989E+03,0.530397E+03,0.534831E+03, &
0.539313E+03,0.543821E+03,0.548355E+03,0.552938E+03,0.557549E+03, &
0.562197E+03,0.566884E+03,0.571598E+03,0.576351E+03,0.581131E+03, &
0.585963E+03,0.590835E+03,0.595722E+03,0.600663E+03,0.605631E+03, &
0.610641E+03,0.615151E+03,0.619625E+03,0.624140E+03,0.628671E+03, &
0.633243E+03,0.637845E+03,0.642465E+03,0.647126E+03,0.651806E+03, &
0.656527E+03,0.661279E+03,0.666049E+03,0.670861E+03,0.675692E+03, &
0.680566E+03,0.685471E+03,0.690396E+03,0.695363E+03,0.700350E+03, &
0.705381E+03,0.710444E+03,0.715527E+03,0.720654E+03,0.725801E+03, &
0.730994E+03,0.736219E+03,0.741465E+03,0.746756E+03,0.752068E+03, &
0.757426E+03,0.762819E+03,0.768231E+03,0.773692E+03,0.779172E+03, &
0.784701E+03,0.790265E+03,0.795849E+03,0.801483E+03,0.807137E+03, &
0.812842E+03,0.818582E+03,0.824343E+03,0.830153E+03,0.835987E+03, &
0.841871E+03,0.847791E+03,0.853733E+03,0.859727E+03,0.865743E+03/
      data (ES(IES),IES=951,1045) / &
0.871812E+03,0.877918E+03,0.884046E+03,0.890228E+03,0.896433E+03, &
0.902690E+03,0.908987E+03,0.915307E+03,0.921681E+03,0.928078E+03, &
0.934531E+03,0.941023E+03,0.947539E+03,0.954112E+03,0.960708E+03, &
0.967361E+03,0.974053E+03,0.980771E+03,0.987545E+03,0.994345E+03, &
0.100120E+04,0.100810E+04,0.101502E+04,0.102201E+04,0.102902E+04, &
0.103608E+04,0.104320E+04,0.105033E+04,0.105753E+04,0.106475E+04, &
0.107204E+04,0.107936E+04,0.108672E+04,0.109414E+04,0.110158E+04, &
0.110908E+04,0.111663E+04,0.112421E+04,0.113185E+04,0.113952E+04, &
0.114725E+04,0.115503E+04,0.116284E+04,0.117071E+04,0.117861E+04, &
0.118658E+04,0.119459E+04,0.120264E+04,0.121074E+04,0.121888E+04, &
0.122709E+04,0.123534E+04,0.124362E+04,0.125198E+04,0.126036E+04, &
0.126881E+04,0.127731E+04,0.128584E+04,0.129444E+04,0.130307E+04, &
0.131177E+04,0.132053E+04,0.132931E+04,0.133817E+04,0.134705E+04, &
0.135602E+04,0.136503E+04,0.137407E+04,0.138319E+04,0.139234E+04, &
0.140156E+04,0.141084E+04,0.142015E+04,0.142954E+04,0.143896E+04, &
0.144845E+04,0.145800E+04,0.146759E+04,0.147725E+04,0.148694E+04, &
0.149672E+04,0.150655E+04,0.151641E+04,0.152635E+04,0.153633E+04, &
0.154639E+04,0.155650E+04,0.156665E+04,0.157688E+04,0.158715E+04, &
0.159750E+04,0.160791E+04,0.161836E+04,0.162888E+04,0.163945E+04/
      data (ES(IES),IES=1046,1140) / &
0.165010E+04,0.166081E+04,0.167155E+04,0.168238E+04,0.169325E+04, &
0.170420E+04,0.171522E+04,0.172627E+04,0.173741E+04,0.174859E+04, &
0.175986E+04,0.177119E+04,0.178256E+04,0.179402E+04,0.180552E+04, &
0.181711E+04,0.182877E+04,0.184046E+04,0.185224E+04,0.186407E+04, &
0.187599E+04,0.188797E+04,0.190000E+04,0.191212E+04,0.192428E+04, &
0.193653E+04,0.194886E+04,0.196122E+04,0.197368E+04,0.198618E+04, &
0.199878E+04,0.201145E+04,0.202416E+04,0.203698E+04,0.204983E+04, &
0.206278E+04,0.207580E+04,0.208887E+04,0.210204E+04,0.211525E+04, &
0.212856E+04,0.214195E+04,0.215538E+04,0.216892E+04,0.218249E+04, &
0.219618E+04,0.220994E+04,0.222375E+04,0.223766E+04,0.225161E+04, &
0.226567E+04,0.227981E+04,0.229399E+04,0.230829E+04,0.232263E+04, &
0.233708E+04,0.235161E+04,0.236618E+04,0.238087E+04,0.239560E+04, &
0.241044E+04,0.242538E+04,0.244035E+04,0.245544E+04,0.247057E+04, &
0.248583E+04,0.250116E+04,0.251654E+04,0.253204E+04,0.254759E+04, &
0.256325E+04,0.257901E+04,0.259480E+04,0.261073E+04,0.262670E+04, &
0.264279E+04,0.265896E+04,0.267519E+04,0.269154E+04,0.270794E+04, &
0.272447E+04,0.274108E+04,0.275774E+04,0.277453E+04,0.279137E+04, &
0.280834E+04,0.282540E+04,0.284251E+04,0.285975E+04,0.287704E+04, &
0.289446E+04,0.291198E+04,0.292954E+04,0.294725E+04,0.296499E+04/
      data (ES(IES),IES=1141,1235) / &
0.298288E+04,0.300087E+04,0.301890E+04,0.303707E+04,0.305529E+04, &
0.307365E+04,0.309211E+04,0.311062E+04,0.312927E+04,0.314798E+04, &
0.316682E+04,0.318577E+04,0.320477E+04,0.322391E+04,0.324310E+04, &
0.326245E+04,0.328189E+04,0.330138E+04,0.332103E+04,0.334073E+04, &
0.336058E+04,0.338053E+04,0.340054E+04,0.342069E+04,0.344090E+04, &
0.346127E+04,0.348174E+04,0.350227E+04,0.352295E+04,0.354369E+04, &
0.356458E+04,0.358559E+04,0.360664E+04,0.362787E+04,0.364914E+04, &
0.367058E+04,0.369212E+04,0.371373E+04,0.373548E+04,0.375731E+04, &
0.377929E+04,0.380139E+04,0.382355E+04,0.384588E+04,0.386826E+04, &
0.389081E+04,0.391348E+04,0.393620E+04,0.395910E+04,0.398205E+04, &
0.400518E+04,0.402843E+04,0.405173E+04,0.407520E+04,0.409875E+04, &
0.412246E+04,0.414630E+04,0.417019E+04,0.419427E+04,0.421840E+04, &
0.424272E+04,0.426715E+04,0.429165E+04,0.431634E+04,0.434108E+04, &
0.436602E+04,0.439107E+04,0.441618E+04,0.444149E+04,0.446685E+04, &
0.449241E+04,0.451810E+04,0.454385E+04,0.456977E+04,0.459578E+04, &
0.462197E+04,0.464830E+04,0.467468E+04,0.470127E+04,0.472792E+04, &
0.475477E+04,0.478175E+04,0.480880E+04,0.483605E+04,0.486336E+04, &
0.489087E+04,0.491853E+04,0.494623E+04,0.497415E+04,0.500215E+04, &
0.503034E+04,0.505867E+04,0.508707E+04,0.511568E+04,0.514436E+04/
      data (ES(IES),IES=1236,1330) / &
0.517325E+04,0.520227E+04,0.523137E+04,0.526068E+04,0.529005E+04, &
0.531965E+04,0.534939E+04,0.537921E+04,0.540923E+04,0.543932E+04, &
0.546965E+04,0.550011E+04,0.553064E+04,0.556139E+04,0.559223E+04, &
0.562329E+04,0.565449E+04,0.568577E+04,0.571727E+04,0.574884E+04, &
0.578064E+04,0.581261E+04,0.584464E+04,0.587692E+04,0.590924E+04, &
0.594182E+04,0.597455E+04,0.600736E+04,0.604039E+04,0.607350E+04, &
0.610685E+04,0.614036E+04,0.617394E+04,0.620777E+04,0.624169E+04, &
0.627584E+04,0.631014E+04,0.634454E+04,0.637918E+04,0.641390E+04, &
0.644887E+04,0.648400E+04,0.651919E+04,0.655467E+04,0.659021E+04, &
0.662599E+04,0.666197E+04,0.669800E+04,0.673429E+04,0.677069E+04, &
0.680735E+04,0.684415E+04,0.688104E+04,0.691819E+04,0.695543E+04, &
0.699292E+04,0.703061E+04,0.706837E+04,0.710639E+04,0.714451E+04, &
0.718289E+04,0.722143E+04,0.726009E+04,0.729903E+04,0.733802E+04, &
0.737729E+04,0.741676E+04,0.745631E+04,0.749612E+04,0.753602E+04, &
0.757622E+04,0.761659E+04,0.765705E+04,0.769780E+04,0.773863E+04, &
0.777975E+04,0.782106E+04,0.786246E+04,0.790412E+04,0.794593E+04, &
0.798802E+04,0.803028E+04,0.807259E+04,0.811525E+04,0.815798E+04, &
0.820102E+04,0.824427E+04,0.828757E+04,0.833120E+04,0.837493E+04, &
0.841895E+04,0.846313E+04,0.850744E+04,0.855208E+04,0.859678E+04/
      data (ES(IES),IES=1331,1425) / &
0.864179E+04,0.868705E+04,0.873237E+04,0.877800E+04,0.882374E+04, &
0.886979E+04,0.891603E+04,0.896237E+04,0.900904E+04,0.905579E+04, &
0.910288E+04,0.915018E+04,0.919758E+04,0.924529E+04,0.929310E+04, &
0.934122E+04,0.938959E+04,0.943804E+04,0.948687E+04,0.953575E+04, &
0.958494E+04,0.963442E+04,0.968395E+04,0.973384E+04,0.978383E+04, &
0.983412E+04,0.988468E+04,0.993534E+04,0.998630E+04,0.100374E+05, &
0.100888E+05,0.101406E+05,0.101923E+05,0.102444E+05,0.102966E+05, &
0.103492E+05,0.104020E+05,0.104550E+05,0.105082E+05,0.105616E+05, &
0.106153E+05,0.106693E+05,0.107234E+05,0.107779E+05,0.108325E+05, &
0.108874E+05,0.109425E+05,0.109978E+05,0.110535E+05,0.111092E+05, &
0.111653E+05,0.112217E+05,0.112782E+05,0.113350E+05,0.113920E+05, &
0.114493E+05,0.115070E+05,0.115646E+05,0.116228E+05,0.116809E+05, &
0.117396E+05,0.117984E+05,0.118574E+05,0.119167E+05,0.119762E+05, &
0.120360E+05,0.120962E+05,0.121564E+05,0.122170E+05,0.122778E+05, &
0.123389E+05,0.124004E+05,0.124619E+05,0.125238E+05,0.125859E+05, &
0.126484E+05,0.127111E+05,0.127739E+05,0.128372E+05,0.129006E+05, &
0.129644E+05,0.130285E+05,0.130927E+05,0.131573E+05,0.132220E+05, &
0.132872E+05,0.133526E+05,0.134182E+05,0.134842E+05,0.135503E+05, &
0.136168E+05,0.136836E+05,0.137505E+05,0.138180E+05,0.138854E+05/
      data (ES(IES),IES=1426,1520) / &
0.139534E+05,0.140216E+05,0.140900E+05,0.141588E+05,0.142277E+05, &
0.142971E+05,0.143668E+05,0.144366E+05,0.145069E+05,0.145773E+05, &
0.146481E+05,0.147192E+05,0.147905E+05,0.148622E+05,0.149341E+05, &
0.150064E+05,0.150790E+05,0.151517E+05,0.152250E+05,0.152983E+05, &
0.153721E+05,0.154462E+05,0.155205E+05,0.155952E+05,0.156701E+05, &
0.157454E+05,0.158211E+05,0.158969E+05,0.159732E+05,0.160496E+05, &
0.161265E+05,0.162037E+05,0.162811E+05,0.163589E+05,0.164369E+05, &
0.165154E+05,0.165942E+05,0.166732E+05,0.167526E+05,0.168322E+05, &
0.169123E+05,0.169927E+05,0.170733E+05,0.171543E+05,0.172356E+05, &
0.173173E+05,0.173993E+05,0.174815E+05,0.175643E+05,0.176471E+05, &
0.177305E+05,0.178143E+05,0.178981E+05,0.179826E+05,0.180671E+05, &
0.181522E+05,0.182377E+05,0.183232E+05,0.184093E+05,0.184955E+05, &
0.185823E+05,0.186695E+05,0.187568E+05,0.188447E+05,0.189326E+05, &
0.190212E+05,0.191101E+05,0.191991E+05,0.192887E+05,0.193785E+05, &
0.194688E+05,0.195595E+05,0.196503E+05,0.197417E+05,0.198332E+05, &
0.199253E+05,0.200178E+05,0.201105E+05,0.202036E+05,0.202971E+05, &
0.203910E+05,0.204853E+05,0.205798E+05,0.206749E+05,0.207701E+05, &
0.208659E+05,0.209621E+05,0.210584E+05,0.211554E+05,0.212524E+05, &
0.213501E+05,0.214482E+05,0.215465E+05,0.216452E+05,0.217442E+05/
      data (ES(IES),IES=1521,1552) / &
0.218439E+05,0.219439E+05,0.220440E+05,0.221449E+05,0.222457E+05, &
0.223473E+05,0.224494E+05,0.225514E+05,0.226542E+05,0.227571E+05, &
0.228606E+05,0.229646E+05,0.230687E+05,0.231734E+05,0.232783E+05, &
0.233839E+05,0.234898E+05,0.235960E+05,0.237027E+05,0.238097E+05, &
0.239173E+05,0.240254E+05,0.241335E+05,0.242424E+05,0.243514E+05, &
0.244611E+05,0.245712E+05,0.246814E+05,0.247923E+05,0.249034E+05, &
0.250152E+05,0.250152E+05/

do I = 1, npnts

  ! Compute the factor that converts from sat vapour pressure in a
  ! pure water system to sat vapour pressure in air, FSUBW.  This formula
  ! is taken from equation A4.7 of Adrian Gill's book: Atmosphere-Ocean
  ! Dynamics.  Note that his formula works in terms of pressure in MB and
  ! temperature in Celsius, so conversion of units leads to the slightly
  ! different equation used here.

  FSUBW = one + 1.0E-8_kind_real * P(I) * &
                          (4.5_kind_real+ 6.0E-4_kind_real * (T(I) - t0c) * (T(I) - t0c))

  ! use the lookup table to find saturated vapour pressure, and stored it in QS.

  TT = max (T_low, T(I))
  TT = min (T_high, TT)

  ATABLE = (TT - T_low + delta_T) / delta_T
  ITABLE = ATABLE
  ATABLE = ATABLE - ITABLE

  QS(I) = (one - ATABLE) * ES(ITABLE) + ATABLE * ES(ITABLE + 1)

  ! Multiply by FSUBW to convert to saturated vapour pressure in air (equation A4.6
  ! of Adrian Gill's book)

  QS(I) = QS(I) * FSUBW

  ! Now form the accurate expression for QS, which is a rearranged
  ! version of equation A4.3 of Gill's book.

  ! Note that at very low pressures we apply a fix, to prevent a
  ! singularity (Qsat tends to 1.0 kg/kg).

  QS(I) = (epsilon * QS(I)) / (max (P(I), QS(I)) - one_minus_epsilon * QS(I))

end do

end subroutine ufo_rttovonedvarcheck_Qsat

!-------------------------------------------------------------------------------
!> Saturation Specific Humidity Scheme: Vapour to Liquid.
!!
!! \details Heritage: Ops_QsatWat.inc
!!
!! Returns a saturation mixing ratio given a temperature and pressure
!! using saturation vapour pressure calculated using the Goff-Gratch
!! Formulaie, adopted by the WMO as taken from Landolt-Bornstein, 1987
!! Numerical Data and Functional Relationships in Science and
!! Technology.  Group V/Vol 4B Meteorology.  Physical and Chemical
!! Properties of Air, P35
!!
!! Values in the lookup table are over water above and below 0 degrees C.
!!
!! Note: For vapour pressure oever water this formula is valid for
!! temperatures between 373K and 223K.  The values for saturated vapour
!! over water in the lookup table below are out of the lower end of
!! this range.  However it is standard WMO practice to use the formula
!! below its accepted range for use with the calculation of dew points
!! in the upper atmosphere
!!
!! Method: <br>
!!   uses lookup tables to find eSAT, calculates qSAT directly from that.
!!
!! \author Met Office
!!
!! \date 09/06/2020: Created
!!
subroutine ufo_rttovonedvarcheck_QsatWat (QS, &
                                          T,  &
                                          P,  &
                                          npnts)

implicit none

! subroutine arguments:
integer                     :: npnts     ! Points (=horizontal dimensions) being processed by qSAT scheme.
real(kind=kind_real)        :: T(npnts)  ! Temperature (K).
real(kind=kind_real)        :: P(npnts)  ! Pressure (Pa).
real(kind=kind_real)        :: QS(npnts) ! Saturation mixing ratio at temperature T and pressure P (KG/KG)

! Local declarations:
real(kind=kind_real), parameter :: one_minus_epsilon = 1.0 - epsilon
real(kind=kind_real), parameter :: T_low = 183.15
real(kind=kind_real), parameter :: T_high = 338.15
real(kind=kind_real), parameter :: delta_T = 0.1
integer, parameter          :: N = ((T_high - T_low + (delta_T * 0.5)) / delta_T) + 1.0  ! gives N = 1551
integer                     :: ITABLE
real                        :: ATABLE
real                        :: FSUBW
real                        :: TT
integer                     :: I
integer                     :: IES
real                        :: ES(0:N + 1)   ! Table of saturation water vapour pressure (PA)
character(len=*), parameter :: RoutineName = "ufo_rttovonedvarcheck_QsatWat"

! Note: 0 element is a repeat of 1st element to cater for special case
!       of low temperatures (.LE.T_low) for which the array index is
!       rounded down due to machine precision.

data (ES(IES),IES=    0, 95) / 0.186905E-01, &
0.186905E-01,0.190449E-01,0.194059E-01,0.197727E-01,0.201462E-01, &
0.205261E-01,0.209122E-01,0.213052E-01,0.217050E-01,0.221116E-01, &
0.225252E-01,0.229463E-01,0.233740E-01,0.238090E-01,0.242518E-01, &
0.247017E-01,0.251595E-01,0.256252E-01,0.260981E-01,0.265795E-01, &
0.270691E-01,0.275667E-01,0.280733E-01,0.285876E-01,0.291105E-01, &
0.296429E-01,0.301835E-01,0.307336E-01,0.312927E-01,0.318611E-01, &
0.324390E-01,0.330262E-01,0.336232E-01,0.342306E-01,0.348472E-01, &
0.354748E-01,0.361117E-01,0.367599E-01,0.374185E-01,0.380879E-01, &
0.387689E-01,0.394602E-01,0.401626E-01,0.408771E-01,0.416033E-01, &
0.423411E-01,0.430908E-01,0.438524E-01,0.446263E-01,0.454124E-01, &
0.462122E-01,0.470247E-01,0.478491E-01,0.486874E-01,0.495393E-01, &
0.504057E-01,0.512847E-01,0.521784E-01,0.530853E-01,0.540076E-01, &
0.549444E-01,0.558959E-01,0.568633E-01,0.578448E-01,0.588428E-01, &
0.598566E-01,0.608858E-01,0.619313E-01,0.629926E-01,0.640706E-01, &
0.651665E-01,0.662795E-01,0.674095E-01,0.685570E-01,0.697219E-01, &
0.709063E-01,0.721076E-01,0.733284E-01,0.745679E-01,0.758265E-01, &
0.771039E-01,0.784026E-01,0.797212E-01,0.810577E-01,0.824164E-01, &
0.837971E-01,0.851970E-01,0.866198E-01,0.880620E-01,0.895281E-01, &
0.910178E-01,0.925278E-01,0.940622E-01,0.956177E-01,0.971984E-01/
data (ES(IES),IES= 96,190) / &
0.988051E-01,0.100433E+00,0.102085E+00,0.103764E+00,0.105467E+00, &
0.107196E+00,0.108953E+00,0.110732E+00,0.112541E+00,0.114376E+00, &
0.116238E+00,0.118130E+00,0.120046E+00,0.121993E+00,0.123969E+00, &
0.125973E+00,0.128009E+00,0.130075E+00,0.132167E+00,0.134296E+00, &
0.136452E+00,0.138642E+00,0.140861E+00,0.143115E+00,0.145404E+00, &
0.147723E+00,0.150078E+00,0.152466E+00,0.154889E+00,0.157346E+00, &
0.159841E+00,0.162372E+00,0.164939E+00,0.167545E+00,0.170185E+00, &
0.172866E+00,0.175584E+00,0.178340E+00,0.181139E+00,0.183977E+00, &
0.186855E+00,0.189773E+00,0.192737E+00,0.195736E+00,0.198783E+00, &
0.201875E+00,0.205007E+00,0.208186E+00,0.211409E+00,0.214676E+00, &
0.217993E+00,0.221355E+00,0.224764E+00,0.228220E+00,0.231728E+00, &
0.235284E+00,0.238888E+00,0.242542E+00,0.246251E+00,0.250010E+00, &
0.253821E+00,0.257688E+00,0.261602E+00,0.265575E+00,0.269607E+00, &
0.273689E+00,0.277830E+00,0.282027E+00,0.286287E+00,0.290598E+00, &
0.294972E+00,0.299405E+00,0.303904E+00,0.308462E+00,0.313082E+00, &
0.317763E+00,0.322512E+00,0.327324E+00,0.332201E+00,0.337141E+00, &
0.342154E+00,0.347234E+00,0.352387E+00,0.357601E+00,0.362889E+00, &
0.368257E+00,0.373685E+00,0.379194E+00,0.384773E+00,0.390433E+00, &
0.396159E+00,0.401968E+00,0.407861E+00,0.413820E+00,0.419866E+00/
data (ES(IES),IES=191,285) / &
0.425999E+00,0.432203E+00,0.438494E+00,0.444867E+00,0.451332E+00, &
0.457879E+00,0.464510E+00,0.471226E+00,0.478037E+00,0.484935E+00, &
0.491920E+00,0.499005E+00,0.506181E+00,0.513447E+00,0.520816E+00, &
0.528279E+00,0.535835E+00,0.543497E+00,0.551256E+00,0.559113E+00, &
0.567081E+00,0.575147E+00,0.583315E+00,0.591585E+00,0.599970E+00, &
0.608472E+00,0.617069E+00,0.625785E+00,0.634609E+00,0.643556E+00, &
0.652611E+00,0.661782E+00,0.671077E+00,0.680487E+00,0.690015E+00, &
0.699656E+00,0.709433E+00,0.719344E+00,0.729363E+00,0.739518E+00, &
0.749795E+00,0.760217E+00,0.770763E+00,0.781454E+00,0.792258E+00, &
0.803208E+00,0.814309E+00,0.825528E+00,0.836914E+00,0.848422E+00, &
0.860086E+00,0.871891E+00,0.883837E+00,0.895944E+00,0.908214E+00, &
0.920611E+00,0.933175E+00,0.945890E+00,0.958776E+00,0.971812E+00, &
0.985027E+00,0.998379E+00,0.101193E+01,0.102561E+01,0.103949E+01, &
0.105352E+01,0.106774E+01,0.108213E+01,0.109669E+01,0.111144E+01, &
0.112636E+01,0.114148E+01,0.115676E+01,0.117226E+01,0.118791E+01, &
0.120377E+01,0.121984E+01,0.123608E+01,0.125252E+01,0.126919E+01, &
0.128604E+01,0.130309E+01,0.132036E+01,0.133782E+01,0.135549E+01, &
0.137339E+01,0.139150E+01,0.140984E+01,0.142839E+01,0.144715E+01, &
0.146616E+01,0.148538E+01,0.150482E+01,0.152450E+01,0.154445E+01/
data (ES(IES),IES=286,380) / &
0.156459E+01,0.158502E+01,0.160564E+01,0.162654E+01,0.164766E+01, &
0.166906E+01,0.169070E+01,0.171257E+01,0.173473E+01,0.175718E+01, &
0.177984E+01,0.180282E+01,0.182602E+01,0.184951E+01,0.187327E+01, &
0.189733E+01,0.192165E+01,0.194629E+01,0.197118E+01,0.199636E+01, &
0.202185E+01,0.204762E+01,0.207372E+01,0.210010E+01,0.212678E+01, &
0.215379E+01,0.218109E+01,0.220873E+01,0.223668E+01,0.226497E+01, &
0.229357E+01,0.232249E+01,0.235176E+01,0.238134E+01,0.241129E+01, &
0.244157E+01,0.247217E+01,0.250316E+01,0.253447E+01,0.256617E+01, &
0.259821E+01,0.263064E+01,0.266341E+01,0.269661E+01,0.273009E+01, &
0.276403E+01,0.279834E+01,0.283302E+01,0.286811E+01,0.290358E+01, &
0.293943E+01,0.297571E+01,0.301236E+01,0.304946E+01,0.308702E+01, &
0.312491E+01,0.316326E+01,0.320208E+01,0.324130E+01,0.328092E+01, &
0.332102E+01,0.336162E+01,0.340264E+01,0.344407E+01,0.348601E+01, &
0.352838E+01,0.357118E+01,0.361449E+01,0.365834E+01,0.370264E+01, &
0.374737E+01,0.379265E+01,0.383839E+01,0.388469E+01,0.393144E+01, &
0.397876E+01,0.402656E+01,0.407492E+01,0.412378E+01,0.417313E+01, &
0.422306E+01,0.427359E+01,0.432454E+01,0.437617E+01,0.442834E+01, &
0.448102E+01,0.453433E+01,0.458816E+01,0.464253E+01,0.469764E+01, &
0.475321E+01,0.480942E+01,0.486629E+01,0.492372E+01,0.498173E+01/
data (ES(IES),IES=381,475) / &
0.504041E+01,0.509967E+01,0.515962E+01,0.522029E+01,0.528142E+01, &
0.534337E+01,0.540595E+01,0.546912E+01,0.553292E+01,0.559757E+01, &
0.566273E+01,0.572864E+01,0.579532E+01,0.586266E+01,0.593075E+01, &
0.599952E+01,0.606895E+01,0.613918E+01,0.621021E+01,0.628191E+01, &
0.635433E+01,0.642755E+01,0.650162E+01,0.657639E+01,0.665188E+01, &
0.672823E+01,0.680532E+01,0.688329E+01,0.696198E+01,0.704157E+01, &
0.712206E+01,0.720319E+01,0.728534E+01,0.736829E+01,0.745204E+01, &
0.753671E+01,0.762218E+01,0.770860E+01,0.779588E+01,0.788408E+01, &
0.797314E+01,0.806318E+01,0.815408E+01,0.824599E+01,0.833874E+01, &
0.843254E+01,0.852721E+01,0.862293E+01,0.871954E+01,0.881724E+01, &
0.891579E+01,0.901547E+01,0.911624E+01,0.921778E+01,0.932061E+01, &
0.942438E+01,0.952910E+01,0.963497E+01,0.974181E+01,0.984982E+01, &
0.995887E+01,0.100690E+02,0.101804E+02,0.102926E+02,0.104063E+02, &
0.105210E+02,0.106367E+02,0.107536E+02,0.108719E+02,0.109912E+02, &
0.111116E+02,0.112333E+02,0.113563E+02,0.114804E+02,0.116056E+02, &
0.117325E+02,0.118602E+02,0.119892E+02,0.121197E+02,0.122513E+02, &
0.123844E+02,0.125186E+02,0.126543E+02,0.127912E+02,0.129295E+02, &
0.130691E+02,0.132101E+02,0.133527E+02,0.134965E+02,0.136415E+02, &
0.137882E+02,0.139361E+02,0.140855E+02,0.142366E+02,0.143889E+02/
data (ES(IES),IES=476,570) / &
0.145429E+02,0.146982E+02,0.148552E+02,0.150135E+02,0.151735E+02, &
0.153349E+02,0.154979E+02,0.156624E+02,0.158286E+02,0.159965E+02, &
0.161659E+02,0.163367E+02,0.165094E+02,0.166838E+02,0.168597E+02, &
0.170375E+02,0.172168E+02,0.173979E+02,0.175806E+02,0.177651E+02, &
0.179513E+02,0.181394E+02,0.183293E+02,0.185210E+02,0.187146E+02, &
0.189098E+02,0.191066E+02,0.193059E+02,0.195065E+02,0.197095E+02, &
0.199142E+02,0.201206E+02,0.203291E+02,0.205397E+02,0.207522E+02, &
0.209664E+02,0.211831E+02,0.214013E+02,0.216221E+02,0.218448E+02, &
0.220692E+02,0.222959E+02,0.225250E+02,0.227559E+02,0.229887E+02, &
0.232239E+02,0.234614E+02,0.237014E+02,0.239428E+02,0.241872E+02, &
0.244335E+02,0.246824E+02,0.249332E+02,0.251860E+02,0.254419E+02, &
0.256993E+02,0.259600E+02,0.262225E+02,0.264873E+02,0.267552E+02, &
0.270248E+02,0.272970E+02,0.275719E+02,0.278497E+02,0.281295E+02, &
0.284117E+02,0.286965E+02,0.289843E+02,0.292743E+02,0.295671E+02, &
0.298624E+02,0.301605E+02,0.304616E+02,0.307650E+02,0.310708E+02, &
0.313803E+02,0.316915E+02,0.320064E+02,0.323238E+02,0.326437E+02, &
0.329666E+02,0.332928E+02,0.336215E+02,0.339534E+02,0.342885E+02, &
0.346263E+02,0.349666E+02,0.353109E+02,0.356572E+02,0.360076E+02, &
0.363606E+02,0.367164E+02,0.370757E+02,0.374383E+02,0.378038E+02/
data (ES(IES),IES=571,665) / &
0.381727E+02,0.385453E+02,0.389206E+02,0.392989E+02,0.396807E+02, &
0.400663E+02,0.404555E+02,0.408478E+02,0.412428E+02,0.416417E+02, &
0.420445E+02,0.424502E+02,0.428600E+02,0.432733E+02,0.436900E+02, &
0.441106E+02,0.445343E+02,0.449620E+02,0.453930E+02,0.458280E+02, &
0.462672E+02,0.467096E+02,0.471561E+02,0.476070E+02,0.480610E+02, &
0.485186E+02,0.489813E+02,0.494474E+02,0.499170E+02,0.503909E+02, &
0.508693E+02,0.513511E+02,0.518376E+02,0.523277E+02,0.528232E+02, &
0.533213E+02,0.538240E+02,0.543315E+02,0.548437E+02,0.553596E+02, &
0.558802E+02,0.564046E+02,0.569340E+02,0.574672E+02,0.580061E+02, &
0.585481E+02,0.590963E+02,0.596482E+02,0.602041E+02,0.607649E+02, &
0.613311E+02,0.619025E+02,0.624779E+02,0.630574E+02,0.636422E+02, &
0.642324E+02,0.648280E+02,0.654278E+02,0.660332E+02,0.666426E+02, &
0.672577E+02,0.678771E+02,0.685034E+02,0.691328E+02,0.697694E+02, &
0.704103E+02,0.710556E+02,0.717081E+02,0.723639E+02,0.730269E+02, &
0.736945E+02,0.743681E+02,0.750463E+02,0.757309E+02,0.764214E+02, &
0.771167E+02,0.778182E+02,0.785246E+02,0.792373E+02,0.799564E+02, &
0.806804E+02,0.814109E+02,0.821479E+02,0.828898E+02,0.836384E+02, &
0.843922E+02,0.851525E+02,0.859198E+02,0.866920E+02,0.874712E+02, &
0.882574E+02,0.890486E+02,0.898470E+02,0.906525E+02,0.914634E+02/
data (ES(IES),IES=666,760) / &
0.922814E+02,0.931048E+02,0.939356E+02,0.947736E+02,0.956171E+02, &
0.964681E+02,0.973246E+02,0.981907E+02,0.990605E+02,0.999399E+02, &
0.100825E+03,0.101718E+03,0.102617E+03,0.103523E+03,0.104438E+03, &
0.105358E+03,0.106287E+03,0.107221E+03,0.108166E+03,0.109115E+03, &
0.110074E+03,0.111039E+03,0.112012E+03,0.112992E+03,0.113981E+03, &
0.114978E+03,0.115981E+03,0.116993E+03,0.118013E+03,0.119041E+03, &
0.120077E+03,0.121122E+03,0.122173E+03,0.123234E+03,0.124301E+03, &
0.125377E+03,0.126463E+03,0.127556E+03,0.128657E+03,0.129769E+03, &
0.130889E+03,0.132017E+03,0.133152E+03,0.134299E+03,0.135453E+03, &
0.136614E+03,0.137786E+03,0.138967E+03,0.140158E+03,0.141356E+03, &
0.142565E+03,0.143781E+03,0.145010E+03,0.146247E+03,0.147491E+03, &
0.148746E+03,0.150011E+03,0.151284E+03,0.152571E+03,0.153862E+03, &
0.155168E+03,0.156481E+03,0.157805E+03,0.159137E+03,0.160478E+03, &
0.161832E+03,0.163198E+03,0.164569E+03,0.165958E+03,0.167348E+03, &
0.168757E+03,0.170174E+03,0.171599E+03,0.173037E+03,0.174483E+03, &
0.175944E+03,0.177414E+03,0.178892E+03,0.180387E+03,0.181886E+03, &
0.183402E+03,0.184930E+03,0.186463E+03,0.188012E+03,0.189571E+03, &
0.191146E+03,0.192730E+03,0.194320E+03,0.195930E+03,0.197546E+03, &
0.199175E+03,0.200821E+03,0.202473E+03,0.204142E+03,0.205817E+03/
data (ES(IES),IES=761,855) / &
0.207510E+03,0.209216E+03,0.210928E+03,0.212658E+03,0.214398E+03, &
0.216152E+03,0.217920E+03,0.219698E+03,0.221495E+03,0.223297E+03, &
0.225119E+03,0.226951E+03,0.228793E+03,0.230654E+03,0.232522E+03, &
0.234413E+03,0.236311E+03,0.238223E+03,0.240151E+03,0.242090E+03, &
0.244049E+03,0.246019E+03,0.248000E+03,0.249996E+03,0.252009E+03, &
0.254037E+03,0.256077E+03,0.258128E+03,0.260200E+03,0.262284E+03, &
0.264384E+03,0.266500E+03,0.268629E+03,0.270779E+03,0.272936E+03, &
0.275110E+03,0.277306E+03,0.279509E+03,0.281734E+03,0.283966E+03, &
0.286227E+03,0.288494E+03,0.290780E+03,0.293083E+03,0.295398E+03, &
0.297737E+03,0.300089E+03,0.302453E+03,0.304841E+03,0.307237E+03, &
0.309656E+03,0.312095E+03,0.314541E+03,0.317012E+03,0.319496E+03, &
0.322005E+03,0.324527E+03,0.327063E+03,0.329618E+03,0.332193E+03, &
0.334788E+03,0.337396E+03,0.340025E+03,0.342673E+03,0.345329E+03, &
0.348019E+03,0.350722E+03,0.353440E+03,0.356178E+03,0.358938E+03, &
0.361718E+03,0.364513E+03,0.367322E+03,0.370160E+03,0.373012E+03, &
0.375885E+03,0.378788E+03,0.381691E+03,0.384631E+03,0.387579E+03, &
0.390556E+03,0.393556E+03,0.396563E+03,0.399601E+03,0.402646E+03, &
0.405730E+03,0.408829E+03,0.411944E+03,0.415083E+03,0.418236E+03, &
0.421422E+03,0.424632E+03,0.427849E+03,0.431099E+03,0.434365E+03/
data (ES(IES),IES=856,950) / &
0.437655E+03,0.440970E+03,0.444301E+03,0.447666E+03,0.451038E+03, &
0.454445E+03,0.457876E+03,0.461316E+03,0.464790E+03,0.468281E+03, &
0.471798E+03,0.475342E+03,0.478902E+03,0.482497E+03,0.486101E+03, &
0.489741E+03,0.493408E+03,0.497083E+03,0.500804E+03,0.504524E+03, &
0.508290E+03,0.512074E+03,0.515877E+03,0.519717E+03,0.523566E+03, &
0.527462E+03,0.531367E+03,0.535301E+03,0.539264E+03,0.543245E+03, &
0.547265E+03,0.551305E+03,0.555363E+03,0.559462E+03,0.563579E+03, &
0.567727E+03,0.571905E+03,0.576102E+03,0.580329E+03,0.584576E+03, &
0.588865E+03,0.593185E+03,0.597514E+03,0.601885E+03,0.606276E+03, &
0.610699E+03,0.615151E+03,0.619625E+03,0.624140E+03,0.628671E+03, &
0.633243E+03,0.637845E+03,0.642465E+03,0.647126E+03,0.651806E+03, &
0.656527E+03,0.661279E+03,0.666049E+03,0.670861E+03,0.675692E+03, &
0.680566E+03,0.685471E+03,0.690396E+03,0.695363E+03,0.700350E+03, &
0.705381E+03,0.710444E+03,0.715527E+03,0.720654E+03,0.725801E+03, &
0.730994E+03,0.736219E+03,0.741465E+03,0.746756E+03,0.752068E+03, &
0.757426E+03,0.762819E+03,0.768231E+03,0.773692E+03,0.779172E+03, &
0.784701E+03,0.790265E+03,0.795849E+03,0.801483E+03,0.807137E+03, &
0.812842E+03,0.818582E+03,0.824343E+03,0.830153E+03,0.835987E+03, &
0.841871E+03,0.847791E+03,0.853733E+03,0.859727E+03,0.865743E+03/
data (ES(IES),IES=951,1045) / &
0.871812E+03,0.877918E+03,0.884046E+03,0.890228E+03,0.896433E+03, &
0.902690E+03,0.908987E+03,0.915307E+03,0.921681E+03,0.928078E+03, &
0.934531E+03,0.941023E+03,0.947539E+03,0.954112E+03,0.960708E+03, &
0.967361E+03,0.974053E+03,0.980771E+03,0.987545E+03,0.994345E+03, &
0.100120E+04,0.100810E+04,0.101502E+04,0.102201E+04,0.102902E+04, &
0.103608E+04,0.104320E+04,0.105033E+04,0.105753E+04,0.106475E+04, &
0.107204E+04,0.107936E+04,0.108672E+04,0.109414E+04,0.110158E+04, &
0.110908E+04,0.111663E+04,0.112421E+04,0.113185E+04,0.113952E+04, &
0.114725E+04,0.115503E+04,0.116284E+04,0.117071E+04,0.117861E+04, &
0.118658E+04,0.119459E+04,0.120264E+04,0.121074E+04,0.121888E+04, &
0.122709E+04,0.123534E+04,0.124362E+04,0.125198E+04,0.126036E+04, &
0.126881E+04,0.127731E+04,0.128584E+04,0.129444E+04,0.130307E+04, &
0.131177E+04,0.132053E+04,0.132931E+04,0.133817E+04,0.134705E+04, &
0.135602E+04,0.136503E+04,0.137407E+04,0.138319E+04,0.139234E+04, &
0.140156E+04,0.141084E+04,0.142015E+04,0.142954E+04,0.143896E+04, &
0.144845E+04,0.145800E+04,0.146759E+04,0.147725E+04,0.148694E+04, &
0.149672E+04,0.150655E+04,0.151641E+04,0.152635E+04,0.153633E+04, &
0.154639E+04,0.155650E+04,0.156665E+04,0.157688E+04,0.158715E+04, &
0.159750E+04,0.160791E+04,0.161836E+04,0.162888E+04,0.163945E+04/
data (ES(IES),IES=1046,1140) / &
0.165010E+04,0.166081E+04,0.167155E+04,0.168238E+04,0.169325E+04, &
0.170420E+04,0.171522E+04,0.172627E+04,0.173741E+04,0.174859E+04, &
0.175986E+04,0.177119E+04,0.178256E+04,0.179402E+04,0.180552E+04, &
0.181711E+04,0.182877E+04,0.184046E+04,0.185224E+04,0.186407E+04, &
0.187599E+04,0.188797E+04,0.190000E+04,0.191212E+04,0.192428E+04, &
0.193653E+04,0.194886E+04,0.196122E+04,0.197368E+04,0.198618E+04, &
0.199878E+04,0.201145E+04,0.202416E+04,0.203698E+04,0.204983E+04, &
0.206278E+04,0.207580E+04,0.208887E+04,0.210204E+04,0.211525E+04, &
0.212856E+04,0.214195E+04,0.215538E+04,0.216892E+04,0.218249E+04, &
0.219618E+04,0.220994E+04,0.222375E+04,0.223766E+04,0.225161E+04, &
0.226567E+04,0.227981E+04,0.229399E+04,0.230829E+04,0.232263E+04, &
0.233708E+04,0.235161E+04,0.236618E+04,0.238087E+04,0.239560E+04, &
0.241044E+04,0.242538E+04,0.244035E+04,0.245544E+04,0.247057E+04, &
0.248583E+04,0.250116E+04,0.251654E+04,0.253204E+04,0.254759E+04, &
0.256325E+04,0.257901E+04,0.259480E+04,0.261073E+04,0.262670E+04, &
0.264279E+04,0.265896E+04,0.267519E+04,0.269154E+04,0.270794E+04, &
0.272447E+04,0.274108E+04,0.275774E+04,0.277453E+04,0.279137E+04, &
0.280834E+04,0.282540E+04,0.284251E+04,0.285975E+04,0.287704E+04, &
0.289446E+04,0.291198E+04,0.292954E+04,0.294725E+04,0.296499E+04/
data (ES(IES),IES=1141,1235) / &
0.298288E+04,0.300087E+04,0.301890E+04,0.303707E+04,0.305529E+04, &
0.307365E+04,0.309211E+04,0.311062E+04,0.312927E+04,0.314798E+04, &
0.316682E+04,0.318577E+04,0.320477E+04,0.322391E+04,0.324310E+04, &
0.326245E+04,0.328189E+04,0.330138E+04,0.332103E+04,0.334073E+04, &
0.336058E+04,0.338053E+04,0.340054E+04,0.342069E+04,0.344090E+04, &
0.346127E+04,0.348174E+04,0.350227E+04,0.352295E+04,0.354369E+04, &
0.356458E+04,0.358559E+04,0.360664E+04,0.362787E+04,0.364914E+04, &
0.367058E+04,0.369212E+04,0.371373E+04,0.373548E+04,0.375731E+04, &
0.377929E+04,0.380139E+04,0.382355E+04,0.384588E+04,0.386826E+04, &
0.389081E+04,0.391348E+04,0.393620E+04,0.395910E+04,0.398205E+04, &
0.400518E+04,0.402843E+04,0.405173E+04,0.407520E+04,0.409875E+04, &
0.412246E+04,0.414630E+04,0.417019E+04,0.419427E+04,0.421840E+04, &
0.424272E+04,0.426715E+04,0.429165E+04,0.431634E+04,0.434108E+04, &
0.436602E+04,0.439107E+04,0.441618E+04,0.444149E+04,0.446685E+04, &
0.449241E+04,0.451810E+04,0.454385E+04,0.456977E+04,0.459578E+04, &
0.462197E+04,0.464830E+04,0.467468E+04,0.470127E+04,0.472792E+04, &
0.475477E+04,0.478175E+04,0.480880E+04,0.483605E+04,0.486336E+04, &
0.489087E+04,0.491853E+04,0.494623E+04,0.497415E+04,0.500215E+04, &
0.503034E+04,0.505867E+04,0.508707E+04,0.511568E+04,0.514436E+04/
data (ES(IES),IES=1236,1330) / &
0.517325E+04,0.520227E+04,0.523137E+04,0.526068E+04,0.529005E+04, &
0.531965E+04,0.534939E+04,0.537921E+04,0.540923E+04,0.543932E+04, &
0.546965E+04,0.550011E+04,0.553064E+04,0.556139E+04,0.559223E+04, &
0.562329E+04,0.565449E+04,0.568577E+04,0.571727E+04,0.574884E+04, &
0.578064E+04,0.581261E+04,0.584464E+04,0.587692E+04,0.590924E+04, &
0.594182E+04,0.597455E+04,0.600736E+04,0.604039E+04,0.607350E+04, &
0.610685E+04,0.614036E+04,0.617394E+04,0.620777E+04,0.624169E+04, &
0.627584E+04,0.631014E+04,0.634454E+04,0.637918E+04,0.641390E+04, &
0.644887E+04,0.648400E+04,0.651919E+04,0.655467E+04,0.659021E+04, &
0.662599E+04,0.666197E+04,0.669800E+04,0.673429E+04,0.677069E+04, &
0.680735E+04,0.684415E+04,0.688104E+04,0.691819E+04,0.695543E+04, &
0.699292E+04,0.703061E+04,0.706837E+04,0.710639E+04,0.714451E+04, &
0.718289E+04,0.722143E+04,0.726009E+04,0.729903E+04,0.733802E+04, &
0.737729E+04,0.741676E+04,0.745631E+04,0.749612E+04,0.753602E+04, &
0.757622E+04,0.761659E+04,0.765705E+04,0.769780E+04,0.773863E+04, &
0.777975E+04,0.782106E+04,0.786246E+04,0.790412E+04,0.794593E+04, &
0.798802E+04,0.803028E+04,0.807259E+04,0.811525E+04,0.815798E+04, &
0.820102E+04,0.824427E+04,0.828757E+04,0.833120E+04,0.837493E+04, &
0.841895E+04,0.846313E+04,0.850744E+04,0.855208E+04,0.859678E+04/
data (ES(IES),IES=1331,1425) / &
0.864179E+04,0.868705E+04,0.873237E+04,0.877800E+04,0.882374E+04, &
0.886979E+04,0.891603E+04,0.896237E+04,0.900904E+04,0.905579E+04, &
0.910288E+04,0.915018E+04,0.919758E+04,0.924529E+04,0.929310E+04, &
0.934122E+04,0.938959E+04,0.943804E+04,0.948687E+04,0.953575E+04, &
0.958494E+04,0.963442E+04,0.968395E+04,0.973384E+04,0.978383E+04, &
0.983412E+04,0.988468E+04,0.993534E+04,0.998630E+04,0.100374E+05, &
0.100888E+05,0.101406E+05,0.101923E+05,0.102444E+05,0.102966E+05, &
0.103492E+05,0.104020E+05,0.104550E+05,0.105082E+05,0.105616E+05, &
0.106153E+05,0.106693E+05,0.107234E+05,0.107779E+05,0.108325E+05, &
0.108874E+05,0.109425E+05,0.109978E+05,0.110535E+05,0.111092E+05, &
0.111653E+05,0.112217E+05,0.112782E+05,0.113350E+05,0.113920E+05, &
0.114493E+05,0.115070E+05,0.115646E+05,0.116228E+05,0.116809E+05, &
0.117396E+05,0.117984E+05,0.118574E+05,0.119167E+05,0.119762E+05, &
0.120360E+05,0.120962E+05,0.121564E+05,0.122170E+05,0.122778E+05, &
0.123389E+05,0.124004E+05,0.124619E+05,0.125238E+05,0.125859E+05, &
0.126484E+05,0.127111E+05,0.127739E+05,0.128372E+05,0.129006E+05, &
0.129644E+05,0.130285E+05,0.130927E+05,0.131573E+05,0.132220E+05, &
0.132872E+05,0.133526E+05,0.134182E+05,0.134842E+05,0.135503E+05, &
0.136168E+05,0.136836E+05,0.137505E+05,0.138180E+05,0.138854E+05/
data (ES(IES),IES=1426,1520) / &
0.139534E+05,0.140216E+05,0.140900E+05,0.141588E+05,0.142277E+05, &
0.142971E+05,0.143668E+05,0.144366E+05,0.145069E+05,0.145773E+05, &
0.146481E+05,0.147192E+05,0.147905E+05,0.148622E+05,0.149341E+05, &
0.150064E+05,0.150790E+05,0.151517E+05,0.152250E+05,0.152983E+05, &
0.153721E+05,0.154462E+05,0.155205E+05,0.155952E+05,0.156701E+05, &
0.157454E+05,0.158211E+05,0.158969E+05,0.159732E+05,0.160496E+05, &
0.161265E+05,0.162037E+05,0.162811E+05,0.163589E+05,0.164369E+05, &
0.165154E+05,0.165942E+05,0.166732E+05,0.167526E+05,0.168322E+05, &
0.169123E+05,0.169927E+05,0.170733E+05,0.171543E+05,0.172356E+05, &
0.173173E+05,0.173993E+05,0.174815E+05,0.175643E+05,0.176471E+05, &
0.177305E+05,0.178143E+05,0.178981E+05,0.179826E+05,0.180671E+05, &
0.181522E+05,0.182377E+05,0.183232E+05,0.184093E+05,0.184955E+05, &
0.185823E+05,0.186695E+05,0.187568E+05,0.188447E+05,0.189326E+05, &
0.190212E+05,0.191101E+05,0.191991E+05,0.192887E+05,0.193785E+05, &
0.194688E+05,0.195595E+05,0.196503E+05,0.197417E+05,0.198332E+05, &
0.199253E+05,0.200178E+05,0.201105E+05,0.202036E+05,0.202971E+05, &
0.203910E+05,0.204853E+05,0.205798E+05,0.206749E+05,0.207701E+05, &
0.208659E+05,0.209621E+05,0.210584E+05,0.211554E+05,0.212524E+05, &
0.213501E+05,0.214482E+05,0.215465E+05,0.216452E+05,0.217442E+05/
data (ES(IES),IES=1521,1552) / &
0.218439E+05,0.219439E+05,0.220440E+05,0.221449E+05,0.222457E+05, &
0.223473E+05,0.224494E+05,0.225514E+05,0.226542E+05,0.227571E+05, &
0.228606E+05,0.229646E+05,0.230687E+05,0.231734E+05,0.232783E+05, &
0.233839E+05,0.234898E+05,0.235960E+05,0.237027E+05,0.238097E+05, &
0.239173E+05,0.240254E+05,0.241335E+05,0.242424E+05,0.243514E+05, &
0.244611E+05,0.245712E+05,0.246814E+05,0.247923E+05,0.249034E+05, &
0.250152E+05,0.250152E+05/

do I = 1, npnts

  ! Compute the factor that converts from sat vapour pressure in a
  ! pure water system to sat vapour pressure in air, FSUBW.  This formula
  ! is taken from equation A4.7 of Adrian Gill's Book: Atmosphere-Ocean
  ! Dynamics.  Note that his formula works in terms of pressure in MB and
  ! temperature in Celsius, so conversion of units leads to the slightly
  ! different equation used here.

  FSUBW = 1.0 + 1.0E-8 * P(I) * (4.5 + 6.0E-4 * (T(I) - t0c) * (T(I) - t0c))

  ! use the lookup table to find saturated vapour pressure, and store it in QS

  TT = max (T_low, T(I))
  TT = min (T_high, TT)

  ATABLE = (TT - T_low + delta_T) / delta_T
  ITABLE = ATABLE
  ATABLE = ATABLE - ITABLE

  QS(I) = (1.0 - ATABLE) * ES(ITABLE) + ATABLE * ES(ITABLE + 1)

  ! Multiple by FSUBW to convert to saturated vapour pressure in air
  ! (equation A4.6 of Adrian Gill's book)

  QS(I) = QS(I) * FSUBW

  ! Now form the accurate expression for QS, which is a rearranged version of
  ! equation A4.3 of Gill's book.
  !
  ! Note that at very low pressures we apply a fix, to prevent a singularity
  ! (Qsat tends to 1.0 kg/kg)

  QS(I) = (epsilon * QS(I)) / (max (P(I), QS(I)) - one_minus_epsilon * QS(I))

end do

end subroutine ufo_rttovonedvarcheck_QsatWat

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
subroutine ufo_rttovonedvarcheck_CheckIteration (geovals,    &
                                      profindex,  &
                                      nlevels_1dvar, &
                                      profile,    &
                                      OutOfRange)

implicit none

! subroutine arguments:
type(ufo_geovals), intent(in)    :: geovals
type(ufo_rttovonedvarcheck_profindex), intent(in) :: profindex
integer, intent(in)              :: nlevels_1dvar
real(kind_real), intent(inout)   :: profile(:)
logical, intent(out)             :: OutOfRange

! Local declarations:
real(kind_real), allocatable :: qsaturated(:)
real(kind_real), allocatable :: scaled_qsaturated(:)
real(kind_real)              :: q2_sat(1)
real(kind_real), allocatable :: Plevels_1DVar(:)
real(kind_real)              :: Pstar_Pa(1)
real(kind_real)              :: Temp(nlevels_1dvar)
real(kind_real)              :: Temp2(1)
real(kind_real)              :: rtbase
integer                      :: nlevels_q
integer                      :: toplevel_q
character(len=*), parameter  :: RoutineName = "ufo_rttovonedvarcheck_CheckIteration"
type(ufo_geoval), pointer    :: geoval
character(len=max_string)    :: varname
logical                      :: useRHwaterForQC
integer                      :: ii

! Setup
OutOfRange = .false.
useRHwaterForQC = .true.

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
  varname = "air_temperature"
  call ufo_geovals_get_var(geovals, varname, geoval)
  Temp = geoval%vals(:, 1) ! K
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
  varname = "air_temperature_at_two_meters_above_surface"
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
  varname = "air_pressure"
  call ufo_geovals_get_var(geovals, varname, geoval)
  if (.not. allocated(Plevels_1DVar)) allocate(Plevels_1DVar(nlevels_q))
  Plevels_1DVar(:) = geoval%vals(:, 1) ! K

  !----
  ! 2.1) Levels
  !----

  ! Note that if number of profile levels is less than number of pressure levels
  ! we assume the levels are from the surface upwards (remember that RTTOV levels
  ! are upside-down)

  if (profindex % q(1) > 0) then

    if (useRHwaterForQC) then
      call ufo_rttovonedvarcheck_QsatWat(qsaturated(1:nlevels_q), & ! out (qsat levels)
                                 Temp(1:nlevels_q),               & ! in  (t levels)
                                 Plevels_1DVar(1:nlevels_q),      & ! in  (p levels)
                                 nlevels_q)                         ! in
    else
      call ufo_rttovonedvarcheck_Qsat(qsaturated(1:nlevels_q), & ! out
                                 Temp(1:nlevels_q),            & ! in
                                 Plevels_1DVar(1:nlevels_q),   & ! in
                                 nlevels_q)                      ! in
    end if

    qsaturated(1:nlevels_q) = log (qsaturated(1:nlevels_q) * 1000.0)
    where (profile(profindex % q(1):profindex % q(2)) > qsaturated(1:nlevels_q))
      profile(profindex % q(1):profindex % q(2)) = qsaturated(1:nlevels_q)
    end where

  end if

  if (profindex % qt(1) > 0) then

    ! Qtotal is not included in the relative humidity quality control changes
    call ufo_rttovonedvarcheck_Qsat (qsaturated(1:nlevels_q), & ! out
                               Temp(1:nlevels_q),             & ! in
                               Plevels_1DVar(1:nlevels_q),    & ! in
                               nlevels_q)                       ! in

    ! scaled_qsaturated is generated as an upper limit on the value of qtot , and
    ! is set to 2*qsat.

    scaled_qsaturated(1:nlevels_q) = log (2 * qsaturated(1:nlevels_q) * 1000.0)
    where (profile(profindex % qt(1):profindex % qt(2)) > scaled_qsaturated(1:nlevels_q))
      profile(profindex % qt(1):profindex % qt(2)) = scaled_qsaturated(1:nlevels_q)
    end where

  end if

  !----
  ! 2.2) Surface
  !----

  if (profindex % q2 > 0) then
    varname = "air_pressure_at_two_meters_above_surface"
    call ufo_geovals_get_var(geovals, varname, geoval)
    Pstar_Pa(1) = geoval%vals(1, 1)
    if (useRHwaterForQC) then
      call ufo_rttovonedvarcheck_QsatWat (q2_sat(1:1), & ! out
                                 Temp2,                & ! in
                                 Pstar_Pa(1:1),        & ! in
                                 1)                      ! in
    else
      call ufo_rttovonedvarcheck_Qsat (q2_sat(1:1), & ! out
                                 Temp2,             & ! in
                                 Pstar_Pa(1:1),     & ! in
                                 1)                   ! in
    end if
    q2_sat(1) = log (q2_sat(1) * 1000.0)
    write(*,*) "profile(profindex % q2) = ",profile(profindex % q2)
    write(*,*) "q2_sat = ",q2_sat(1)
    if (profile(profindex % q2) > q2_sat(1)) then
      profile(profindex % q2) = q2_sat(1)
    end if
  end if

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
  OutOfRange )     ! out

implicit none

type(ufo_geovals), intent(in)    :: geovals
type(ufo_rttovonedvarcheck_profindex), intent(in) :: profindex
integer, intent(in)              :: nlevels_1dvar
logical, intent(out)             :: OutOfRange

! Local variables:
real(kind_real) :: LWP
real(kind_real) :: IWP
real(kind_real) :: dp
real(kind_real) :: meanql
real(kind_real) :: meanqi
integer         :: i
integer         :: nlevels_q, toplevel_q
character(len=*), parameter :: RoutineName = "ufo_rttovonedvarcheck_CheckCloudyIteration"

real(kind_real), parameter   :: MaxLWP = 2.0
real(kind_real), parameter   :: MaxIWP = 3.0
real(kind_real)              :: Plevels_1DVar(nlevels_1dvar)
type(ufo_geoval), pointer    :: geoval
real(kind_real)              :: clw(nlevels_1dvar)
real(kind_real)              :: ciw(nlevels_1dvar)
character(len=max_string)    :: varname

!-------------------------------------------------------------------------------

!initialise
OutOfRange = .false.
IWP = 0.0
LWP = 0.0

! Get pressure from geovals
varname = "air_pressure"
call ufo_geovals_get_var(geovals, varname, geoval)
Plevels_1DVar(:) = geoval%vals(:, 1) ! K

! Get clw from geovals
varname = "mass_content_of_cloud_liquid_water_in_atmosphere_layer"
call ufo_geovals_get_var(geovals, varname, geoval)
clw = geoval%vals(:, 1)

! Get ciw from geovals
varname = "mass_content_of_cloud_ice_in_atmosphere_layer"
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

if (any(ciw(:) > 0.0) .or. &
    any(clw(:) > 0.0)) then

!1.1 compute iwp, lwp

  do i=1, nlevels_1dvar-1

    dp =  Plevels_1DVar(i) - Plevels_1DVar(i+1)

    ! Calculate layer mean from CloudIce on levels
    meanqi = 0.5 * &
      (ciw(i) + ciw(i+1))
    if (meanqi > 0.0) then
      IWP = IWP + dp * meanqi
    end if

    ! Calculate layer mean from CLW on levels
    meanql = 0.5 * (clw(i) + clw(i+1))
    if (meanql > 0.0) then
      LWP = LWP + dp * meanql
    end if

  end do

  IWP = IWP / grav
  LWP = LWP / grav


!2.1 test if lwp iwp exceeds thresholds

  if ((IWP > MaxIWP) .or. (LWP > MaxLWP)) then
    write(*,*) "lwp or iwp exceeds thresholds"
    OutOfRange = .true.
    write(*,*) "lwp and iwp = ",LWP,IWP
  else
    write(*,*) "lwp and iwp less than thresholds"
    write(*,*) "lwp and iwp = ",LWP,IWP
  end if

end if

end subroutine ufo_rttovonedvarcheck_CheckCloudyIteration

!-------------------------------------------------------------------------------
!> Do cholesky decomposition
!!
!! \details Heritage: Ops_Cholesky.inc
!!
!! Solves the Linear equation UQ=V for Q where U is a symmetric positive definite
!! matrix and U and Q are vectors of length N.  The method follows that in Golub
!! and Van Loan although this is pretty standard.
!!
!! if U is not positive definite this will be detected by the program and flagged
!! as an error.  U is assumed to be symmetric as only the upper triangle is in
!! fact used.
!!
!! \author Met Office
!!
!! \date 09/06/2020: Created
!!
subroutine ufo_rttovonedvarcheck_Cholesky (U,         &
                         V,         &
                         N,         &
                         Q,         &
                         ErrorCode)

implicit none

! subroutine arguments:
integer, intent(in)          :: n
real(kind_real), intent(in)  :: U(n,n)
real(kind_real), intent(in)  :: V(n)
real(kind_real), intent(out) :: Q(n)
integer, intent(out)         :: ErrorCode

! Local declarations:
character(len=*), parameter  :: RoutineName = "ufo_rttovonedvarcheck_Cholesky"
real(kind_real), parameter   :: Tolerance = tiny (0.0) * 100.0
character(len=max_string)    :: ErrorMessage
integer                      :: j
integer                      :: k
real(kind_real)              :: G(n,n)   ! The Cholesky Triangle Matrix
real(kind_real)              :: X(n)     ! Temporary array used in calculating G

ErrorCode = 0

! Determine the Cholesky triangle matrix.

do j = 1, n
  X(j:n) = U(j:n,j)
  if (j /= 1) then
    do k = 1, j - 1
      X(j:n) = X(j:n) - G(j,k) * G(j:n,k)
    end do
  end if
  if (X(j) <= Tolerance) then
    ErrorCode = 1
    Errormessage = ' :U matrix is not positive definite'
    write(*,*) RoutineName,ErrorMessage
    goto 9999
  end if
  G(J:N,J) = X(J:N) / sqrt (X(J))
end do

! Solve Gx=v for x by forward substitution

X = V
X(1) = X(1) / G(1,1)
do j = 2, n
  X(j) = (X(j) - dot_product(G(j,1:j - 1), X(1:j - 1))) / G(j,j)
end do

! Solve G^T.q=x for q by backward substitution

Q = x
Q(n) = Q(n) / G(n,n)
do j = n - 1, 1, -1
  Q(j) = (Q(j) - dot_product(G(j + 1:n,j), Q(j + 1:n))) / G(j,j)
end do

9999 continue

end subroutine ufo_rttovonedvarcheck_Cholesky

! ----------------------------------------------------------

end module ufo_rttovonedvarcheck_minimize_utils_mod
