! (C) Copyright 2020 Met Office
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module to provide code shared between nonlinear and tlm/adm radiance calculations

module ufo_rttovonedvarcheck_get_jacobian_mod

use fckit_configuration_module, only: fckit_configuration
use iso_c_binding
use kinds
use ufo_geovals_mod
use ufo_radiancerttov_tlad_mod
use ufo_rttovonedvarcheck_utils_mod
use ufo_rttovonedvarcheck_minimize_utils_mod, only: ufo_rttovonedvarcheck_Qsplit
use ufo_rttovonedvarcheck_profindex_mod, only: profindex_type

implicit none
private

! subroutines - all listed for complete
public ufo_rttovonedvarcheck_get_jacobian

contains

!------------------------------------------------------------------------------

subroutine ufo_rttovonedvarcheck_get_jacobian(geovals, ob_info, obsdb, &
                                           channels, conf, &
                                           profindex, prof_x, &
                                           hofx, H_matrix)

implicit none

! subroutine arguments
type(ufo_geovals), intent(in)               :: geovals
type(Obinfo_type), intent(in)               :: ob_info
type(c_ptr), value, intent(in)              :: obsdb
integer(c_int), intent(in)                  :: channels(:)
type(fckit_configuration), intent(in)       :: conf
type(profindex_type), intent(in)          :: profindex
real(kind_real), intent(in)                 :: prof_x(:)    ! x vector for 1D-var
real(kind_real), intent(out)                :: hofx(:)
real(kind_real), intent(out)                :: H_matrix(:,:)

! local variables
type(ufo_radiancerttov_tlad)                :: rttov_tlad   ! structure for holding original rttov_k setup data

select case (trim(ob_info % forward_mod_name))
  case ("RTTOV")
    call rttov_tlad % setup(conf)
    call ufo_rttovonedvarcheck_GetHmatrixRTTOV(geovals, ob_info, obsdb, &
                                            channels(:), rttov_tlad, &
                                            profindex, prof_x(:), &
                                            hofx(:), H_matrix) ! out
    call rttov_tlad % delete()

   case default
      write(*,*) "No suitable forward model exiting"
      stop
end select

end  subroutine ufo_rttovonedvarcheck_get_jacobian

!------------------------------------------------------------------------------

subroutine ufo_rttovonedvarcheck_GetHmatrixRTTOV(geovals, ob_info, obsdb, &
                                       channels, rttov_data, &
                                       profindex, prof_x, &
                                       hofx, H_matrix)

implicit none

! subroutine arguments
type(ufo_geovals), intent(in)               :: geovals
type(Obinfo_type), intent(in)               :: ob_info
type(c_ptr), value, intent(in)              :: obsdb
integer(c_int), intent(in)                  :: channels(:)
type(ufo_radiancerttov_tlad), intent(inout) :: rttov_data      ! structure for running rttov_k
type(profindex_type), intent(in)          :: profindex
real(kind_real), intent(in)                 :: prof_x(:)    ! x vector for 1D-var
real(kind_real), intent(out)                :: hofx(:)
real(kind_real), intent(out)                :: H_matrix(:,:)

! Local arguments
integer :: nchans, nlevels, nq_levels
integer :: i
logical :: RTTOV_GasunitConv = .false.
real(kind_real),allocatable :: q_kgkg(:)
real(kind_real)             :: s2m_kgkg
type(ufo_geoval), pointer    :: geoval
real(kind_real), allocatable :: temperature(:)
real(kind_real), allocatable :: pressure(:)
real(kind_real), allocatable :: dq_dqt(:)
real(kind_real), allocatable :: dql_dqt(:)
real(kind_real), allocatable :: dqi_dqt(:)

! call rttov k code?
!call rttov_data%settraj(geovals, obsdb, channels, obs_info=ob_info, Hx=H_matrix, BT=hofx)
call rttov_data%settraj(geovals, obsdb, channels, obs_info=ob_info, BT=hofx)

nchans = size(channels)
nlevels = size(rttov_data % profiles_k(1) % t)

! ----------------
! output new hofx
! ---------------

!hofx(:) = rttov_data % bt(:)

!----------------------------------------------------------------
! 1.1) Temperature - invert as RTTOV level 1 as top of atmosphere and
!               1Dvar profile is from the surface up. 
!----------------------------------------------------------------

if (profindex % t(1) > 0) then
  do i = 1, nchans
    H_matrix(i,profindex % t(1):profindex % t(2)) = rttov_data % profiles_k(i) % t(nlevels:1:-1)
  end do
end if

!------
! 1.2) Water vapour
!------

! Water Vapour Jacobians must be converted from
! kg/kg to ln(g/kg) - the unit conversion cancels, then:
! dy/d(ln q) = dy/dq * q(kg/kg)

if (profindex % q(1) > 0) then

  nq_levels = profindex % q(2)-profindex % q(1)+1
  allocate(q_kgkg(nq_levels))

  q_kgkg(:) = exp(prof_x(profindex % q(1):profindex % q(2))) / 1000.0_kind_real

  do i = 1, nchans
    H_matrix(i,profindex % q(1):profindex % q(2)) = &
      rttov_data % profiles_k(i) % q(nlevels:1:-1) * q_kgkg(:)
  end do

  deallocate(q_kgkg)

end if

!------
! 1.3) Total water
!------

! For the sake of this first implementation this will not include liquid
! and ice water content just water vapour which is consistent with the
! profile loaded from GeoVaLs.

if (profindex % qt(1) > 0) then

  allocate(q_kgkg(nlevels))
  allocate(temperature(nlevels))
  allocate(pressure(nlevels))
  allocate(dq_dqt(nlevels))
  allocate(dql_dqt(nlevels))
  allocate(dqi_dqt(nlevels))

  ! Get qtotal from profile
  q_kgkg(:) = exp(prof_x(profindex % qt(1):profindex % qt(2))) / 1000.0_kind_real

  ! Get temperature and pressure from geovals
  call ufo_geovals_get_var(geovals, "air_temperature", geoval)
  temperature(:) = geoval%vals(:, 1) ! K
  call ufo_geovals_get_var(geovals, "air_pressure", geoval)
  pressure(:) = geoval%vals(:, 1)    ! Pa

  ! Calculate the gradients with respect to qtotal
  call ufo_rttovonedvarcheck_Qsplit (0,               & ! in
                          temperature,     & ! in
                          pressure,        & ! in
                          nlevels,         & ! in
                          q_kgkg(:),       & ! in
                          dq_dqt(:),       & ! out
                          dql_dqt(:),      & ! out
                          dqi_dqt(:))        ! out

  ! Calculate jacobnian wrt humidity and clw
  do i = 1, nchans
    H_matrix(i,profindex % qt(1):profindex % qt(2)) = &
            (rttov_data % profiles_k(i) % q(nlevels:1:-1) * dq_dqt(:) + &
             rttov_data % profiles_k(i) % clw(nlevels:1:-1) * dql_dqt(:) ) * q_kgkg(:)
  end do

  ! Clean up
  deallocate(q_kgkg)
  deallocate(temperature)
  deallocate(pressure)
  deallocate(dq_dqt)
  deallocate(dql_dqt)
  deallocate(dqi_dqt)

end if

!----
! 2.) Single-valued variables
!----

! 2.1) Surface Temperature

if (profindex % t2 > 0) then
  do i = 1, nchans
    H_matrix(i,profindex % t2) = rttov_data % profiles_k(i) % s2m % t
  end do
end if

! 2.2) Water vapour

if (profindex % q2 > 0) then

  s2m_kgkg = exp(prof_x(profindex % q2)) / 1000.0_kind_real
  do i = 1, nchans
    H_matrix(i,profindex % q2) = rttov_data % profiles_k(i) % s2m % q * s2m_kgkg
  end do

end if

! 2.3) Surface pressure

if (profindex % pstar > 0) then
  do i = 1, nchans
    H_matrix(i,profindex % pstar) = rttov_data % profiles_k(i) % s2m % p
  end do
end if

! 2.4) Skin temperature

if (profindex % tstar > 0) then
  do i = 1, nchans
    H_matrix(i,profindex % tstar) = rttov_data % profiles_k(i) % skin % t
  end do
end if

call ufo_rttovonedvarcheck_PrintHmatrix( &
  nchans,   &       ! in
  size(prof_x),  &  ! in
  Channels, &       ! in
  H_matrix, &       ! in
  profindex )       ! in

end  subroutine ufo_rttovonedvarcheck_GetHmatrixRTTOV

!---------------------------------------------------------------------------

subroutine ufo_rttovonedvarcheck_PrintHmatrix( &
  nchans,   &       ! in
  nprofelements, &  ! in
  Channels, &       ! in
  H_matrix, &       ! in
  profindex )       ! in

implicit none

!Subroutine arguments:
integer,                intent(IN) :: nchans
integer,                intent(IN) :: nprofelements
integer(c_int),         intent(IN) :: channels(nchans)
real(kind_real),        intent(IN) :: H_matrix(nchans,nprofelements)
type(profindex_type), intent(in)   :: profindex

! Local variables:
integer :: i
character(LEN=10) :: int_fmt
character(LEN=12) :: real_fmt
character(LEN=3) :: txt_nchans
character(LEN=*), parameter :: RoutineName = "ufo_rttovonedvarcheck_PrintHmatrix"
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

end module ufo_rttovonedvarcheck_get_jacobian_mod
