! (C) Copyright 2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module to provide code shared between nonlinear and tlm/adm radiance calculations

module ufo_onedvarfortran_forward_model_mod

use fckit_configuration_module, only: fckit_configuration
use iso_c_binding
use kinds
use ufo_geovals_mod
use ufo_onedvarfortran_utils_mod
use ufo_radiancerttov_tlad_mod

implicit none

private

! subroutines - all listed for complete
public ufo_onedvarfortran_ForwardModel
public ufo_onedvarfortran_GetHmatrixRTTOV

contains

!------------------------------------------------------------------------------

subroutine ufo_onedvarfortran_ForwardModel(geovals, ob_info, obsdb, &
                                         channels, conf, &
                                         profindex, prof_x, &
                                         hofx, H_matrix)

implicit none

! subroutine arguments
type(ufo_geovals), intent(in)               :: geovals
type(Obinfo_type), intent(in)               :: ob_info
type(c_ptr), value, intent(in)              :: obsdb
integer(c_int), intent(in)                  :: channels(:)
type(c_ptr), value, intent(in)              :: conf
type(Profileinfo_type), intent(in)          :: profindex
real(kind_real), intent(in)                 :: prof_x(:)    ! x vector for 1D-var
real(kind_real), intent(out)                :: hofx(:)
real(kind_real), intent(out)                :: H_matrix(:,:)

! local variables
type(ufo_radiancerttov_tlad)                :: rttov_tlad      ! structure for holding original rttov_k setup data
type(fckit_configuration)                   :: f_conf

f_conf = fckit_configuration(conf)

select case (trim(ob_info%forward_mod_name))
  case ("RTTOV")
    call rttov_tlad % setup(f_conf)
    call ufo_onedvarfortran_GetHmatrixRTTOV(geovals, ob_info, obsdb, &
                                 channels(:), rttov_tlad, &
                                 profindex, prof_x(:), &
                                 hofx(:), H_matrix)
    call rttov_tlad % delete()

   case default
      write(*,*) "No suitable forward model exiting"
      stop
end select

end subroutine ufo_onedvarfortran_ForwardModel

!------------------------------------------------------------------------------

subroutine ufo_onedvarfortran_GetHmatrixRTTOV(geovals, ob_info, obsdb, &
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
type(Profileinfo_type), intent(in)          :: profindex
real(kind_real), intent(in)                 :: prof_x(:)    ! x vector for 1D-var
real(kind_real), intent(out)                :: hofx(:)
real(kind_real), intent(out)                :: H_matrix(:,:)

! Local arguments
integer :: nchans, nlevels, nq_levels
integer :: i
logical :: RTTOV_GasUnitConv = .false.
real(kind_real),allocatable :: q_kgkg(:)
real(kind_real)             :: s2m_kgkg

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
    H_matrix(i,profindex % t(1):profindex % t(2)) = rttov_data % profiles_k(1) % t(nlevels:1:-1)
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

  q_kgkg(:) = exp(prof_x(profindex % q(1):profindex % q(2))) / 1000.0

  do i = 1, nchans
    H_matrix(i,profindex % q(1):profindex % q(2)) = &
      rttov_data % profiles_k(1) % q(nlevels:1:-1) * q_kgkg(:)
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

  nq_levels = profindex % qt(2)-profindex % qt(1)+1
  allocate(q_kgkg(nq_levels))

  q_kgkg(:) = exp(prof_x(profindex % qt(1):profindex % qt(2))) / 1000.0

  do i = 1, nchans
    H_matrix(i,profindex % qt(1):profindex % qt(2)) = &
      rttov_data % profiles_k(1) % q(nlevels:1:-1) * q_kgkg(:)
  end do

  deallocate(q_kgkg)

end if

!----
! 2.) Single-valued variables
!----

! 2.1) Surface Temperature

if (profindex % t2 > 0) then
  H_matrix(:,profindex % t2) = rttov_data % profiles_k(1) % s2m % t
end if

! 2.2) Water vapour

if (profindex % q2 > 0) then

  s2m_kgkg = exp(prof_x(profindex % q2)) / 1000.0
  H_matrix(:,profindex % q2) = rttov_data % profiles_k(1) % s2m % q * s2m_kgkg

end if

! 2.3) Surface pressure

if (profindex % pstar > 0) then
  H_matrix(:,profindex % pstar) = rttov_data % profiles_k(1) % s2m % p
end if

! 2.4) Skin temperature

if (profindex % tstar > 0) then
  H_matrix(:,profindex % tstar) = rttov_data % profiles_k(1) % skin % t
end if

end subroutine ufo_onedvarfortran_GetHmatrixRTTOV

!---------------------------------------------------------------------------

end module ufo_onedvarfortran_forward_model_mod
