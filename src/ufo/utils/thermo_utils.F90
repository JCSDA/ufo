! (C) Copyright 2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Fortran module for thermodynamic computations and conversions for use in UFO 

module thermo_utils_mod

use kinds, only: kind_real

implicit none
public

contains

! ------------------------------------------------------------------------------

subroutine calc_theta(t_in, p_in, t_out)
  ! compute potential (virtual) temperature from a given (virtual) temperature
  ! and pressure value, assumes the standard definition of theta at 1000 hPa
  ! units must be in K and Pa!
  use ufo_constants_mod, only: rd_over_cp
  implicit none
  real(kind_real), intent(in) :: t_in, p_in
  real(kind_real), intent(out) :: t_out

  t_out = t_in * (1.0e5_kind_real / p_in) ** rd_over_cp

end subroutine calc_theta

! ------------------------------------------------------------------------------

subroutine gsi_tp_to_qs( t, p, es, qs)
  ! calculate saturation specific humidity for a given
  ! temperature and pressure
  ! based on subroutin DA_TP_To_Qs in GSI
   use ufo_constants_mod, only: t0c, rd_over_rv, es_w_0

   implicit none
   real(kind_real), intent(in) :: t                ! Temperature.
   real(kind_real), intent(in) :: p             ! Pressure.
   real(kind_real), intent(out) :: es          ! Sat. vapour pressure.
   real(kind_real), intent(out) :: qs              ! Sat. specific humidity.

!  Saturation Vapour Pressure Constants(Rogers & Yau, 1989)
   real(kind_real), parameter    :: es_beta = 17.67_kind_real
   real(kind_real), parameter    :: es_gamma = 243.5_kind_real

   real(kind_real) :: omeps
   real(kind_real) :: t_c  ! T in degreesC.

   omeps = 1.0_kind_real - rd_over_rv
   t_c = t - t0c

   es = es_w_0 * exp( es_beta * t_c / ( t_c + es_gamma ) )

   qs = rd_over_rv * es / ( p - omeps * es )

   return
end subroutine gsi_tp_to_qs

! ------------------------------------------------------------------------------

end module thermo_utils_mod
