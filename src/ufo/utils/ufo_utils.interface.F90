!-------------------------------------------------------------------------------
! (C) Crown Copyright 2021 Met Office
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
!-------------------------------------------------------------------------------

module ufo_cloudcostfunction_mod_c

use iso_c_binding
use kinds
use ufo_utils_mod, only: Ops_SatRad_Qsplit, Ops_QsatWat, Ops_Qsat

implicit none

private

contains

!-------------------------------------------------------------------------------
! Met Office Qsplit routine
! output_type = 1: Partition total humidity into water vapour, liquid water and ice
! output_type != 1: Compute derivatives dq/dqtotal, dql/dqtotal, dqi/dqtotal
subroutine ufo_ops_satrad_qsplit_c(output_type, nvals, p_in, t_in, qtotal_in, &
                                   q_out, ql_out, qi_out, use_qt_split_rain) &
                                   bind(C, name='ufo_ops_satrad_qsplit_f90')

implicit none
integer(c_int), intent(in) :: output_type
integer(c_int), intent(in) :: nvals
real(c_float), intent(in)  :: p_in(nvals)
real(c_float), intent(in)  :: t_in(nvals)
real(c_float), intent(in)  :: qtotal_in(nvals)
real(c_float), intent(out) :: q_out(nvals)
real(c_float), intent(out) :: ql_out(nvals)
real(c_float), intent(out) :: qi_out(nvals)
logical(c_bool), intent(in) :: use_qt_split_rain

real(kind_real) :: q_out_kind_real(nvals)
real(kind_real) :: ql_out_kind_real(nvals)
real(kind_real) :: qi_out_kind_real(nvals)

call Ops_SatRad_Qsplit(output_type, dble(p_in), dble(t_in), dble(qtotal_in), &
                       q_out_kind_real, ql_out_kind_real, qi_out_kind_real, &
                       logical(use_qt_split_rain))

q_out = real(q_out_kind_real)
ql_out = real(ql_out_kind_real)
qi_out = real(qi_out_kind_real)

end subroutine ufo_ops_satrad_qsplit_c

!-------------------------------------------------------------------------------
! Met Office Ops_QsatWat routine
! Returns a saturation mixing ratio (kg/kg) given a temperature and pressure
! using saturation vapour pressure calculated using Goff-Gratch formulae.
! Values in the lookup table are over water above and below 0 degrees C.
subroutine ufo_ops_satrad_qsatwat(qs_out, t_in, p_in, nvals) &
                                  bind(C, name='ufo_ops_satrad_qsatwat_f90')

implicit none
integer(c_int), intent(in) :: nvals
real(c_float), intent(out) :: qs_out(nvals) ! saturation mixing ratio (kg/kg)
real(c_float), intent(in)  :: t_in(nvals)   ! temperature (K)
real(c_float), intent(in)  :: p_in(nvals)   ! pressure (Pa)

real(kind_real) :: qs_out_kind_real(nvals)

call Ops_QsatWat(qs_out_kind_real, dble(t_in), dble(p_in), nvals)

qs_out = real(qs_out_kind_real)

end subroutine ufo_ops_satrad_qsatwat

!-------------------------------------------------------------------------------
! Met Office Ops_Qsat routine
! Returns a saturation mixing ratio (kg/kg) given a temperature and pressure
! using saturation vapour pressure calculated using Goff-Gratch formulae.
! Values in the lookup table are over water above 0 degrees C and over
! ice below this temperature.
subroutine ufo_ops_qsat(qs_out, t_in, p_in, nvals) &
     bind(C, name='ufo_ops_qsat_f90')

implicit none
integer(c_int), intent(in) :: nvals
real(c_float), intent(out) :: qs_out(nvals) ! saturation mixing ratio (kg/kg)
real(c_float), intent(in)  :: t_in(nvals)   ! temperature (K)
real(c_float), intent(in)  :: p_in(nvals)   ! pressure (Pa)

real(kind_real) :: qs_out_kind_real(nvals)

call Ops_Qsat(qs_out_kind_real, dble(t_in), dble(p_in), nvals)

qs_out = real(qs_out_kind_real)

end subroutine ufo_ops_qsat

!-------------------------------------------------------------------------------
end module ufo_cloudcostfunction_mod_c
