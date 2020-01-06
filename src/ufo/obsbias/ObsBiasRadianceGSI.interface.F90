! (C) Copyright 2019 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module to handle radiancecrtm observations

module ufo_obsbiasradiancegsi_utils_c

  use iso_c_binding
  use clw_mod, only : calc_clw

  implicit none
  private

  ! ------------------------------------------------------------------------------

contains

! ------------------------------------------------------------------------------

subroutine calc_clw_c(nadir, tb_obs, tsim, ich, nchanl, &
                      no85GHz, amsua, ssmi, ssmis, amsre, &
                      atms, amsr2, gmi, saphir, &
                      tsavg5, sfc_speed, zasat, clw) &
                 bind(c,name='calc_clw_f90')
implicit none

integer(c_int),      intent(in   ) :: nadir
integer(c_int),      intent(in   ) :: nchanl
real(c_float),       intent(in   ) :: tb_obs(nchanl), tsim(nchanl)
integer(c_int),      intent(in   ) :: ich(nchanl)
logical(c_bool),     intent(in   ) :: no85GHz, amsua, ssmi, ssmis, amsre,&
                                      atms, amsr2, gmi, saphir
real(c_float),       intent(in   ) :: tsavg5,sfc_speed,zasat
real(c_float),       intent(  out) :: clw

!> Local variables

real(c_float)  :: tpwc, gwp
integer(c_int) :: kraintype, ierrret

call calc_clw(nadir, tb_obs, tsim, ich, nchanl, &
              no85GHz, amsua, ssmi, ssmis, amsre, &
              atms, amsr2, gmi, saphir, &
              tsavg5, sfc_speed, zasat, &
              clw, tpwc, gwp, kraintype, ierrret)

end subroutine calc_clw_c

! ------------------------------------------------------------------------------

end module ufo_obsbiasradiancegsi_utils_c
