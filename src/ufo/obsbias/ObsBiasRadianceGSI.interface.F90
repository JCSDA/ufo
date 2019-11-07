! (C) Copyright 2019 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module to handle radiancecrtm observations

module ufo_obsbiasradiancegsi_utils_c

  use iso_c_binding
  use clw_mod, only : calc_clw
  use tlap_mod, only : calc_tlap

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

subroutine calc_tlap_c(newpc4pred, nsig, nchanl, &
                       ptau5, tsavg5, tvp, tlapmean, tlap) &
                  bind(c,name='calc_tlap_f90')
implicit none

logical(c_bool)                       ,intent(in   ) :: newpc4pred
integer(c_int)                        ,intent(in   ) :: nsig, nchanl
real(c_float), dimension(nsig,nchanl) ,intent(in   ) :: ptau5
real(c_float),                         intent(in   ) :: tsavg5
real(c_float), dimension(nsig)        ,intent(in   ) :: tvp
real(c_float), dimension(nchanl)      ,intent(in   ) :: tlapmean
real(c_float), dimension(nchanl)      ,intent(  out) :: tlap

call calc_tlap(newpc4pred, nsig, nchanl, &
               ptau5, tsavg5, tvp, tlapmean, tlap) 

end subroutine calc_tlap_c

! ------------------------------------------------------------------------------

end module ufo_obsbiasradiancegsi_utils_c
