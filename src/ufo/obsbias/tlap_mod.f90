!
! (C) Copyright 2019 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Fortran module adopted from GSI to calculate laps rate

module tlap_mod

use iso_c_binding, only : c_bool

implicit none
! set default to private
  private
! set routines used externally to public
  public :: calc_tlap

contains


 subroutine calc_tlap(newpc4pred, nsig, nchanl, &
                      ptau5, tsavg5, tvp, tlapmean, tlap)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:   calc_tlap
!   prgmmr: Xin ZHang           org: JCSDA                date: 2019-10-01
!
! abstract: adopt from GSI setuprad.f90.
!
!
!  input argument list:
!     ich       - channel number array
!     nchanl    - number of channels    
!
!   output argument list:
!     tlap      - tlap
!
! attributes:
!   language: f90
!   machine:  Linux/Gfortran
!
!$$$
  use kinds, only: r_kind => kind_float, i_kind => kind_int
  use constants, only: r0_01

  implicit none

  logical(c_bool)                      ,intent(in   ) :: newpc4pred
  integer(i_kind)                      ,intent(in   ) :: nsig, nchanl
  real(r_kind), dimension(nsig,nchanl) ,intent(in   ) :: ptau5
  real(r_kind),                         intent(in   ) :: tsavg5
  real(r_kind), dimension(nsig)        ,intent(in   ) :: tvp
  real(r_kind), dimension(nchanl)      ,intent(in   ) :: tlapmean
  real(r_kind), dimension(nchanl)      ,intent(  out) :: tlap

! Declare local variables
  real(r_kind)                    :: ptau5deriv, ptau5derivmax
  real(r_kind), dimension(nchanl) :: tlapchn
  integer(i_kind)                 :: k, i

  do i=1,nchanl

    tlapchn(i)= (ptau5(2,i)-ptau5(1,i))*(tsavg5-tvp(2))
    do k=2,nsig-1
       tlapchn(i)=tlapchn(i)+&
          (ptau5(k+1,i)-ptau5(k,i))*(tvp(k-1)-tvp(k+1))
    end do
    if (.not. newpc4pred) tlapchn(i) = r0_01*tlapchn(i)
    tlap(i) = tlapchn(i)-tlapmean(i)

  end do

  return
 end subroutine calc_tlap


end module tlap_mod
