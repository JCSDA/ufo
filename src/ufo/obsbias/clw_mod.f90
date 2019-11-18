!
! (C) Copyright 2019 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Fortran module adopted from GSI to estimates cloud liquid water for micro.

module clw_mod

use iso_c_binding, only : c_bool
use kinds, only: r_kind => kind_float, i_kind => kind_int

implicit none
! set default to private
  private
! set routines used externally to public
  public :: calc_clw

contains


 subroutine calc_clw(nadir,tb_obs,tsim,ich,nchanl,no85GHz,amsua,ssmi,ssmis,amsre,atms, &   
          amsr2,gmi,saphir,tsavg5,sfc_speed,zasat,clw,tpwc,gwp,kraintype,ierrret)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:   calc_clw    estimates cloud liquid water for micro. QC
!   prgmmr: derber           org: np23                date: 1995-07-06
!
! abstract: estimates cloud liquid water for microwave quality control and
!        bias correction.
!
!  input argument list:
!     nadir     - scan position
!     tb_obs    - observed brightness temperatures
!     tsim      - simulated brightness temperatures             
!     ich       - channel number array
!     nchanl    - number of channels    
!     no85ghz   - flag for instrument with no 85ghz channel   
!     amsua     - flag for amsua data
!     ssmi      - flag for ssmi  data
!     ssmis     - flag for ssmis data
!     amsre     - flag for amsre data
!     atms      - flag for atms data
!     amsr2     - flag for amsr2 data
!     gmi       - flag for gmi data
!     saphir    - flag for saphir data
!     tsavg5    - Surface temperature value
!     sfc_speed - surface wind speed (10m)
!     zasat     - satellite zenith angle
!
!   output argument list:
!     clw       - cloud liquid water
!     gwp       - graupel water path                                                   
!     tpwc      - total column water vapor                                           
!     kraintype - rain type
!     ierrret   - return flag
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$
!  use radinfo, only: ang_rad,cbias,air_rad,predx,adp_anglebc
  use constants, only: zero,one,amsua_clw_d1,amsua_clw_d2,t0c,r1000

  integer(i_kind)                   ,intent(in   ) :: nadir,nchanl
  real(r_kind),dimension(nchanl)    ,intent(in   ) :: tb_obs,tsim
  integer(i_kind),dimension(nchanl) ,intent(in   ) :: ich
  logical(c_bool)                   ,intent(in   ) :: no85GHz,amsre,ssmi,ssmis,amsua,atms,amsr2,gmi,saphir
  real(r_kind)                      ,intent(in   ) :: tsavg5,sfc_speed,zasat
  real(r_kind)                      ,intent(  out) :: clw,tpwc,gwp
  integer(i_kind)                   ,intent(  out) :: kraintype,ierrret


! Declare local parameters
  real(r_kind),parameter:: r284=284.0_r_kind
  real(r_kind),parameter:: r285=285.0_r_kind

! Declare local variables
  real(r_kind) tbcx1,tbcx2
  logical adp_anglebc ! logical to turn off or on the variational radiance angle bias correction

! For simplicity
  adp_anglebc = .true. ! true.=turn on angle bias correction

  if (amsua .or. atms) then

    ! We want to reject sea ice points that may be frozen.  The sea freezes
    ! around -1.9C but we set the threshold at 1C to be safe.
     if(tsavg5>t0c-one .and. tb_obs(1) > zero .and. tb_obs(2) > zero) then 
        if (adp_anglebc) then
           tbcx1=tsim(1) !+cbias(nadir,ich(1))*ang_rad(ich(1))+predx(1,ich(1))*air_rad(ich(1))
           tbcx2=tsim(2) !+cbias(nadir,ich(2))*ang_rad(ich(2))+predx(1,ich(2))*air_rad(ich(2))
        else
           tbcx1=tsim(1) !+cbias(nadir,ich(1))*ang_rad(ich(1))
           tbcx2=tsim(2) !+cbias(nadir,ich(2))*ang_rad(ich(2))
        end if
        if (tbcx1 <=r284 .and. tbcx2<=r284 .and. tb_obs(1) > zero &
            .and. tb_obs(2) > zero) then 
             clw=amsua_clw_d1*(tbcx1-tb_obs(1))/(r285-tbcx1)+ &
                 amsua_clw_d2*(tbcx2-tb_obs(2))/(r285-tbcx2)
             ierrret = 0
        else
             ierrret = 1
        endif
     else
        clw = r1000
        ierrret = 1
     end if
     
     if (.not. adp_anglebc) clw = max(zero,clw)

  endif

  return
 end subroutine calc_clw

end module clw_mod
