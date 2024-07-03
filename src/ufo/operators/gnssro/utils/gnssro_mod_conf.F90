!==========================================================================
module gnssro_mod_conf
!==========================================================================
use fckit_configuration_module, only: fckit_configuration
use iso_c_binding
use kinds
use obsspace_mod
use ufo_constants_mod, only: mean_earth_rad
use gnssro_mod_constants, only: MAXVARLEN

implicit none

private
public   :: gnssro_conf
public   :: gnssro_conf_setup

type gnssro_conf
  integer(c_int)     :: use_compress
  integer(c_int)     :: n_horiz
  integer(c_int)     :: sr_steps
  integer(c_int)     :: modeltop
  integer(c_int)     :: nlevadd
  integer(c_int)     :: ngrd
  character(len=MAXVARLEN)      :: super_ref_qc
  character(len=:), allocatable :: str
  real(kind_real)    :: res
  real(kind_real)    :: top_2d
  character(len=20)  :: ro_type
  real(kind_real)    :: dtheta
  character(len=20)  :: vertlayer
  character(len=20)  :: output_diags
  character(len=20)  :: modeltopconfig
  character(len=20)  :: GSI_version
end type gnssro_conf

!--------- ropp2d location default parameters-----------------
integer(c_int),  parameter, public :: use_compress = 0   !0 or 1
integer(c_int),  parameter, public :: n_horiz_2d = 31   !should be odd number
real(kind_real), parameter, public :: res_2d     = 40.0 !km
real(kind_real), parameter, public :: top_2d     = 20.0 !km; maximum height the 2d operator is applied
character(len=20), parameter, public :: ro_type  = "spaceborne" !takes values spaceborne/airborne
contains
!-------------------------------

subroutine gnssro_conf_setup(roconf, f_conf)
implicit none
type(gnssro_conf), intent(inout)      :: roconf
type(fckit_configuration), intent(in) :: f_conf

character(len=:), allocatable :: str

roconf%use_compress = use_compress !0: Coefficient to calculate dry refractivity contribution for ideal gas
                       !1: non ideal gas, use adjusted coefficients for refractivity computation
if (f_conf%has("use_compress")) call f_conf%get_or_die("use_compress",roconf%use_compress)
roconf%n_horiz = n_horiz_2d
if (f_conf%has("n_horiz")) call f_conf%get_or_die("n_horiz",roconf%n_horiz)
roconf%res = res_2d
if (f_conf%has("res")) call f_conf%get_or_die("res",roconf%res)
roconf%top_2d = top_2d
if (f_conf%has("top_2d")) call f_conf%get_or_die("top_2d",roconf%top_2d)
roconf%top_2d        = roconf%top_2d*1000.0     ! km to m
roconf%ro_type = ro_type
if (f_conf%has("ro_type"))  then
  call f_conf%get_or_die("ro_type",str)
  roconf%ro_type=trim(str)
endif
!if (f_conf%has("ro_type")) call f_conf%get_or_die("ro_type",roconf%ro_type)
roconf%dtheta        = roconf%res/mean_earth_rad
roconf%vertlayer = "full"
if (f_conf%has("vertlayer")) then
   call f_conf%get_or_die("vertlayer",str)
   roconf%vertlayer = str
end if
roconf%super_ref_qc = "NBAM"
if (f_conf%has("super_ref_qc"))  then
  call f_conf%get_or_die("super_ref_qc",str)
  roconf%super_ref_qc=trim(str)
endif
roconf%sr_steps = 2
if (f_conf%has("sr_steps")) call f_conf%get_or_die("sr_steps",roconf%sr_steps)
roconf%output_diags="false"
if (f_conf%has("output_diags"))  then
  call f_conf%get_or_die("output_diags",str)
  roconf%output_diags=trim(str)
endif
roconf%GSI_version="EMC"
if (f_conf%has("GSI_version"))  then
  call f_conf%get_or_die("GSI_version",str)
  roconf%GSI_version=trim(str)
endif
roconf%modeltopconfig="false"
if (f_conf%has("modeltop")) then  !an integer specifying roughly the model top in km
  call f_conf%get_or_die("modeltop", roconf%modeltop)
  if (roconf%modeltop .gt. 0)   roconf%modeltopconfig="true"
endif
if (f_conf%has("nlevadd")) call f_conf%get_or_die("nlevadd", roconf%nlevadd)
if (f_conf%has("ngrd")) call f_conf%get_or_die("ngrd", roconf%ngrd)
end subroutine gnssro_conf_setup

end module gnssro_mod_conf
