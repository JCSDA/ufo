! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module for seaicethickness tl/ad observation operator

module ufo_seaicethickness_tlad_mod

 use fckit_configuration_module, only: fckit_configuration
 use iso_c_binding
 use kinds

 use ufo_geovals_mod, only: ufo_geovals, ufo_geoval, ufo_geovals_get_var
 use ufo_basis_tlad_mod, only: ufo_basis_tlad
 use ufo_vars_mod
 use obsspace_mod
 use missing_values_mod
 use obs_variables_mod
 use ufo_utils_mod, only: cmp_strings

 implicit none
 private

 integer, parameter :: max_string=800

 !> Fortran derived type for the tl/ad observation operator
 type, extends(ufo_basis_tlad), public :: ufo_seaicethickness_tlad
 private
  type(obs_variables), public :: obsvars
  character(max_string) :: thickness_sim_option
  type(ufo_geoval) :: icethick !< ice thickness (traj)
  type(ufo_geoval) :: icefrac  !< ice fraction  (traj)
  type(ufo_geoval) :: snowthick!< snow thickness(traj) 
  real(kind=kind_real) :: rho_ice  = 905.0 !< [kg/m3]
  real(kind=kind_real) :: rho_snow = 330.0 !< [kg/m3]
  real(kind=kind_real) :: rho_water= 1000.0!< [kg/m3]

 contains
  procedure :: setup  => ufo_seaicethickness_tlad_setup
  procedure :: delete  => ufo_seaicethickness_tlad_delete
  procedure :: settraj => ufo_seaicethickness_tlad_settraj
  procedure :: simobs_tl  => ufo_seaicethickness_simobs_tl
  procedure :: simobs_ad  => ufo_seaicethickness_simobs_ad
 end type ufo_seaicethickness_tlad

contains

! ------------------------------------------------------------------------------
subroutine ufo_seaicethickness_tlad_setup(self, f_conf)
implicit none
class(ufo_seaicethickness_tlad), intent(inout) :: self
type(fckit_configuration),       intent(in)    :: f_conf
real(kind=kind_real) :: rho_ice, rho_snow, rho_water
integer :: ivar, nvars
character(max_string)  :: err_msg

nvars = self%obsvars%nvars()
if (nvars /= 1) then
  write(err_msg,*) 'ufo_seaicethickness_tlad_setup error: only variables size 1 supported!'
  call abor1_ftn(err_msg)
endif

! Set thickness-simulate option from ymal file
!self%thickness_sim_option = self%obsvars%variable(1)

end subroutine ufo_seaicethickness_tlad_setup

! ------------------------------------------------------------------------------
subroutine ufo_seaicethickness_tlad_delete(self)
implicit none
class(ufo_seaicethickness_tlad), intent(inout) :: self

self%ltraj = .false.

end subroutine ufo_seaicethickness_tlad_delete

! ------------------------------------------------------------------------------
subroutine ufo_seaicethickness_tlad_settraj(self, geovals, obss)
implicit none
class(ufo_seaicethickness_tlad), intent(inout) :: self
type(ufo_geovals),       intent(in)    :: geovals
type(c_ptr), value,      intent(in)    :: obss

character(len=*), parameter :: myname_="ufo_seaicethick_tlad_settraj"

type(ufo_geoval), pointer :: icethick, icefrac, snowthick

! check if sea ice thickness variables is in geovals and get it
call ufo_geovals_get_var(geovals, var_seaicethick, icethick)

! check if sea ice fraction variables is in geovals and get it
call ufo_geovals_get_var(geovals, var_seaicefrac, icefrac)

if (cmp_strings(self%obsvars%variable(1), "seaIceFreeboard")) then
   call ufo_geovals_get_var(geovals, var_seaicesnowthick, snowthick)
   self%snowthick= snowthick
endif

self%icethick = icethick
self%icefrac  = icefrac
self%ltraj    = .true.

end subroutine ufo_seaicethickness_tlad_settraj

! ------------------------------------------------------------------------------
subroutine ufo_seaicethickness_simobs_tl(self, geovals, hofx, obss)
implicit none
class(ufo_seaicethickness_tlad), intent(in)    :: self
type(ufo_geovals),       intent(in)    :: geovals
real(c_double),          intent(inout) :: hofx(:)
type(c_ptr), value,      intent(in)    :: obss

character(len=*), parameter :: myname_="ufo_seaicethick_simobs_tl"
character(max_string) :: err_msg
real(c_double) :: missing

integer :: iobs, icat, ncat
type(ufo_geoval), pointer :: icethick_d, icefrac_d, snowthick
real(kind=kind_real) :: rho_wiw, rho_wsw

!> Set missing value
missing = missing_value(missing)

! check if trajectory was set
if (.not. self%ltraj) then
  write(err_msg,*) myname_, ' trajectory wasnt set!'
  call abor1_ftn(err_msg)
endif

! check if nlocs is consistent in geovals & hofx
if (geovals%nlocs /= size(hofx,1)) then
  write(err_msg,*) myname_, ' error: nlocs inconsistent!'
  call abor1_ftn(err_msg)
endif

! check if sea ice fraction variable is in geovals and get it
call ufo_geovals_get_var(geovals, var_seaicefrac, icefrac_d)

! check if sea ice thickness variable is in geovals and get it
call ufo_geovals_get_var(geovals, var_seaicethick, icethick_d)

if (cmp_strings(self%obsvars%variable(1), "seaIceFreeboard")) then
   rho_wiw = (self%rho_water-self%rho_ice)/self%rho_water
   rho_wsw = (-self%rho_snow)/self%rho_water  
endif

! sea ice thickness obs operator
ncat = icefrac_d%nval
hofx = 0.0

select case (trim(self%obsvars%variable(1)))
case ("seaIceFreeboard")
   do iobs = 1, size(hofx,1)
      do icat = 1, ncat
         ! check for missing input values
         if (self%icefrac%vals(icat,iobs) == missing .or. &
             self%icethick%vals(icat,iobs) == missing .or. &
             self%snowthick%vals(icat,iobs) == missing ) then
            hofx(iobs) = missing
            exit
         end if

         hofx(iobs) = hofx(iobs) +                                         &
                      rho_wiw * self%icefrac%vals(icat,iobs) * icethick_d%vals(icat,iobs) + &
                      rho_wiw * icefrac_d%vals(icat,iobs) * self%icethick%vals(icat,iobs) + &
                      rho_wsw * icefrac_d%vals(icat,iobs) * self%snowthick%vals(icat,iobs)
      enddo
   enddo
case ("iceThickness")
   do iobs = 1, size(hofx,1)
      do icat = 1, ncat
         ! check for missing input values
         if (self%icefrac%vals(icat,iobs) == missing .or. &
             self%icethick%vals(icat,iobs) == missing) then
            hofx(iobs) = missing
            exit
         end if

         hofx(iobs) = hofx(iobs) +                                         &
                      self%icefrac%vals(icat,iobs) * icethick_d%vals(icat,iobs) + &
                      icefrac_d%vals(icat,iobs) * self%icethick%vals(icat,iobs)
      enddo
   enddo
case default
  write(err_msg,*) myname_, ' error: no match seaice thickness_option!'
  call abor1_ftn(err_msg)
end select

end subroutine ufo_seaicethickness_simobs_tl

! ------------------------------------------------------------------------------
subroutine ufo_seaicethickness_simobs_ad(self, geovals, hofx, obss)
implicit none
class(ufo_seaicethickness_tlad), intent(in)    :: self
type(ufo_geovals),       intent(inout) :: geovals
real(c_double),          intent(in)    :: hofx(:)
type(c_ptr), value,      intent(in)    :: obss

character(len=*), parameter :: myname_="ufo_seaicethick_simobs_ad"
character(max_string) :: err_msg

integer :: iobs, icat, ncat
type(ufo_geoval), pointer :: icefrac_d, icethick_d
real(c_double) :: missing
real(kind=kind_real) :: rho_wiw, rho_wsw

!> Set missing value
missing = missing_value(missing)

! check if trajectory was set
if (.not. self%ltraj) then
  write(err_msg,*) myname_, ' trajectory wasnt set!'
  call abor1_ftn(err_msg)
endif

! check if nlocs is consistent in geovals & hofx
if (geovals%nlocs /= size(hofx,1)) then
  write(err_msg,*) myname_, ' error: nlocs inconsistent!'
  call abor1_ftn(err_msg)
endif

if (cmp_strings(self%obsvars%variable(1), "seaIceFreeboard")) then
   rho_wiw = (self%rho_water-self%rho_ice)/self%rho_water
   rho_wsw = (-self%rho_snow)/self%rho_water   
endif

if (.not. geovals%linit ) geovals%linit=.true.

! Get sea-ice fraction & thickness geovals
call ufo_geovals_get_var(geovals, var_seaicefrac, icefrac_d)
call ufo_geovals_get_var(geovals, var_seaicethick, icethick_d)

ncat = self%icethick%nval
if (ncat < 1) then
  write(err_msg,*) myname_, ' unknown number of categories'
  call abor1_ftn(err_msg)
endif

! backward sea ice thickness obs operator

select case (trim(self%obsvars%variable(1)))
case ("seaIceFreeboard")
   do iobs = 1, size(hofx,1)
      if (hofx(iobs) /= missing) then   
         do icat = 1, ncat
            icefrac_d%vals(icat,iobs)  = icefrac_d%vals(icat,iobs)&
                                         + rho_wiw*self%icethick%vals(icat,iobs) * hofx(iobs)&
                                         + rho_wsw*self%snowthick%vals(icat,iobs) * hofx(iobs)
            icethick_d%vals(icat,iobs) = icethick_d%vals(icat,iobs)&
                                         + rho_wiw*self%icefrac%vals(icat,iobs) * hofx(iobs)
         end do
      end if
   enddo
case ("iceThickness")
   do iobs = 1, size(hofx,1)
      if (hofx(iobs) /= missing) then
         do icat = 1, ncat
            icefrac_d%vals(icat,iobs)  = icefrac_d%vals(icat,iobs) + self%icethick%vals(icat,iobs) * hofx(iobs)
            icethick_d%vals(icat,iobs) = icethick_d%vals(icat,iobs) + self%icefrac%vals(icat,iobs) * hofx(iobs)
         end do
      end if
   enddo
case default
  write(err_msg,*) myname_, ' error: no match seaice thickness_option!'
  call abor1_ftn(err_msg)
end select

end subroutine ufo_seaicethickness_simobs_ad

! ------------------------------------------------------------------------------

end module ufo_seaicethickness_tlad_mod
