! (C) Copyright 2021.
!
! This software is developed by NOAA/NWS/EMC under the Apache 2.0 license
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module for insitupm tl/ad observation operator

module ufo_insitupm_tlad_mod

 use iso_c_binding
 use kinds
 use ufo_geovals_mod, only: ufo_geovals, ufo_geoval, ufo_geovals_get_var
 use ufo_vars_mod
 use ufo_constants_mod, only: rd
 use PM_cmaq_mod
 use oops_variables_mod
 use obs_variables_mod

 implicit none
 private

 !> Fortran derived type for the tl/ad observation operator
 type, public :: ufo_insitupm_tlad
 private
  integer :: nlocs, nlayers, nvars, ntracers
  character(len=maxvarlen)     :: model, v_coord
  logical :: scalefactor
  type(obs_variables), public :: obsvars
  type(oops_variables), public :: geovars
  integer(kind=c_int), dimension(:), allocatable :: tracer_modes_cmaq(:)
  real(kind=kind_real), dimension(:,:), allocatable :: prs(:,:)
  real(kind=kind_real), dimension(:,:), allocatable :: ts(:,:)
  real(kind=c_double), dimension(:,:,:), allocatable :: facs(:,:,:) 
  real(kind=c_double), dimension(:), allocatable :: wf(:)
  integer(kind=c_int), dimension(:), allocatable :: wi(:) 
 contains
  procedure :: setup  => ufo_insitupm_tlad_setup_
  procedure :: cleanup  => ufo_insitupm_tlad_cleanup_
  procedure :: settraj => ufo_insitupm_tlad_settraj_
  procedure :: simobs_tl  => ufo_insitupm_simobs_tl_
  procedure :: simobs_ad  => ufo_insitupm_simobs_ad_
  final :: destructor
 end type ufo_insitupm_tlad

contains

! ------------------------------------------------------------------------------
subroutine ufo_insitupm_tlad_setup_(self, f_conf)
use fckit_configuration_module, only: fckit_configuration
implicit none
class(ufo_insitupm_tlad), intent(inout) :: self
type(fckit_configuration), intent(in)  :: f_conf

!Locals
integer :: iq
character(len=:), allocatable :: tracer_variables(:)
character(len=MAXVARLEN) :: err_msg
character(kind=c_char,len=:), allocatable :: model_name, coord_name

  call f_conf%get_or_die("tracer_geovals",tracer_variables)
  self%ntracers = f_conf%get_size("tracer_geovals")

  do iq = 1, self%ntracers
     call self%geovars%push_back(tracer_variables(iq))          ! aerosol species
  enddo

  call f_conf%get_or_die("vertical_coordinate", coord_name)     ! vertical coordinate for interpolation
  self%v_coord = coord_name

  call f_conf%get_or_die("model", model_name)                   ! model name
  self%model = model_name

  ! To be edited (for non-CMAQ models)
  if(self%model .ne. "CMAQ") then
    write(err_msg, *) "not CMAQ, not supported by this operator yet"
    call abor1_ftn(err_msg)
  end if

  ! CMAQ related
  if(self%model .eq. "CMAQ") then

  call f_conf%get_or_die("use_scalefac_cmaq",self%scalefactor)  ! determine whether mode-specific scaling factors will be used

   if(self%scalefactor) then
    if(f_conf%has("tracer_modes_cmaq")) then
    call f_conf%get_or_die("tracer_modes_cmaq",self%tracer_modes_cmaq)  ! CMAQ aerosol modes

    ! check the lengths of aerosol species and aerosol mode lists, stop if inconsistent
     if(f_conf%get_size("tracer_modes_cmaq") .ne. self%ntracers) then
     write(err_msg, *) "mode information missing for some CMAQ aeorsol species"
     call abor1_ftn(err_msg)
     end if
    else
    write(err_msg, *) "CMAQ aerosol mode information not provided"
    call abor1_ftn(err_msg)
    end if

   end if

  end if

  deallocate(tracer_variables, model_name, coord_name)

 ! Size of variables (number of obs type)
  self%nvars = self%obsvars%nvars()

end subroutine ufo_insitupm_tlad_setup_

! ------------------------------------------------------------------------------
subroutine ufo_insitupm_tlad_cleanup_(self)
  implicit none
  class(ufo_insitupm_tlad), intent(inout) :: self

  if (allocated(self%prs))    deallocate(self%prs)
  if (allocated(self%ts))    deallocate(self%ts)
  if (allocated(self%facs))    deallocate(self%facs)
  if (allocated(self%wi))    deallocate(self%wi)
  if (allocated(self%wf))    deallocate(self%wf)

end subroutine ufo_insitupm_tlad_cleanup_

! ------------------------------------------------------------------------------
subroutine ufo_insitupm_tlad_settraj_(self, geovals, obss)
use vert_interp_mod, only: vert_interp_weights
use obsspace_mod

implicit none
class(ufo_insitupm_tlad), intent(inout) :: self
type(ufo_geovals),       intent(in)    :: geovals
type(c_ptr), value,      intent(in)    :: obss

! Local variables
type(ufo_geoval), pointer :: fac1_profile
type(ufo_geoval), pointer :: fac2_profile
type(ufo_geoval), pointer :: fac3_profile
type(ufo_geoval), pointer :: prs_profile
type(ufo_geoval), pointer :: t_profile
type(ufo_geoval), pointer :: hgt_profile
type(ufo_geoval), pointer :: elevation

real(kind_real), dimension(:,:),   allocatable :: hgt   ! height (agl) profiles at obs loc
real(kind_real), dimension(:,:),   allocatable :: hgtasl   ! height (asl) profiles at obs loc
real(kind_real), dimension(:,:),   allocatable :: elev  ! elevation at obs loc

character(len=MAXVARLEN) :: err_msg
real(kind_real), dimension(:), allocatable :: obss_metadata

integer :: ilayer, iloc

 ! Make sure nothing already allocated
 call self%cleanup()

 ! Get number of locations
 self%nlocs = obsspace_get_nlocs(obss)

 ! Get prs and ts from model interp at obs loc (from geovals)
 call ufo_geovals_get_var(geovals, var_prs, prs_profile)
 self%nlayers = prs_profile%nval                ! number of model layers
 allocate(self%prs(self%nlayers,self%nlocs))
 self%prs = prs_profile%vals

 ! Get ts from geovals
 allocate(self%ts(self%nlayers,self%nlocs))
 call ufo_geovals_get_var(geovals, var_ts, t_profile)
 self%ts = t_profile%vals

 ! Allocate arrays for obsspace and interpolation weights
 allocate(obss_metadata(self%nlocs))
 allocate(self%wi(self%nlocs))
 allocate(self%wf(self%nlocs))

 if(self%v_coord .eq. "height_asl") then  
 ! Obs station elevation and model heights (asl) at obs loc
 call obsspace_get_db(obss, "MetaData", "stationElevation", obss_metadata)

 call ufo_geovals_get_var(geovals, var_geomz, hgt_profile)
 allocate(hgt(self%nlayers,self%nlocs))
 hgt = hgt_profile%vals

 allocate(elev(1,self%nlocs))
 call ufo_geovals_get_var(geovals, var_sfc_geomz, elevation)
 elev = elevation%vals

 allocate(hgtasl(self%nlayers,self%nlocs))
 do ilayer = 1, self%nlayers
 hgtasl(ilayer,:) = hgt(ilayer,:) + elev(1,:)
 enddo

 ! Calculate the vertical interpolation weights
 do iloc = 1, self%nlocs
 call vert_interp_weights(self%nlayers, obss_metadata(iloc), hgtasl(:,iloc), &
                          self%wi(iloc), self%wf(iloc))
 end do

 else if(self%v_coord .eq. "log_pressure") then   !log scale
 ! Obs air pressure at obs loc
 call obsspace_get_db(obss, "MetaData", "pressure", obss_metadata)

 ! Calculate the vertical interpolation weights
 do iloc = 1, self%nlocs
 call vert_interp_weights(self%nlayers, log(obss_metadata(iloc)), & 
                          log(self%prs(:,iloc)), self%wi(iloc), self%wf(iloc))
 end do
  
 else
 write(err_msg, *) "coordinate for vertical interpolation not supported"
 call abor1_ftn(err_msg)
 end if

  ! To be edited (for non-CMAQ models)
 if(self%model .ne. "CMAQ") then
   write(err_msg, *) "not CMAQ, not supported by this operator yet"
   call abor1_ftn(err_msg)
 end if

 ! CMAQ related
 if(self%model .eq. "CMAQ") then

  if(self%scalefactor) then
 !Get scaling factors
  allocate(self%facs(3, self%nlayers, self%nlocs))
  call ufo_geovals_get_var(geovals, var_pm25at, fac1_profile)
  self%facs(1,:,:) = fac1_profile%vals

  call ufo_geovals_get_var(geovals, var_pm25ac, fac2_profile)
  self%facs(2,:,:) = fac2_profile%vals

  call ufo_geovals_get_var(geovals, var_pm25co, fac3_profile)
  self%facs(3,:,:) = fac3_profile%vals
  end if

 end if

 deallocate(hgt, hgtasl, elev, obss_metadata)

end subroutine ufo_insitupm_tlad_settraj_

! ------------------------------------------------------------------------------
subroutine ufo_insitupm_simobs_tl_(self, geovals, obss, nvars, nlocs, hofx)

implicit none
class(ufo_insitupm_tlad), intent(in)    :: self
type(ufo_geovals),       intent(in)    :: geovals
type(c_ptr), value,      intent(in)    :: obss
integer,                 intent(in)    :: nvars, nlocs
real(c_double),          intent(inout) :: hofx(nvars,nlocs)

integer :: iq

real(c_double), dimension(:,:,:), allocatable :: qm_tl
type(ufo_geoval), pointer :: aer_profile

character(len=MAXVARLEN) :: geovar
character(len=MAXVARLEN) :: err_msg

 allocate(qm_tl(self%ntracers, self%nlayers, nlocs))
 do iq = 1, self%ntracers
    geovar = self%geovars%variable(iq)                      
    call ufo_geovals_get_var(geovals, geovar, aer_profile)
    qm_tl(iq,:,:) = aer_profile%vals                         
    qm_tl(iq,:,:) = qm_tl(iq,:,:) * self%prs / self%ts / rd          

 enddo

  ! To be edited (for non-CMAQ models)
  if(self%model .ne. "CMAQ") then
    write(err_msg, *) "not CMAQ, not supported by this operator yet"
    call abor1_ftn(err_msg)
  end if

  ! CMAQ related
 if(self%model .eq. "CMAQ") then
  if(self%scalefactor) then
  call get_PM_cmaq_tl(self%nlayers, self%nvars, self%nlocs, self%ntracers,    &
                      self%tracer_modes_cmaq, self%facs,   &
                      qm_tl, self%wi, self%wf, hofx)
  else
  call get_PM_cmaq_tl(self%nlayers, self%nvars, self%nlocs, self%ntracers,    &
                      qm_tl=qm_tl, wi=self%wi, wf=self%wf, pm_tl=hofx)
  end if
 end if

 deallocate(qm_tl)

end subroutine ufo_insitupm_simobs_tl_

! ------------------------------------------------------------------------------
subroutine ufo_insitupm_simobs_ad_(self, geovals, obss, nvars, nlocs, hofx)

implicit none
class(ufo_insitupm_tlad), intent(in)    :: self
type(ufo_geovals),       intent(inout) :: geovals
integer,                 intent(in)    :: nvars, nlocs
real(c_double),          intent(in)    :: hofx(nvars, nlocs)
type(c_ptr), value,      intent(in)    :: obss

integer :: iq

real(c_double), dimension(:,:,:), allocatable :: qm_ad

type(ufo_geoval), pointer :: aer_profile
character(len=MAXVARLEN) :: geovar
character(len=MAXVARLEN) :: err_msg

 allocate(qm_ad(self%ntracers, self%nlayers, nlocs)) 

  ! To be edited (for non-CMAQ models)
  if(self%model .ne. "CMAQ") then
    write(err_msg, *) "not CMAQ, not supported by this operator yet"
    call abor1_ftn(err_msg)
  end if

  ! CMAQ related
 if(self%model .eq. "CMAQ") then
  if(self%scalefactor) then  

  call get_PM_cmaq_ad(self%nlayers, self%nvars, self%nlocs, self%ntracers, &
                      self%tracer_modes_cmaq, &
                      self%facs, self%wi, self%wf, hofx, qm_ad)
  else
  call get_PM_cmaq_ad(self%nlayers, self%nvars, self%nlocs, self%ntracers,    &
                      wi=self%wi, wf=self%wf, pm_ad=hofx, qm_ad=qm_ad)
  end if
 end if

 do iq = self%ntracers,1,-1

   geovar = self%geovars%variable(iq)                   
   call ufo_geovals_get_var(geovals, geovar, aer_profile)
   if (.not. allocated(aer_profile%vals)) then
       aer_profile%nprofiles = self%nlocs
       aer_profile%nval  = self%nlayers
       allocate(aer_profile%vals(aer_profile%nval, aer_profile%nprofiles))
       aer_profile%vals(:,:) = 0.0_kind_real
   endif

   qm_ad(iq,:,:) = qm_ad(iq,:,:) * self%prs / self%ts / rd
   aer_profile%vals = qm_ad(iq,:,:)

 enddo

 deallocate(qm_ad)

end subroutine ufo_insitupm_simobs_ad_

! ------------------------------------------------------------------------------

subroutine  destructor(self)
  type(ufo_insitupm_tlad), intent(inout)  :: self

  call self%cleanup()

  if (allocated(self%tracer_modes_cmaq)) deallocate(self%tracer_modes_cmaq)

end subroutine destructor

end module ufo_insitupm_tlad_mod
