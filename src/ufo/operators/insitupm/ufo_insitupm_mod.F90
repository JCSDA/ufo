! (C) Copyright 2021.
!
! This software is developed by NOAA/NWS/EMC under the Apache 2.0 license
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module for insitupm observation operator

module ufo_insitupm_mod

 use iso_c_binding
 use kinds
 use PM_cmaq_mod
 use oops_variables_mod
 use obs_variables_mod
 use ufo_vars_mod
 use ufo_basis_mod, only: ufo_basis

 implicit none
 private

!> Fortran derived type for the observation type
 type, public :: ufo_insitupm
 private
   type(oops_variables), public :: geovars
   type(obs_variables), public :: obsvars
   integer                      :: ntracers
   character(len=maxvarlen)     :: model, v_coord        
   logical                      :: scalefactor
   integer, allocatable :: tracer_modes_cmaq(:) 
 contains
   procedure :: setup  => ufo_insitupm_setup
   procedure :: simobs => ufo_insitupm_simobs
   final :: destructor
 end type ufo_insitupm

!> Default variables required from model: air pressure, temperature, CMAQ scaling factors, height and surface height
 character(len=maxvarlen), dimension(2), parameter :: varindefault = (/var_prs, var_ts/)
 character(len=maxvarlen), dimension(3), parameter :: varincmaq = (/var_pm25at, var_pm25ac, var_pm25co/)
 character(len=maxvarlen), dimension(2), parameter :: varvertical = (/var_geomz, var_sfc_geomz/)
contains

! ------------------------------------------------------------------------------
subroutine ufo_insitupm_setup(self, f_conf)
use fckit_configuration_module, only: fckit_configuration
implicit none
class(ufo_insitupm), intent(inout)     :: self
type(fckit_configuration), intent(in) :: f_conf

!Locals
integer :: iq, nvars
character(len=:), allocatable :: tracer_variables(:)
character(len=MAXVARLEN) :: err_msg
character(kind=c_char,len=:), allocatable :: model_name, coord_name

  ! Fill in geovars: variables we need from the model
  ! Let users choose specific aerosol species (defined in the yaml file) needed in the PM2.5/total PM calculation.
  ! Followed by slots for prs, ts, heights, and mode-specific scaling factors. 

  call f_conf%get_or_die("tracer_geovals",tracer_variables)
  self%ntracers = f_conf%get_size("tracer_geovals")

  do iq = 1, self%ntracers
     call self%geovars%push_back(tracer_variables(iq))         ! aerosol species
  enddo

  ! Size of variables (number of obs type (= 1 for this case))
  nvars = self%obsvars%nvars()

  call self%geovars%push_back(varindefault)                    ! pressure and ts (needed for unit conversion to ug/m3) 

  call f_conf%get_or_die("vertical_coordinate", coord_name)    ! vertical coordinate for interpolation
  self%v_coord = coord_name

  if(self%v_coord .eq. "height_asl") then
  call self%geovars%push_back(varvertical)                     ! height and surface height
  endif

  call f_conf%get_or_die("model", model_name)                  ! model name 
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

    call self%geovars%push_back(varincmaq)                    ! CMAQ scaling factors for three modes
   end if   

  end if  
   
  deallocate(tracer_variables, model_name, coord_name)

end subroutine ufo_insitupm_setup

! ------------------------------------------------------------------------------
subroutine destructor(self)
implicit none
type(ufo_insitupm), intent(inout) :: self

if (allocated(self%tracer_modes_cmaq)) deallocate(self%tracer_modes_cmaq)

end subroutine destructor

! ------------------------------------------------------------------------------
subroutine ufo_insitupm_simobs(self, geovals, obss, nvars, nlocs, hofx)
use ufo_constants_mod, only: rd
use vert_interp_mod, only: vert_interp_weights
use ufo_geovals_mod
use obsspace_mod

implicit none
class(ufo_insitupm), intent(in)   :: self
integer, intent(in)               :: nvars, nlocs
type(ufo_geovals),  intent(in)    :: geovals
real(c_double),     intent(inout) :: hofx(nvars, nlocs)
type(c_ptr), value, intent(in)    :: obss

! Local variables
type(ufo_geoval), pointer :: aer_profile
type(ufo_geoval), pointer :: fac1_profile
type(ufo_geoval), pointer :: fac2_profile
type(ufo_geoval), pointer :: fac3_profile
type(ufo_geoval), pointer :: prs_profile
type(ufo_geoval), pointer :: t_profile
type(ufo_geoval), pointer :: hgt_profile
type(ufo_geoval), pointer :: elevation

integer :: nlayers, iq, ilayer, iloc

real(c_double), dimension(:,:,:), allocatable :: qm    ! aerosol mass mix ratio(ug/kg) that becomes concentration (*prs/t/rd) profiles at obs loc
real(kind_real), dimension(:,:),  allocatable :: ts    ! temperature profiles at obs loc
real(kind_real), dimension(:,:),  allocatable :: prs   ! air pressure profiles at obs loc
real(c_double), dimension(:,:,:), allocatable :: facs  ! aerosol scaling factor profiles pm25at, pm25ac, pm25co

real(kind_real), dimension(:,:),   allocatable :: hgt   ! height (agl) profiles at obs loc
real(kind_real), dimension(:,:),   allocatable :: hgtasl   ! height (asl) profiles at obs loc
real(kind_real), dimension(:,:),   allocatable :: elev  ! elevation at obs loc

real(kind_real), dimension(:), allocatable :: obss_metadata
real(c_double), allocatable :: wf(:)
integer(c_int), allocatable :: wi(:)

character(len=MAXVARLEN) :: geovar
character(len=MAXVARLEN) :: err_msg
character(len=MAXVARLEN), dimension(:), allocatable:: tracer_name


  ! Get prs and ts from model interp at obs loc (from geovals)
  call ufo_geovals_get_var(geovals, var_prs, prs_profile)
  nlayers = prs_profile%nval                            ! number of model layers

  allocate(prs(nlayers,nlocs))
  prs = prs_profile%vals

  ! Get ts from geovals
  allocate(ts(nlayers,nlocs))
  call ufo_geovals_get_var(geovals, var_ts, t_profile)
  ts = t_profile%vals

  ! Allocate arrays for obsspace and interpolation weights
  allocate(obss_metadata(nlocs))
  allocate(wi(nlocs))
  allocate(wf(nlocs))

  if(self%v_coord .eq. "height_asl") then  
  ! Obs station elevation and model heights (asl) at obs loc
  call obsspace_get_db(obss, "MetaData", "stationElevation", obss_metadata)

  call ufo_geovals_get_var(geovals, var_geomz, hgt_profile)
  allocate(hgt(nlayers,nlocs))
  hgt = hgt_profile%vals

  allocate(elev(1,nlocs))
  call ufo_geovals_get_var(geovals, var_sfc_geomz, elevation)
  elev = elevation%vals

  allocate(hgtasl(nlayers,nlocs))
  do ilayer = 1, nlayers
  hgtasl(ilayer,:) = hgt(ilayer,:) + elev(1,:)
  enddo

  ! Calculate the vertical interpolation weights
  do iloc = 1, nlocs
  call vert_interp_weights(nlayers, obss_metadata(iloc), &
                           hgtasl(:,iloc), wi(iloc), wf(iloc))
  end do

  else if(self%v_coord .eq. "log_pressure") then  !log scale
  ! Obs air pressure at obs loc
  call obsspace_get_db(obss, "MetaData", "pressure", obss_metadata)

  ! Calculate the vertical interpolation weights
  do iloc = 1, nlocs 
  call vert_interp_weights(nlayers, log(obss_metadata(iloc)), & 
                           log(prs(:,iloc)), wi(iloc), wf(iloc))
  end do
  
  else
  write(err_msg, *) "coordinate for vertical interpolation not supported"
  call abor1_ftn(err_msg)
  end if

  ! Get aerosol profiles interpolated at obs loc
  allocate(qm(self%ntracers, nlayers, nlocs))
  allocate(tracer_name(self%ntracers))
  do iq = 1, self%ntracers
     geovar = self%geovars%variable(iq)                   !self%geovars contains tracers 
     tracer_name(iq) = geovar
     call ufo_geovals_get_var(geovals, geovar, aer_profile)
     qm(iq,:,:) = aer_profile%vals             ! aerosol mass mixing ratio
     qm(iq,:,:) = qm(iq,:,:) * prs / ts / rd   ! aerosol concentration (ug/m3) 
  enddo

  ! To be edited (for non-CMAQ models)
  if(self%model .ne. "CMAQ") then
    write(err_msg, *) "not CMAQ, not supported by this operator yet"
    call abor1_ftn(err_msg)
  end if

  ! CMAQ related
  if(self%model .eq. "CMAQ") then

  hofx = 0.0

   if(self%scalefactor) then
  ! Get scaling factors from geovals
   allocate(facs(3, nlayers, nlocs))
   call ufo_geovals_get_var(geovals, var_pm25at, fac1_profile)
   facs(1,:,:) = fac1_profile%vals           
   call ufo_geovals_get_var(geovals, var_pm25ac, fac2_profile)
   facs(2,:,:) = fac2_profile%vals  
   call ufo_geovals_get_var(geovals, var_pm25co, fac3_profile)
   facs(3,:,:) = fac3_profile%vals  
 
   call get_PM_cmaq(nlayers, nvars, nlocs, self%ntracers, &
                    self%tracer_modes_cmaq, &
                    facs, qm, wi, wf, hofx)
   else
   call get_PM_cmaq(nlayers, nvars, nlocs, self%ntracers, & 
                    qm=qm, wi=wi, wf=wf, pm=hofx)
   end if

  end if
 
  ! cleanup memory
  ! --------
  deallocate(qm, ts, prs, tracer_name, wi, wf, obss_metadata)
  if (allocated(facs))    deallocate(facs)
  if (allocated(hgt))    deallocate(hgt)
  if (allocated(hgtasl))    deallocate(hgtasl)
  if (allocated(elev))    deallocate(elev)

end subroutine ufo_insitupm_simobs


! ------------------------------------------------------------------------------

end module ufo_insitupm_mod
