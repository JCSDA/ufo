! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module to handle tl/ad for aod observations

module ufo_aodcrtm_tlad_mod

 use fckit_configuration_module, only: fckit_configuration
 use iso_c_binding
 use kinds
 use missing_values_mod

 use ufo_geovals_mod, only: ufo_geovals, ufo_geoval, ufo_geovals_get_var
 use ufo_vars_mod
 use ufo_crtm_utils_mod
 use crtm_module
 use obsspace_mod

 implicit none
 private

 !> Fortran derived type for aod trajectory
 type, public :: ufo_aodcrtm_tlad
 private
  character(len=MAXVARLEN), public, allocatable :: varin(:)  ! variablesrequested from the model
  integer, allocatable                          :: channels(:)
  type(crtm_conf) :: conf
  integer :: n_Profiles
  integer :: n_Layers
  integer :: n_Channels
  type(CRTM_Atmosphere_type), allocatable :: atm_K(:,:)
  type(CRTM_Surface_type), allocatable :: sfc_K(:,:)
  REAL(kind_real), allocatable  :: scaling_factor(:,:)  
  logical :: ltraj
 contains
  procedure :: setup  => ufo_aodcrtm_tlad_setup
  procedure :: delete  => ufo_aodcrtm_tlad_delete
  procedure :: settraj => ufo_aodcrtm_tlad_settraj
  procedure :: simobs_tl  => ufo_aodcrtm_simobs_tl
  procedure :: simobs_ad  => ufo_aodcrtm_simobs_ad
 end type ufo_aodcrtm_tlad

contains

! ------------------------------------------------------------------------------

subroutine ufo_aodcrtm_tlad_setup(self, f_confOper, channels)

implicit none
class(ufo_aodcrtm_tlad),   intent(inout) :: self
type(fckit_configuration), intent(in)    :: f_confOper
integer(c_int),               intent(in)    :: channels(:)  !List of channels to use

type(fckit_configuration) :: f_confOpts
integer :: nvars_in

CHARACTER(len=MAXVARLEN), ALLOCATABLE :: var_aerosols(:)

 call f_confOper%get_or_die("obs options",f_confOpts)

 call crtm_conf_setup(self%conf, f_confOpts, f_confOper)

 CALL assign_aerosol_names(self%conf%aerosol_option,var_aerosols)

 nvars_in = size(var_aerosols)
 allocate(self%varin(nvars_in))
 self%varin(1:nvars_in) = var_aerosols

 ! save channels
 allocate(self%channels(size(channels)))
 self%channels(:) = channels(:)

end subroutine ufo_aodcrtm_tlad_setup

! ------------------------------------------------------------------------------

subroutine ufo_aodcrtm_tlad_delete(self)

implicit none
class(ufo_aodcrtm_tlad), intent(inout) :: self

 self%ltraj = .false.
 call crtm_conf_delete(self%conf)

 if (allocated(self%atm_k)) then
   call CRTM_Atmosphere_Destroy(self%atm_K)
   deallocate(self%atm_k)
 endif

 if (allocated(self%sfc_k)) then
   call CRTM_Surface_Destroy(self%sfc_K)
   deallocate(self%sfc_k)
 endif

 IF (ALLOCATED(self%scaling_factor)) THEN
   deallocate(self%scaling_factor)
 endif


end subroutine ufo_aodcrtm_tlad_delete

! ------------------------------------------------------------------------------

SUBROUTINE ufo_aodcrtm_tlad_settraj(self, geovals, obss)

implicit none

class(ufo_aodcrtm_tlad), intent(inout) :: self
type(ufo_geovals),        intent(in)    :: geovals
type(c_ptr), value,       intent(in)    :: obss

! Local Variables
character(*), parameter :: PROGRAM_NAME = 'ufo_aodcrtm_tlad_mod.F90'
character(len=MAXVARLEN) :: def_aero_mod
character(255) :: message, version
integer        :: err_stat, alloc_stat
INTEGER        :: n,l,m
type(ufo_geoval), pointer :: temp

! Define the "non-demoninational" arguments
type(CRTM_ChannelInfo_type)             :: chinfo(self%conf%n_Sensors)
type(CRTM_Geometry_type),   allocatable :: geo(:)

! Define the FORWARD variables
type(CRTM_Atmosphere_type), allocatable :: atm(:)
type(CRTM_Surface_type),    allocatable :: sfc(:)
type(CRTM_RTSolution_type), allocatable :: rts(:,:)

! Define the K-MATRIX variables
type(CRTM_RTSolution_type), allocatable :: rts_K(:,:)


 ! Get number of profile and layers from geovals
 ! ---------------------------------------------
 self%n_Profiles = geovals%nlocs
 call ufo_geovals_get_var(geovals, var_ts, temp)
 self%n_Layers = temp%nval
 nullify(temp)

 ! Program header
 ! --------------
 ! call CRTM_Version( Version )
 ! call Program_Message( PROGRAM_NAME, &
 !                      'Check/example program for the CRTM Forward and K-Matrix (setTraj) functions using '//&
 !                      trim(self%conf%ENDIAN_type)//' coefficient datafiles', &
 !                      'CRTM Version: '//TRIM(Version) )


 ! Initialise all the sensors at once
 ! ----------------------------------
 !** NOTE: CRTM_Init points to the various binary files needed for CRTM.  See the
 !**       CRTM_Lifecycle.f90 for more details.

 ! write( *,'(/5x,"Initializing the CRTM (setTraj) ...")' )
 call define_aerosol_model(self%conf%AerosolCoeff_File, def_aero_mod)
 err_stat = CRTM_Init( self%conf%SENSOR_ID, &
            chinfo, &
            File_Path           = trim(self%conf%COEFFICIENT_PATH), &
            NC_File_Path        = trim(self%conf%NC_COEFFICIENT_PATH), &
            Aerosol_Model       = trim(def_aero_mod), &
            AerosolCoeff_Format = trim(self%conf%AerosolCoeff_Format), &
            AerosolCoeff_File   = trim(self%conf%AerosolCoeff_File), &
            Quiet=.TRUE.)
 if ( err_stat /= SUCCESS ) THEN
   message = 'Error initializing CRTM (setTraj)'
   call Display_Message( PROGRAM_NAME, message, FAILURE )
   stop
 end if

 ! Loop over all sensors. Not necessary if we're calling CRTM for each sensor
 ! ----------------------------------------------------------------------------
 Sensor_Loop:do n = 1, self%conf%n_Sensors


   ! Determine the number of channels for the current sensor
   ! -------------------------------------------------------
   self%N_Channels = CRTM_ChannelInfo_n_Channels(chinfo(n))


   ! Allocate the ARRAYS
   ! -------------------
   allocate( geo( self%n_Profiles )                         , &
             atm( self%n_Profiles )                         , &
             sfc( self%n_Profiles )                         , &
             rts( self%N_Channels, self%n_Profiles )        , &
             self%atm_K( self%N_Channels, self%n_Profiles ) , &
             self%sfc_K( self%N_Channels, self%n_Profiles ) , &
             self%scaling_factor(self%n_Layers,self%n_Profiles ) , &
             rts_K( self%N_Channels, self%n_Profiles )      , &
             STAT = alloc_stat                                )
   if ( alloc_stat /= 0 ) THEN
      message = 'Error allocating structure arrays (setTraj)'
      call Display_Message( PROGRAM_NAME, message, FAILURE )
      stop
   end if


   ! Create the input FORWARD structure (atm)
   ! ----------------------------------------
   call CRTM_Atmosphere_Create( atm, self%n_Layers, self%conf%n_Absorbers, self%conf%n_Clouds, self%conf%n_Aerosols )
   if ( ANY(.NOT. CRTM_Atmosphere_Associated(atm)) ) THEN
      message = 'Error allocating CRTM Forward Atmosphere structure (setTraj)'
      CALL Display_Message( PROGRAM_NAME, message, FAILURE )
      STOP
   END IF

   ! Create output K-MATRIX structure (atm)
   ! --------------------------------------
   call CRTM_Atmosphere_Create( self%atm_K, self%n_Layers, self%conf%n_Absorbers, self%conf%n_Clouds, self%conf%n_Aerosols )
   if ( ANY(.NOT. CRTM_Atmosphere_Associated(self%atm_K)) ) THEN
      message = 'Error allocating CRTM K-matrix Atmosphere structure (setTraj)'
      CALL Display_Message( PROGRAM_NAME, message, FAILURE )
      STOP
   END IF

   !Assign the data from the GeoVaLs
   !--------------------------------
   CALL Load_Atm_Data(self%N_PROFILES,self%N_LAYERS,geovals,atm,self%conf)

   IF (TRIM(self%conf%aerosol_option) /= "") &
        &CALL load_aerosol_data(self%n_profiles, self%n_layers, geovals,&
        &self%conf, self%varin, trim(def_aero_mod), atm)

   CALL CRTM_RTSolution_Create(rts, self%n_layers )
   CALL CRTM_RTSolution_Create(rts_k, self%n_layers )

   ! Zero the K-matrix OUTPUT structures
   ! -----------------------------------
   call CRTM_Atmosphere_Zero( self%atm_K )

   CALL calculate_aero_layer_factor(atm,self%scaling_factor)

   ! Inintialize the K-matrix INPUT so that the results are daero/dx
   ! -------------------------------------------------------------

   FORALL (l=1:self%n_channels,m=1:self%n_profiles) rts_k(l,m)%layer_optical_depth = one

   ! Call the K-matrix model
   ! -----------------------
   err_stat = CRTM_AOD_K( atm         , &  ! FORWARD  Input
                             rts_K       , &  ! K-MATRIX Input
                             chinfo(n:n) , &  ! Input
                             rts         , &  ! FORWARD  Output
                             self%atm_K    )  ! K-MATRIX Output
   if ( err_stat /= SUCCESS ) THEN
      message = 'Error calling CRTM (setTraj) K-Matrix Model for '//TRIM(self%conf%SENSOR_ID(n))
      call Display_Message( PROGRAM_NAME, message, FAILURE )
      stop
   end if

   ! Deallocate the structures
   ! -------------------------
   call CRTM_Atmosphere_Destroy(atm)
   call CRTM_RTSolution_Destroy(rts_K)
   call CRTM_RTSolution_Destroy(rts)



   ! Deallocate all arrays
   ! ---------------------
   deallocate(geo, atm, sfc, rts, rts_K, STAT = alloc_stat)
   if ( alloc_stat /= 0 ) THEN
      message = 'Error deallocating structure arrays (setTraj)'
      call Display_Message( PROGRAM_NAME, message, FAILURE )
      stop
   end if

 end do Sensor_Loop


 ! Destroy CRTM instance
 ! ---------------------
 ! write( *, '( /5x, "Destroying the CRTM (setTraj)..." )' )
 err_stat = CRTM_Destroy( chinfo )
 if ( err_stat /= SUCCESS ) THEN
    message = 'Error destroying CRTM (setTraj)'
    call Display_Message( PROGRAM_NAME, message, FAILURE )
    stop
 end if


 ! Set flag that the tracectory was set
 ! ------------------------------------
 self%ltraj = .true.

end subroutine ufo_aodcrtm_tlad_settraj

! ------------------------------------------------------------------------------

SUBROUTINE ufo_aodcrtm_simobs_tl(self, geovals, obss, nvars, nlocs, hofx)

implicit none
class(ufo_aodcrtm_tlad), intent(in) :: self
type(ufo_geovals),        intent(in) :: geovals
type(c_ptr), value,    intent(in)    :: obss
integer,                  intent(in)    :: nvars, nlocs
real(c_double),        intent(inout) :: hofx(nvars, nlocs)

character(len=*), parameter :: myname_="ufo_aodcrtm_simobs_tl"
character(max_string) :: err_msg
integer :: jprofile, jchannel, jlevel, jaero
type(ufo_geoval), pointer :: var_p

CHARACTER(len=MAXVARLEN), ALLOCATABLE :: var_aerosols(:)


 ! Initial checks
 ! --------------

 ! Check if trajectory was set
 if (.not. self%ltraj) then
   write(err_msg,*) myname_, ' trajectory wasnt set!'
   call abor1_ftn(err_msg)
 endif

 ! Check if nlocs is consistent in geovals & hofx
 if (geovals%nlocs /= self%n_Profiles) then
   write(err_msg,*) myname_, ' error: nlocs inconsistent!'
   call abor1_ftn(err_msg)
 endif

 CALL assign_aerosol_names(self%conf%aerosol_option,var_aerosols)

 IF (SIZE(var_aerosols) /= self%conf%n_aerosols) THEN
    WRITE(err_msg,*) myname_, ' error: n_aerosols inconsistent!'
    call abor1_ftn(err_msg)
 ENDIF


 call ufo_geovals_get_var(geovals, var_aerosols(1), var_p)

 ! Check model levels is consistent in geovals & crtm
 if (var_p%nval /= self%n_Layers) then
   write(err_msg,*) myname_, ' error: layers inconsistent!'
   call abor1_ftn(err_msg)
 endif


 ! Initialize hofx
 ! ---------------
 hofx(:,:) = 0.0_kind_real

 ! Multiply by Jacobian and add to hofx
 do jprofile = 1, self%n_profiles

   do jchannel = 1, size(self%channels)
     DO jaero = 1, self%conf%n_aerosols
        CALL ufo_geovals_get_var(geovals, var_aerosols(jaero), var_p)
        DO jlevel = 1, var_p%nval
           hofx(jchannel, jprofile) = hofx(jchannel, jprofile) + &
                                      self%atm_k(self%channels(jchannel),jprofile)%aerosol(jaero)%concentration(jlevel) *  &
                                      var_p%vals(jlevel,jprofile) * self%scaling_factor(jlevel,jprofile) * self%conf%unit_coef
        ENDDO
     ENDDO
   enddo
 enddo

end subroutine ufo_aodcrtm_simobs_tl

! ------------------------------------------------------------------------------

SUBROUTINE ufo_aodcrtm_simobs_ad(self, geovals, obss, nvars, nlocs, hofx)

implicit none
class(ufo_aodcrtm_tlad), intent(in) :: self
type(ufo_geovals),     intent(inout) :: geovals
type(c_ptr), value,    intent(in)    :: obss
integer,                  intent(in)    :: nvars, nlocs
real(c_double),           intent(in) :: hofx(nvars, nlocs)

character(len=*), parameter :: myname_="ufo_aodcrtm_simobs_ad"
character(max_string) :: err_msg
integer :: jprofile, jchannel, jlevel
type(ufo_geoval), pointer :: var_p
real(c_double) :: missing

CHARACTER(len=MAXVARLEN), ALLOCATABLE :: var_aerosols(:)

INTEGER :: jaero

 ! Initial checks
 ! --------------

 ! Check if trajectory was set
 if (.not. self%ltraj) then
   write(err_msg,*) myname_, ' trajectory wasnt set!'
   call abor1_ftn(err_msg)
 endif

 ! Check if nlocs is consistent in geovals & hofx
 if (geovals%nlocs /= self%n_Profiles) then
   write(err_msg,*) myname_, ' error: nlocs inconsistent!'
   call abor1_ftn(err_msg)
 endif

 ! Set missing value
 missing = missing_value(missing)

 CALL assign_aerosol_names(self%conf%aerosol_option,var_aerosols)

 DO jaero=1,self%conf%n_aerosols

    CALL ufo_geovals_get_var(geovals, var_aerosols(jaero), var_p)

! Multiply by Jacobian and add to hofx (adjoint)
    DO jprofile = 1, self%n_Profiles
       DO jchannel = 1, size(self%channels)
          if (hofx(jchannel, jprofile) /= missing) then
            DO jlevel = 1, var_p%nval
               var_p%vals(jlevel,jprofile) = var_p%vals(jlevel,jprofile) + &
                    self%atm_k(self%channels(jchannel),jprofile)%aerosol(jaero)%concentration(jlevel) * hofx(jchannel, jprofile)
            ENDDO
          end if
       ENDDO
    ENDDO

    FORALL (jlevel=1:var_p%nval,jprofile=1:self%n_profiles) &
        var_p%vals(jlevel,jprofile) = var_p%vals(jlevel,jprofile) * self%scaling_factor(jlevel,jprofile) * self%conf%unit_coef

 ENDDO

 ! Once all geovals set replace flag
 ! ---------------------------------
 if (.not. geovals%linit ) geovals%linit=.true.


end subroutine ufo_aodcrtm_simobs_ad

! ------------------------------------------------------------------------------

END MODULE ufo_aodcrtm_tlad_mod
