! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module to handle tl/ad for radiancecrtm observations

module ufo_radiancecrtm_tlad_mod
 use,intrinsic :: iso_c_binding
 use,intrinsic :: iso_fortran_env
 use crtm_module

 use fckit_configuration_module, only: fckit_configuration
 use fckit_mpi_module,   only: fckit_mpi_comm
 use iso_c_binding
 use kinds
 use missing_values_mod

 use obsspace_mod
 use obsdatavector_mod
 use ufo_geovals_mod, only: ufo_geovals, ufo_geoval, ufo_geovals_get_var
 use ufo_vars_mod
 use ufo_crtm_utils_mod

 use ufo_constants_mod, only: deg2rad
 use ufo_crtm_passive_mod
 use ufo_crtm_active_mod

 implicit none
 private

 !> Fortran derived type for radiancecrtm trajectory
 type, public :: ufo_radiancecrtm_tlad
 private
  character(len=MAXVARLEN), public, allocatable :: varin(:)  ! variables requested from the model
  integer, allocatable                          :: channels(:)
  type(crtm_conf) :: conf
  type(crtm_conf) :: conf_traj
  integer :: n_Profiles
  integer :: n_Layers
  integer :: n_Channels
  type(CRTM_Atmosphere_type), allocatable :: atm_K(:,:)
  type(CRTM_Surface_type), allocatable :: sfc_K(:,:)
  logical :: ltraj
  type(CRTM_Options_type), allocatable :: Options(:)
  logical :: use_qc_flags
 contains
  procedure :: setup  => ufo_radiancecrtm_tlad_setup
  procedure :: delete  => ufo_radiancecrtm_tlad_delete
  procedure :: settraj => ufo_radiancecrtm_tlad_settraj
  procedure :: simobs_tl  => ufo_radiancecrtm_simobs_tl
  procedure :: simobs_ad  => ufo_radiancecrtm_simobs_ad
 end type ufo_radiancecrtm_tlad

 character(len=maxvarlen), dimension(1), parameter :: varin_default = &
                            (/var_ts/)

contains

! ------------------------------------------------------------------------------

subroutine ufo_radiancecrtm_tlad_setup(self, f_confOper, channels, comm)

implicit none
class(ufo_radiancecrtm_tlad), intent(inout) :: self
type(fckit_configuration),    intent(in)    :: f_confOper
integer(c_int),               intent(in)    :: channels(:)  !List of channels to use
type(fckit_mpi_comm),         intent(in)    :: comm

integer :: nvars_in
integer :: ind, jspec
type(fckit_configuration) :: f_confOpts,f_confLinOper
character(max_string) :: err_msg

 call f_confOper%get_or_die("obs options",f_confOpts)
 call crtm_conf_setup(self%conf_traj, f_confOpts, f_confOper, comm)
 call f_confOper%get_or_die("UseQCFlagsToSkipHofX",self%use_qc_flags)  

 if ( f_confOper%has("linear obs operator") ) then
    call f_confOper%get_or_die("linear obs operator",f_confLinOper)
    call crtm_conf_setup(self%conf, f_confOpts, f_confLinOper, comm)
 else
    call crtm_conf_setup(self%conf, f_confOpts, f_confOper, comm)
 end if

 ! request from the model var_ts +
 ! 1 * n_Absorbers
 ! 1 * n_Clouds (mass content only)
 nvars_in = size(varin_default) + self%conf%n_Absorbers + self%conf%n_Clouds + self%conf%n_Surfaces
 allocate(self%varin(nvars_in))
 self%varin(1:size(varin_default)) = varin_default
 ind = size(varin_default) + 1

 !Use list of Absorbers and Clouds from conf
 do jspec = 1, self%conf%n_Absorbers
   self%varin(ind) = self%conf%Absorbers(jspec)
   ind = ind + 1
 end do
 do jspec = 1, self%conf%n_Clouds
   self%varin(ind) = self%conf%Clouds(jspec,1)
   ind = ind + 1
 end do
 do jspec = 1, self%conf%n_Surfaces
   self%varin(ind) = self%conf%Surfaces(jspec)
   ind = ind + 1
 end do
 if ( (ufo_vars_getindex(self%varin, var_sfc_wspeed) > 0 .or. &
       ufo_vars_getindex(self%varin, var_sfc_wdir) > 0) .and. &
      trim(self%conf_traj%sfc_wind_geovars) /= "vector") then
   write(err_msg,*) 'ufo_radiancecrtm_tlad_setup error: sfc_wind_geovars not supported in tlad --> ', self%conf_traj%sfc_wind_geovars
   call abor1_ftn(err_msg)
 end if

 ! save channels
 allocate(self%channels(size(channels)))
 self%channels(:) = channels(:)

end subroutine ufo_radiancecrtm_tlad_setup

! ------------------------------------------------------------------------------

subroutine ufo_radiancecrtm_tlad_delete(self)

implicit none
class(ufo_radiancecrtm_tlad), intent(inout) :: self

 self%ltraj = .false.

 call crtm_conf_delete(self%conf)
 call crtm_conf_delete(self%conf_traj)

 if (allocated(self%atm_k)) then
   call CRTM_Atmosphere_Destroy(self%atm_K)
   deallocate(self%atm_k)
 endif

 if (allocated(self%sfc_k)) then
   call CRTM_Surface_Destroy(self%sfc_K)
   deallocate(self%sfc_k)
 endif

 if (allocated(self%Options)) then
    call CRTM_Options_Destroy(self%Options)
    deallocate(self%Options)
 endif

end subroutine ufo_radiancecrtm_tlad_delete

! ------------------------------------------------------------------------------

subroutine ufo_radiancecrtm_tlad_settraj(self, geovals, obss, hofxdiags)
use fckit_mpi_module,   only: fckit_mpi_comm
use fckit_log_module,   only: fckit_log
use ieee_arithmetic,    only: ieee_is_nan
use ufo_utils_mod,      only: cmp_strings
use CRTM_SpcCoeff, only: SC, &
                         SpcCoeff_IsMicrowaveSensor , &
                         SpcCoeff_IsInfraredSensor  , &
                         SpcCoeff_IsVisibleSensor   , &
                         SpcCoeff_IsUltravioletSensor

implicit none

class(ufo_radiancecrtm_tlad), intent(inout) :: self
type(ufo_geovals),        intent(in)    :: geovals
type(c_ptr), value,       intent(in)    :: obss
type(ufo_geovals),        intent(inout) :: hofxdiags    !non-h(x) diagnostics

! Local Variables
character(*), parameter :: PROGRAM_NAME = 'ufo_radiancecrtm_tlad_settraj'
character(255) :: message, version
character(max_string) :: err_msg
integer        :: err_stat, alloc_stat
integer        :: n
type(ufo_geoval), pointer :: temp
integer :: jvar, ichannel, jchannel, jprofile, jlevel, jspec
real(c_double) :: missing
type(fckit_mpi_comm)  :: f_comm

! Define the "non-demoninational" arguments
type(CRTM_ChannelInfo_type)             :: chinfo(self%conf_traj%n_Sensors)
type(CRTM_Geometry_type),   allocatable :: geo(:)

! Define the FORWARD variables
type(CRTM_Atmosphere_type), allocatable :: atm(:)
type(CRTM_Surface_type),    allocatable :: sfc(:)
type(CRTM_RTSolution_type), allocatable :: rts(:,:)

! Define the K-MATRIX variables
type(CRTM_RTSolution_type), allocatable :: rts_K(:,:)

!for gmi
type(CRTM_Geometry_type),   allocatable :: geo_hf(:)
type(CRTM_Atmosphere_type), allocatable :: atm_Ka(:,:)
type(CRTM_Surface_type),    allocatable :: sfc_Ka(:,:)
type(CRTM_RTSolution_type), allocatable :: rtsa(:,:)
type(CRTM_RTSolution_type), allocatable :: rts_Ka(:,:)

character(len=MAXVARLEN),  dimension(hofxdiags%nvar) :: &
                          ystr_diags, xstr_diags
integer :: str_pos(4), ch_diags(hofxdiags%nvar),numNaN
integer :: lch
! Used to parse hofxdiags
character(len=MAXVARLEN) :: varstr
character(10), parameter :: jacobianstr = "_jacobian_"
integer(c_size_t) :: nvars, nlocs
logical :: qc_ff

! set a local boolean variable for whether we are in vis or ultraviolet channels
logical :: Is_Vis_or_UV = .false.

 call obsspace_get_comm(obss, f_comm)

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
 !                       'UFO interface for the CRTM Forward and K-Matrix (setTraj) functions using '//&
 !                       trim(self%conf_traj%ENDIAN_type)//' coefficient datafiles', &
 !                       'CRTM Version: '//TRIM(Version) )


 ! Initialise all the sensors at once
 ! ----------------------------------
 !** NOTE: CRTM_Init points to the various binary files needed for CRTM.  See the
 !**       CRTM_Lifecycle.f90 for more details.

 ! write( *,'(/5x,"Initializing the CRTM (setTraj) ...")' )
 err_stat = CRTM_Init( self%conf_traj%SENSOR_ID                                      , &
                       chinfo                                                        , &
                       File_Path           = trim(self%conf_traj%COEFFICIENT_PATH)   , &
                       NC_File_Path        = trim(self%conf_traj%NC_COEFFICIENT_PATH), &
                       Aerosol_Model       = trim(self%conf_traj%Aerosol_Model)      , &
                       AerosolCoeff_Format = trim(self%conf_traj%AerosolCoeff_Format), &
                       AerosolCoeff_File   = trim(self%conf_traj%AerosolCoeff_File)  , &
                       Cloud_Model         = trim(self%conf_traj%Cloud_Model)        , &
                       CloudCoeff_Format   = trim(self%conf_traj%CloudCoeff_Format)  , &
                       CloudCoeff_File     = trim(self%conf_traj%CloudCoeff_File)    , &
                       IRwaterCoeff_File   = trim(self%conf_traj%IRwaterCoeff_File)  , &
                       IRlandCoeff_File    = trim(self%conf_traj%IRlandCoeff_File)   , &
                       IRsnowCoeff_File    = trim(self%conf_traj%IRsnowCoeff_File)   , &
                       IRiceCoeff_File     = trim(self%conf_traj%IRiceCoeff_File)    , &
                       VISwaterCoeff_File  = trim(self%conf_traj%VISwaterCoeff_File) , &
                       VISlandCoeff_File   = trim(self%conf_traj%VISlandCoeff_File)  , &
                       VISsnowCoeff_File   = trim(self%conf_traj%VISsnowCoeff_File)  , &
                       VISiceCoeff_File    = trim(self%conf_traj%VISiceCoeff_File)   , &
                       MWwaterCoeff_File   = trim(self%conf_traj%MWwaterCoeff_File)  , &
                       Quiet               = .TRUE.)

 message = 'Error initializing CRTM (setTraj)'
 call crtm_comm_stat_check(err_stat, PROGRAM_NAME, message, f_comm)

 ! Loop over all sensors. Not necessary if we're calling CRTM for each sensor
 ! ----------------------------------------------------------------------------
 Sensor_Loop:do n = 1, self%conf_traj%n_Sensors


   ! Pass channel list to CRTM
   ! -------------------------
   err_stat = CRTM_ChannelInfo_Subset(chinfo(n), self%channels, reset=.false.)
   message = 'Error subsetting channels'
   call crtm_comm_stat_check(err_stat, PROGRAM_NAME, message, f_comm)

   ! Determine the number of channels for the current sensor
   ! -------------------------------------------------------
   self%n_Channels = CRTM_ChannelInfo_n_Channels(chinfo(n))


   ! Allocate the ARRAYS
   ! -------------------
   allocate( geo( self%n_Profiles )                         , &
             atm( self%n_Profiles )                         , &
             sfc( self%n_Profiles )                         , &
             rts( self%n_Channels, self%n_Profiles )        , &
             self%atm_K( self%n_Channels, self%n_Profiles ) , &
             self%sfc_K( self%n_Channels, self%n_Profiles ) , &
             rts_K( self%n_Channels, self%n_Profiles )      , &
             self%Options( self%n_Profiles )                , &
             STAT = alloc_stat                                )
   message = 'Error allocating structure arrays (setTraj)'
   call crtm_comm_stat_check(alloc_stat, PROGRAM_NAME, message, f_comm)

   ! Create the input FORWARD structure (atm)
   ! ----------------------------------------
   call CRTM_Atmosphere_Create( atm, self%n_Layers, self%conf_traj%n_Absorbers, self%conf_traj%n_Clouds, self%conf_traj%n_Aerosols )
   if ( ANY(.NOT. CRTM_Atmosphere_Associated(atm)) ) THEN
      message = 'Error allocating CRTM Forward Atmosphere structure (setTraj)'
      CALL Display_Message( PROGRAM_NAME, message, FAILURE )
      STOP
   END IF

   if (self%n_Layers > 0) CALL CRTM_RTSolution_Create(rts, self%n_Layers )

   ! Create the input FORWARD structure (sfc)
   ! ----------------------------------------
   call CRTM_Surface_Create(sfc, self%n_Channels)
   IF ( ANY(.NOT. CRTM_Surface_Associated(sfc)) ) THEN
      message = 'Error allocating CRTM Surface structure (setTraj)'
      CALL Display_Message( PROGRAM_NAME, message, FAILURE )
      STOP
   END IF


   ! Create output K-MATRIX structure (atm)
   ! --------------------------------------
   call CRTM_Atmosphere_Create( self%atm_K, self%n_Layers, self%conf_traj%n_Absorbers, &
                                self%conf_traj%n_Clouds, self%conf_traj%n_Aerosols )
   if ( ANY(.NOT. CRTM_Atmosphere_Associated(self%atm_K)) ) THEN
      message = 'Error allocating CRTM K-matrix Atmosphere structure (setTraj)'
      CALL Display_Message( PROGRAM_NAME, message, FAILURE )
      STOP
   END IF


   ! Create output K-MATRIX structure (sfc)
   ! --------------------------------------
   call CRTM_Surface_Create(self%sfc_K, self%n_Channels)
   IF ( ANY(.NOT. CRTM_Surface_Associated(self%sfc_K)) ) THEN
      message = 'Error allocating CRTM K-matrix Surface structure (setTraj)'
      CALL Display_Message( PROGRAM_NAME, message, FAILURE )
      STOP
   END IF

   ! Some special treatment for Vis or UV as compared to IR or MW.
   ! ----------------------------------------
   if (SpcCoeff_IsVisibleSensor(SC(n)) .or. SpcCoeff_IsUltravioletSensor(SC(n))) then
      Is_Vis_or_UV = .true.
   else
      Is_Vis_or_UV = .false.
   endif

   !Assign the data from the GeoVaLs
   !--------------------------------
   call Load_Atm_Data(self%N_PROFILES,self%N_LAYERS,geovals,atm,self%conf_traj, SC(n)%Is_Active_Sensor)
   call Load_Sfc_Data(self%N_PROFILES,self%n_Channels,self%channels,geovals,sfc,chinfo,obss,self%conf_traj, &
                      SC(n)%Is_Active_Sensor, Is_Vis_or_UV)
   if (cmp_strings(self%conf%SENSOR_ID(n),'gmi_gpm')) then
      allocate( geo_hf( self%n_Profiles ))
      call Load_Geom_Data(obss,geo,geo_hf,self%conf%SENSOR_ID(n))
   else
      call Load_Geom_Data(obss,geo)
   endif

   ! Zero the K-matrix OUTPUT structures
   ! -----------------------------------
   call CRTM_Atmosphere_Zero( self%atm_K )
   call CRTM_Surface_Zero( self%sfc_K )

   ! Inintialize the K-matrix INPUT so that the results are dTb/dx or dR/dx
   ! -------------------------------------------------------------
   if (SC(n)%Is_Active_Sensor) then
       do jchannel = 1, self%n_Channels
          do jprofile = 1, self%n_Profiles
             do jlevel = 1, self%n_Layers
                rts_K(jchannel,jprofile)%Reflectivity(jlevel)            = ZERO
                rts_K(jchannel,jprofile)%Reflectivity_Attenuated(jlevel) = ONE
             end do
          end do
       end do

       rts_K%Radiance                = ZERO
       rts_K%Brightness_Temperature  = ZERO
   else if (Is_Vis_or_UV) then
       rts_K%Radiance                = ONE
       rts_K%Brightness_Temperature  = ZERO
   else
       rts_K%Radiance                = ZERO
       rts_K%Brightness_Temperature  = ONE
   end if


   call ufo_crtm_skip_profiles(self%n_Profiles,self%n_Channels,self%channels,obss,atm,sfc, &
                               SC(n)%Is_Active_Sensor, Is_Vis_or_UV, self%Options)

   ! Call the K-matrix model
   ! -----------------------
   err_stat = CRTM_K_Matrix( atm         , &  ! FORWARD  Input
                             sfc         , &  ! FORWARD  Input
                             rts_K       , &  ! K-MATRIX Input
                             geo         , &  ! Input
                             chinfo(n:n) , &  ! Input
                             self%atm_K  , &  ! K-MATRIX Output
                             self%sfc_K  , &  ! K-MATRIX Output
                             rts         , &  ! FORWARD  Output
                             self%Options  )  ! Input
   message = 'Error calling CRTM (setTraj) K-Matrix Model for '//TRIM(self%conf_traj%SENSOR_ID(n))
   call crtm_comm_stat_check(err_stat, PROGRAM_NAME, message, f_comm)
   if (cmp_strings(self%conf%SENSOR_ID(n),'gmi_gpm')) then
      allocate( atm_Ka( self%n_Channels, self%n_Profiles ),               &
                sfc_Ka( self%n_Channels, self%n_Profiles ),   &
                rts_Ka( self%n_Channels, self%n_Profiles ),   &
                rtsa( self%n_Channels, self%n_Profiles ),     &
                STAT = alloc_stat )
      message = 'Error allocating K structure arrays rtsa, atm_Ka ......'
      call crtm_comm_stat_check(alloc_stat, PROGRAM_NAME, message, f_comm)
      !! save resutls for gmi channels 1-9.
      atm_Ka = self%atm_K
      sfc_Ka = self%sfc_K
      rts_Ka = rts_K
      rtsa   = rts
      ! Zero the K-matrix OUTPUT structures
      ! -----------------------------------
      call CRTM_Atmosphere_Zero( self%atm_K )
      call CRTM_Surface_Zero( self%sfc_K )
      ! Inintialize the K-matrix INPUT so that the results are dTb/dx
      ! -------------------------------------------------------------
      rts_K%Radiance                = ZERO
      rts_K%Brightness_Temperature  = ONE

      ! Call the K-matrix model
      ! -----------------------
      err_stat = CRTM_K_Matrix( atm         , &  ! FORWARD  Input
                                sfc         , &  ! FORWARD  Input
                                rts_K       , &  ! K-MATRIX Input
                                geo_hf        , &  ! Input
                                chinfo(n:n) , &  ! Input
                                self%atm_K  , &  ! K-MATRIX Output
                                self%sfc_K  , &  ! K-MATRIX Output
                                rts         , &  ! FORWARD  Output
                                self%Options  )  ! Input
      message = 'Error calling CRTM (setTraj, geo_hf) K-Matrix Model for '&
                //TRIM(self%conf_traj%SENSOR_ID(n))
      call crtm_comm_stat_check(err_stat, PROGRAM_NAME, message, f_comm)
      !! replace data for gmi channels 1-9 by early results calculated with geo.
      do lch = 1, size(self%channels)
         if ( self%channels(lch) <= 9 ) then
            self%atm_K(lch,:) = atm_Ka(lch,:)
            self%sfc_K(lch,:) = sfc_Ka(lch,:)
            rts_K(lch,:) = rts_Ka(lch,:)
            rts(lch,:)   = rtsa(lch,:)
         endif
      enddo
      deallocate(atm_Ka,sfc_Ka,rts_Ka,rtsa)
   endif ! cmp_strings(self%conf%SENSOR_ID(n),'gmi_gpm')

   !call CRTM_RTSolution_Inspect(rts)

   ! check for NaN values in atm_k
   numNaN = 0
   do jprofile = 1, self%n_Profiles
      do jchannel = 1, size(self%channels)
         do jlevel = 1, self%atm_K(jchannel,jprofile)%n_layers
            if (ieee_is_nan(self%atm_K(jchannel,jprofile)%Temperature(jlevel))) then
               self%Options(jprofile)%Skip_Profile = .TRUE.
               numNaN = numNaN + 1
               write(message,*) numNaN, 'th NaN in Jacobian Profiles'
               call fckit_log%info(message)
               cycle
            end if
         end do
      end do
   end do

   !! Parse hofxdiags%variables into independent/dependent variables and channel
   !! assumed formats:
   !!   jacobian var -->     <ystr>_jacobian_<xstr>_<chstr>
   !!   non-jacobian var --> <ystr>_<chstr>

   ch_diags = -9999
   do jvar = 1, hofxdiags%nvar
      varstr = hofxdiags%variables(jvar)
      str_pos(4) = len_trim(varstr)
      if (str_pos(4) < 1) cycle
      str_pos(3) = index(varstr,"_",back=.true.)        !final "_" before channel
      read(varstr(str_pos(3)+1:str_pos(4)),*, err=999) ch_diags(jvar)
 999  str_pos(1) = index(varstr,jacobianstr) - 1        !position before jacobianstr
      if (str_pos(1) == 0) then
         write(err_msg,*) 'ufo_crtm_passive_tlad_diags: _jacobian_ must be // &
                           & preceded by dependent variable in config: ', &
                           & hofxdiags%variables(jvar)
         call abor1_ftn(err_msg)
      else if (str_pos(1) > 0) then
         !Diagnostic is a Jacobian member (dy/dx)
         ystr_diags(jvar) = varstr(1:str_pos(1))
         str_pos(2) = str_pos(1) + len(jacobianstr) + 1 !begin xstr_diags
         str_pos(4) = str_pos(3) - str_pos(2)
         xstr_diags(jvar)(1:str_pos(4)) = varstr(str_pos(2):str_pos(3)-1)
         xstr_diags(jvar)(str_pos(4)+1:) = ""
      else !null
         !Diagnostic is a dependent variable (y)
         xstr_diags(jvar) = ""
         ystr_diags(jvar)(1:str_pos(3)-1) = varstr(1:str_pos(3)-1)
         ystr_diags(jvar)(str_pos(3):) = ""
         if (ch_diags(jvar) < 0) ystr_diags(jvar) = varstr
      end if
   end do

   nvars = hofxdiags%nvar
   nlocs = geovals%nlocs
   ! put CRTM output into hofxdiags
   if (SC(n)%Is_Active_Sensor) then
      call ufo_crtm_active_diag(rts, &
                                rts_K, &
                                atm, &
                                self%atm_K, &
                                self%sfc_K, &
                                self%conf, &
                                n, &
                                self%Options, &
                                self%channels, &
                                geovals, &
                                obss, &
                                nvars, &
                                nlocs, &
                                self%n_Profiles, &
                                self%n_Layers, &
                                xstr_diags, &
                                ystr_diags, &
                                ch_diags, &
                                hofxdiags,&
                                err_stat)
   else
      call ufo_crtm_passive_diag(rts, &
                                 rts_K, &
                                 atm, &
                                 self%atm_K, &
                                 self%sfc_K, &
                                 self%conf, &
                                 n, &
                                 self%Options, &
                                 self%channels, &
                                 geovals, &
                                 obss, &
                                 nvars, &
                                 nlocs, &
                                 self%n_Profiles, &
                                 self%n_Layers, &
                                 xstr_diags, &
                                 ystr_diags, &
                                 ch_diags, &
                                 hofxdiags,&
                                 err_stat)


!       call ufo_crtm_passive_tlad_diags(rts, rts_K, atm, obss, &
!                                        self%conf, &
!                                        self%Options, self%channels, &
!                                        self%n_layers, self%n_profiles, &
!                                        hofxdiags, xstr_diags, ystr_diags, &
!                                        ch_diags, err_stat)
   end if

   ! check for error from either passive or active
   if (err_stat > 0) then
       write(err_msg,*) 'ufo_radiancecrtm_tlad_settraj error: failed to put simulated diagnostics into hofxdiags'
       call abor1_ftn(err_msg)
    end if

   ! Deallocate the structures
   ! -------------------------
   call CRTM_Geometry_Destroy(geo)
   call CRTM_Atmosphere_Destroy(atm)
   call CRTM_RTSolution_Destroy(rts_K)
   call CRTM_RTSolution_Destroy(rts)
   call CRTM_Surface_Destroy(sfc)


   ! Deallocate all arrays
   ! ---------------------
   deallocate(geo, atm, sfc, rts, rts_K, STAT = alloc_stat)
   if(allocated(geo_hf)) deallocate(geo_hf)
   message = 'Error deallocating structure arrays (setTraj)'
   call crtm_comm_stat_check(alloc_stat, PROGRAM_NAME, message, f_comm)

 end do Sensor_Loop


 ! Destroy CRTM instance
 ! ---------------------
 ! write( *, '( /5x, "Destroying the CRTM (setTraj)..." )' )
 err_stat = CRTM_Destroy( chinfo )
 message = 'Error destroying CRTM (setTraj)'
 call crtm_comm_stat_check(err_stat, PROGRAM_NAME, message, f_comm)

 ! Set flag that the tracectory was set
 ! ------------------------------------
 self%ltraj = .true.

end subroutine ufo_radiancecrtm_tlad_settraj

! ------------------------------------------------------------------------------

subroutine ufo_radiancecrtm_simobs_tl(self, geovals, obss, nvars, nlocs, hofx, qc_flags_ptr)

implicit none
class(ufo_radiancecrtm_tlad), intent(in)    :: self
type(ufo_geovals),        intent(in)    :: geovals
type(c_ptr), value,       intent(in)    :: obss
integer,                  intent(in)    :: nvars, nlocs
real(c_double),           intent(inout) :: hofx(nvars, nlocs)
type(c_ptr), value, intent(in) :: qc_flags_ptr
type(obsdatavector_int) :: qc_flags
character(len=*), parameter :: myname_="ufo_radiancecrtm_simobs_tl"
character(max_string) :: err_msg
integer :: jprofile, jchannel, jlevel, jspec, ispec
type(ufo_geoval), pointer :: geoval_d

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

 qc_flags%data_ptr = qc_flags_ptr
 ! Initialize hofx
 ! ---------------
 hofx(:,:) = 0.0_kind_real

 ! Temperature
 ! -----------

 ! Get t from geovals
 call ufo_geovals_get_var(geovals, var_ts, geoval_d)

 ! Check model levels is consistent in geovals & crtm
 if (geoval_d%nval /= self%n_Layers) then
   write(err_msg,*) myname_, ' error: layers inconsistent!'
   call abor1_ftn(err_msg)
 endif

 if ( .not.  self%use_qc_flags) then
   ! Multiply by Jacobian and add to hofx
   do jprofile = 1, self%n_Profiles
     if (.not.self%Options(jprofile)%Skip_Profile) then
       do jchannel = 1, size(self%channels) 
         do jlevel = 1, geoval_d%nval
           hofx(jchannel, jprofile) = hofx(jchannel, jprofile) + &
                      self%atm_K(jchannel,jprofile)%Temperature(jlevel) * &
                      geoval_d%vals(jlevel,jprofile)
         enddo
       enddo
     end if
   enddo
 else
   ! Multiply by Jacobian and add to hofx
   do jprofile = 1, self%n_Profiles
     if (.not.self%Options(jprofile)%Skip_Profile) then
       do jchannel = 1, size(self%channels) 
         if ( qc_flags%get(int(jchannel,int64),int(jprofile,int64)) > 1 ) then
           cycle 
         end if
         do jlevel = 1, geoval_d%nval
           hofx(jchannel, jprofile) = hofx(jchannel, jprofile) + &
                      self%atm_K(jchannel,jprofile)%Temperature(jlevel) * &
                      geoval_d%vals(jlevel,jprofile)
         enddo
       enddo
     end if
   enddo
 end if

 ! Absorbers
 ! ---------

 if ( .not.  self%use_qc_flags) then

   do jspec = 1, self%conf%n_Absorbers
     ! Get Absorber from geovals
     call ufo_geovals_get_var(geovals, self%conf%Absorbers(jspec), geoval_d)

     ispec = ufo_vars_getindex(self%conf_traj%Absorbers, self%conf%Absorbers(jspec))

     ! Multiply by Jacobian and add to hofx
     do jprofile = 1, self%n_Profiles
       if (.not.self%Options(jprofile)%Skip_Profile) then
         do jchannel = 1, size(self%channels)
           do jlevel = 1, geoval_d%nval
             hofx(jchannel, jprofile) = hofx(jchannel, jprofile) + &
                        self%atm_K(jchannel,jprofile)%Absorber(jlevel,ispec) * &
                        geoval_d%vals(jlevel,jprofile)
           enddo
         enddo
       end if
     enddo
   end do

   ! Clouds (mass content only)
   ! --------------------------

   do jspec = 1, self%conf%n_Clouds
     ! Get Cloud from geovals
     call ufo_geovals_get_var(geovals, self%conf%Clouds(jspec,1), geoval_d)

     ispec = ufo_vars_getindex(self%conf_traj%Clouds(:,1), self%conf%Clouds(jspec,1))

     ! Multiply by Jacobian and add to hofx
     do jprofile = 1, self%n_Profiles
       if (.not.self%Options(jprofile)%Skip_Profile) then
         do jchannel = 1, size(self%channels)
           do jlevel = 1, geoval_d%nval
             hofx(jchannel, jprofile) = hofx(jchannel, jprofile) + &
                        self%atm_K(jchannel,jprofile)%Cloud(ispec)%Water_Content(jlevel) * &
                        geoval_d%vals(jlevel,jprofile)
           enddo
         enddo
       end if
     enddo
   end do

 else
   do jspec = 1, self%conf%n_Absorbers
     ! Get Absorber from geovals
     call ufo_geovals_get_var(geovals, self%conf%Absorbers(jspec), geoval_d)

     ispec = ufo_vars_getindex(self%conf_traj%Absorbers, self%conf%Absorbers(jspec))

     ! Multiply by Jacobian and add to hofx
     do jprofile = 1, self%n_Profiles
       if (.not.self%Options(jprofile)%Skip_Profile) then
         do jchannel = 1, size(self%channels)
           if ( qc_flags%get(int(jchannel,int64),int(jprofile,int64)) > 1) cycle
           do jlevel = 1, geoval_d%nval
             hofx(jchannel, jprofile) = hofx(jchannel, jprofile) + &
                        self%atm_K(jchannel,jprofile)%Absorber(jlevel,ispec) * &
                        geoval_d%vals(jlevel,jprofile)
           enddo
         enddo
       end if
     enddo
   end do

   ! Clouds (mass content only)
   ! --------------------------

   do jspec = 1, self%conf%n_Clouds
     ! Get Cloud from geovals
     call ufo_geovals_get_var(geovals, self%conf%Clouds(jspec,1), geoval_d)

     ispec = ufo_vars_getindex(self%conf_traj%Clouds(:,1), self%conf%Clouds(jspec,1))

     ! Multiply by Jacobian and add to hofx
     do jprofile = 1, self%n_Profiles
       if (.not.self%Options(jprofile)%Skip_Profile) then
         do jchannel = 1, size(self%channels)
           if ( qc_flags%get(int(jchannel,int64),int(jprofile,int64)) > 1) cycle
           do jlevel = 1, geoval_d%nval
             hofx(jchannel, jprofile) = hofx(jchannel, jprofile) + &
                        self%atm_K(jchannel,jprofile)%Cloud(ispec)%Water_Content(jlevel) * &
                        geoval_d%vals(jlevel,jprofile)
           enddo
         enddo
       end if
     enddo
   end do
 end if


 ! Surface Variables
 ! --------------------------

 do jspec = 1, self%conf%n_Surfaces
   ! Get Surface from geovals
   call ufo_geovals_get_var(geovals, self%conf%Surfaces(jspec), geoval_d)

   select case(self%conf%Surfaces(jspec))

      ! Surface Temperature Where Sea / Water Temperature
      case(var_sfc_wtmp)

         ! Multiply by Jacobian and add to hofx
         do jprofile = 1, self%n_Profiles
            if (.not.self%Options(jprofile)%Skip_Profile) then
               do jchannel = 1, size(self%channels)
                  jlevel = 1
                  hofx(jchannel, jprofile) = hofx(jchannel, jprofile) + &
                        self%sfc_K(jchannel,jprofile)%water_temperature * &
                        geoval_d%vals(jlevel,jprofile)
               enddo
            end if
         enddo

      ! Surface Temperature Where Land / Land Temperature
      case(var_sfc_ltmp)

         ! Multiply by Jacobian and add to hofx
         do jprofile = 1, self%n_Profiles
            if (.not.self%Options(jprofile)%Skip_Profile) then
               do jchannel = 1, size(self%channels)
                  jlevel = 1
                  hofx(jchannel, jprofile) = hofx(jchannel, jprofile) + &
                        self%sfc_K(jchannel,jprofile)%land_temperature * &
                        geoval_d%vals(jlevel,jprofile)
               enddo
            end if
         enddo

      ! Surface Temperature Where Ice / Ice Temperature
      case(var_sfc_itmp)

         ! Multiply by Jacobian and add to hofx
         do jprofile = 1, self%n_Profiles
            if (.not.self%Options(jprofile)%Skip_Profile) then
               do jchannel = 1, size(self%channels)
                  jlevel = 1
                  hofx(jchannel, jprofile) = hofx(jchannel, jprofile) + &
                        self%sfc_K(jchannel,jprofile)%ice_temperature * &
                        geoval_d%vals(jlevel,jprofile)
               enddo
            end if
         enddo

      ! Surface Temperature Where Snow / Snow Temperature
      case(var_sfc_stmp)

         ! Multiply by Jacobian and add to hofx
         do jprofile = 1, self%n_Profiles
            if (.not.self%Options(jprofile)%Skip_Profile) then
               do jchannel = 1, size(self%channels)
                  jlevel = 1
                  hofx(jchannel, jprofile) = hofx(jchannel, jprofile) + &
                        self%sfc_K(jchannel,jprofile)%snow_temperature * &
                        geoval_d%vals(jlevel,jprofile)
               enddo
            end if
         enddo

      ! Surface Wind Speed / Wind Speed
      case(var_sfc_wspeed)

         ! Multiply by Jacobian and add to hofx
         do jprofile = 1, self%n_Profiles
            if (.not.self%Options(jprofile)%Skip_Profile) then
               do jchannel = 1, size(self%channels)
                  jlevel = 1
                  hofx(jchannel, jprofile) = hofx(jchannel, jprofile) + &
                        self%sfc_K(jchannel,jprofile)%wind_speed * &
                        geoval_d%vals(jlevel,jprofile)
               enddo
            end if
         enddo

      ! Surface Wind From Direction / Wind Direction
      case(var_sfc_wdir)

         ! Multiply by Jacobian and add to hofx
         do jprofile = 1, self%n_Profiles
            if (.not.self%Options(jprofile)%Skip_Profile) then
               do jchannel = 1, size(self%channels)
                  jlevel = 1
                  hofx(jchannel, jprofile) = hofx(jchannel, jprofile) + &
                        self%sfc_K(jchannel,jprofile)%wind_direction * &
                        geoval_d%vals(jlevel,jprofile)
               enddo
            end if
         enddo

      ! Sea Surface Salinity / Salinity
      case(var_sfc_sss)

         ! Multiply by Jacobian and add to hofx
         do jprofile = 1, self%n_Profiles
            if (.not.self%Options(jprofile)%Skip_Profile) then
               do jchannel = 1, size(self%channels)
                  jlevel = 1
                  hofx(jchannel, jprofile) = hofx(jchannel, jprofile) + &
                        self%sfc_K(jchannel,jprofile)%salinity * &
                        geoval_d%vals(jlevel,jprofile)
               enddo
            end if
         enddo

   end select
 end do

end subroutine ufo_radiancecrtm_simobs_tl


! ------------------------------------------------------------------------------

subroutine ufo_radiancecrtm_simobs_ad(self, geovals, obss, nvars, nlocs, hofx, qc_flags_ptr)

implicit none
class(ufo_radiancecrtm_tlad), intent(in)    :: self
type(ufo_geovals),        intent(inout) :: geovals
type(c_ptr), value,       intent(in)    :: obss
integer,                  intent(in)    :: nvars, nlocs
real(c_double),           intent(in)    :: hofx(nvars, nlocs)
type(c_ptr), value,       intent(in) :: qc_flags_ptr
type(obsdatavector_int) :: qc_flags
character(len=*), parameter :: myname_="ufo_radiancecrtm_simobs_ad"
character(max_string) :: err_msg
integer :: jprofile, jchannel, jlevel, jspec, ispec
type(ufo_geoval), pointer :: geoval_d
real(c_double) :: missing


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

 qc_flags%data_ptr = qc_flags_ptr
 ! Temperature
 ! -----------

 ! Get t from geovals
 call ufo_geovals_get_var(geovals, var_ts, geoval_d)

 if (.not. self%use_qc_flags) then

   ! Multiply by Jacobian and add to hofx (adjoint)
   do jprofile = 1, self%n_Profiles
     if (.not.self%Options(jprofile)%Skip_Profile) then
       do jchannel = 1, size(self%channels)
         if (hofx(jchannel, jprofile) /= missing) then
           do jlevel = 1, geoval_d%nval
               geoval_d%vals(jlevel,jprofile) = geoval_d%vals(jlevel,jprofile) + &
                                          self%atm_K(jchannel,jprofile)%Temperature(jlevel) * &
                                          hofx(jchannel, jprofile)
           enddo
         endif
       enddo
     end if
   enddo

   ! Absorbers
   ! ---------

   do jspec = 1, self%conf%n_Absorbers
   ! Get Absorber from geovals
     call ufo_geovals_get_var(geovals, self%conf%Absorbers(jspec), geoval_d)

     ispec = ufo_vars_getindex(self%conf_traj%Absorbers, self%conf%Absorbers(jspec))

     ! Multiply by Jacobian and add to hofx (adjoint)
     do jprofile = 1, self%n_Profiles
       if (.not.self%Options(jprofile)%Skip_Profile) then
         do jchannel = 1, size(self%channels)
           if (hofx(jchannel, jprofile) /= missing) then
             do jlevel = 1, geoval_d%nval
               geoval_d%vals(jlevel,jprofile) = geoval_d%vals(jlevel,jprofile) + &
                                            self%atm_K(jchannel,jprofile)%Absorber(jlevel,ispec) * &
                                            hofx(jchannel, jprofile)
             enddo
           endif
         enddo
       end if
     enddo
   end do

   ! Clouds (mass content only)
   ! --------------------------

   do jspec = 1, self%conf%n_Clouds
     ! Get Cloud from geovals
     call ufo_geovals_get_var(geovals, self%conf%Clouds(jspec,1), geoval_d)

     ispec = ufo_vars_getindex(self%conf_traj%Clouds(:,1), self%conf%Clouds(jspec,1))

     ! Multiply by Jacobian and add to hofx (adjoint)
     do jprofile = 1, self%n_Profiles
       if (.not.self%Options(jprofile)%Skip_Profile) then
         do jchannel = 1, size(self%channels)
           if (hofx(jchannel, jprofile) /= missing) then
             do jlevel = 1, geoval_d%nval
               geoval_d%vals(jlevel,jprofile) = geoval_d%vals(jlevel,jprofile) + &
                                            self%atm_K(jchannel,jprofile)%Cloud(ispec)%Water_Content(jlevel) * &
                                            hofx(jchannel, jprofile)
             enddo
           endif
         enddo
       end if
     enddo
   end do
 else

   ! Multiply by Jacobian and add to hofx (adjoint)
   do jprofile = 1, self%n_Profiles
     if (.not.self%Options(jprofile)%Skip_Profile) then
       do jchannel = 1, size(self%channels)
         if (qc_flags%get(int(jchannel,int64),int(jprofile,int64))>1) cycle
         if (hofx(jchannel, jprofile) /= missing) then
           do jlevel = 1, geoval_d%nval
               geoval_d%vals(jlevel,jprofile) = geoval_d%vals(jlevel,jprofile) + &
                                          self%atm_K(jchannel,jprofile)%Temperature(jlevel) * &
                                          hofx(jchannel, jprofile)
           enddo
         endif
       enddo
     end if
   enddo

   ! Absorbers
   ! ---------

   do jspec = 1, self%conf%n_Absorbers
   ! Get Absorber from geovals
     call ufo_geovals_get_var(geovals, self%conf%Absorbers(jspec), geoval_d)

     ispec = ufo_vars_getindex(self%conf_traj%Absorbers, self%conf%Absorbers(jspec))

     ! Multiply by Jacobian and add to hofx (adjoint)
     do jprofile = 1, self%n_Profiles
       if (.not.self%Options(jprofile)%Skip_Profile) then
         do jchannel = 1, size(self%channels)
           if (qc_flags%get(int(jchannel,int64),int(jprofile,int64))>1) cycle
           if (hofx(jchannel, jprofile) /= missing) then
             do jlevel = 1, geoval_d%nval
               geoval_d%vals(jlevel,jprofile) = geoval_d%vals(jlevel,jprofile) + &
                                            self%atm_K(jchannel,jprofile)%Absorber(jlevel,ispec) * &
                                            hofx(jchannel, jprofile)
             enddo
           endif
         enddo
       end if
     enddo
   end do

   ! Clouds (mass content only)
   ! --------------------------

   do jspec = 1, self%conf%n_Clouds
     ! Get Cloud from geovals
     call ufo_geovals_get_var(geovals, self%conf%Clouds(jspec,1), geoval_d)

     ispec = ufo_vars_getindex(self%conf_traj%Clouds(:,1), self%conf%Clouds(jspec,1))

     ! Multiply by Jacobian and add to hofx (adjoint)
     do jprofile = 1, self%n_Profiles
       if (.not.self%Options(jprofile)%Skip_Profile) then
         do jchannel = 1, size(self%channels)
           if (qc_flags%get(int(jchannel,int64),int(jprofile,int64))>1) cycle
           if (hofx(jchannel, jprofile) /= missing) then
             do jlevel = 1, geoval_d%nval
               geoval_d%vals(jlevel,jprofile) = geoval_d%vals(jlevel,jprofile) + &
                                            self%atm_K(jchannel,jprofile)%Cloud(ispec)%Water_Content(jlevel) * &
                                            hofx(jchannel, jprofile)
             enddo
           endif
         enddo
       end if
     enddo
   end do
 end if 

 ! Surface Variables
 ! --------------------------
 do jspec = 1, self%conf%n_Surfaces
   ! Get Cloud from geovals
   call ufo_geovals_get_var(geovals, self%conf%Surfaces(jspec), geoval_d)

   select case(self%conf%Surfaces(jspec))

      ! Surface Temperature Where Sea / Water Temperature
      case(var_sfc_wtmp)

         ! Multiply by Jacobian and add to hofx
         do jprofile = 1, self%n_Profiles
            if (.not.self%Options(jprofile)%Skip_Profile) then
               do jchannel = 1, size(self%channels)
                  if (hofx(jchannel, jprofile) /= missing) then
                     jlevel = 1
                     geoval_d%vals(jlevel, jprofile) = geoval_d%vals(jlevel,jprofile) + &
                          self%sfc_K(jchannel,jprofile)%water_temperature * &
                          hofx(jchannel,jprofile)
                  endif
               enddo
            end if
         enddo

      ! Surface Temperature Where Land / Land Temperature
      case(var_sfc_ltmp)

         ! Multiply by Jacobian and add to hofx
         do jprofile = 1, self%n_Profiles
            if (.not.self%Options(jprofile)%Skip_Profile) then
               do jchannel = 1, size(self%channels)
                  if (hofx(jchannel, jprofile) /= missing) then
                     jlevel = 1
                     geoval_d%vals(jlevel, jprofile) = geoval_d%vals(jlevel,jprofile) + &
                          self%sfc_K(jchannel,jprofile)%land_temperature * &
                          hofx(jchannel,jprofile)
                  endif
               enddo
            end if
         enddo

      ! Surface Temperature Where Ice / Ice Temperature
      case(var_sfc_itmp)

         ! Multiply by Jacobian and add to hofx
         do jprofile = 1, self%n_Profiles
            if (.not.self%Options(jprofile)%Skip_Profile) then
               do jchannel = 1, size(self%channels)
                  if (hofx(jchannel, jprofile) /= missing) then
                     jlevel = 1
                     geoval_d%vals(jlevel, jprofile) = geoval_d%vals(jlevel,jprofile) + &
                          self%sfc_K(jchannel,jprofile)%ice_temperature * &
                          hofx(jchannel,jprofile)
                  endif
               enddo
            end if
         enddo

      ! Surface Temperature Where Snow / Snow Temperature
      case(var_sfc_stmp)

         ! Multiply by Jacobian and add to hofx
         do jprofile = 1, self%n_Profiles
            if (.not.self%Options(jprofile)%Skip_Profile) then
               do jchannel = 1, size(self%channels)
                  if (hofx(jchannel, jprofile) /= missing) then
                     jlevel = 1
                     geoval_d%vals(jlevel, jprofile) = geoval_d%vals(jlevel,jprofile) + &
                          self%sfc_K(jchannel,jprofile)%snow_temperature * &
                          hofx(jchannel,jprofile)
                  endif
               enddo
            end if
         enddo

      ! Surface Wind Speed / Wind Speed
      case(var_sfc_wspeed)

         ! Multiply by Jacobian and add to hofx
         do jprofile = 1, self%n_Profiles
            if (.not.self%Options(jprofile)%Skip_Profile) then
               do jchannel = 1, size(self%channels)
                  if (hofx(jchannel, jprofile) /= missing) then
                     jlevel = 1
                     geoval_d%vals(jlevel, jprofile) = geoval_d%vals(jlevel,jprofile) + &
                          self%sfc_K(jchannel,jprofile)%wind_speed * &
                          hofx(jchannel,jprofile)
                  endif
               enddo
            end if
         enddo

      ! Surface Wind From Direction / Wind Direction
      case(var_sfc_wdir)

         ! Multiply by Jacobian and add to hofx
         do jprofile = 1, self%n_Profiles
            if (.not.self%Options(jprofile)%Skip_Profile) then
               do jchannel = 1, size(self%channels)
                  if (hofx(jchannel, jprofile) /= missing) then
                     jlevel = 1
                     geoval_d%vals(jlevel, jprofile) = geoval_d%vals(jlevel,jprofile) + &
                          self%sfc_K(jchannel,jprofile)%wind_direction * &
                          hofx(jchannel,jprofile)
                  endif
               enddo
            end if
         enddo

      ! Sea Surface Salinity / Salinity
      case(var_sfc_sss)

         ! Multiply by Jacobian and add to hofx
         do jprofile = 1, self%n_Profiles
            if (.not.self%Options(jprofile)%Skip_Profile) then
               do jchannel = 1, size(self%channels)
                  if (hofx(jchannel, jprofile) /= missing) then
                     jlevel = 1
                     geoval_d%vals(jlevel, jprofile) = geoval_d%vals(jlevel,jprofile) + &
                          self%sfc_K(jchannel,jprofile)%salinity * &
                          hofx(jchannel,jprofile)
                  endif
               enddo
            end if
         enddo

   end select

 enddo


end subroutine ufo_radiancecrtm_simobs_ad

! ------------------------------------------------------------------------------

end module ufo_radiancecrtm_tlad_mod
