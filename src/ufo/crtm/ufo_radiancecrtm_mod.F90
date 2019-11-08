! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module to handle radiancecrtm observations

module ufo_radiancecrtm_mod

 use fckit_configuration_module, only: fckit_configuration
 use iso_c_binding
 use kinds

 use ufo_geovals_mod, only: ufo_geovals, ufo_geoval, ufo_geovals_get_var
 use ufo_vars_mod
 use ufo_crtm_utils_mod
 use crtm_module
 use obsspace_mod

 implicit none
 private

 !> Fortran derived type for radiancecrtm trajectory
 type, public :: ufo_radiancecrtm
 private
   character(len=MAXVARLEN), public, allocatable :: varin(:)  ! variables requested from the model
   integer, allocatable                          :: channels(:)
   type(crtm_conf) :: conf
 contains
   procedure :: setup  => ufo_radiancecrtm_setup
   procedure :: delete => ufo_radiancecrtm_delete
   procedure :: simobs => ufo_radiancecrtm_simobs
 end type ufo_radiancecrtm

 character(len=maxvarlen), dimension(21), parameter :: varin_default = &
                            (/var_ts, var_prs, var_prsi,                                  &
                              var_sfc_wfrac, var_sfc_lfrac, var_sfc_ifrac, var_sfc_sfrac, &
                              var_sfc_wtmp,  var_sfc_ltmp,  var_sfc_itmp,  var_sfc_stmp,  &
                              var_sfc_vegfrac, var_sfc_wspeed, var_sfc_wdir, var_sfc_lai, &
                              var_sfc_soilm, var_sfc_soilt, var_sfc_landtyp,              &
                              var_sfc_vegtyp, var_sfc_soiltyp, var_sfc_sdepth/)

contains

! ------------------------------------------------------------------------------

subroutine ufo_radiancecrtm_setup(self, f_confOper, channels)

implicit none
class(ufo_radiancecrtm),   intent(inout) :: self
type(fckit_configuration), intent(in)    :: f_confOper
integer(c_int),            intent(in)    :: channels(:)  !List of channels to use

integer :: nvars_in
integer :: ind, jspec
character(len=max_string) :: err_msg
character(len=:), allocatable :: str
type(fckit_configuration) :: f_confOpts

 call f_confOper%get_or_die("ObsOptions",f_confOpts)

 call crtm_conf_setup(self%conf,f_confOpts,f_confOper)
 if ( ufo_vars_getindex(self%conf%Absorbers, var_mixr) < 1 ) then
   write(err_msg,*) 'ufo_radiancecrtm_setup error: H2O must be included in CRTM Absorbers'
   call abor1_ftn(err_msg)
 end if
 if ( ufo_vars_getindex(self%conf%Absorbers, var_oz) < 1 ) then
   write(err_msg,*) 'ufo_radiancecrtm_setup error: O3 must be included in CRTM Absorbers'
   call abor1_ftn(err_msg)
 end if

 ! request from the model all the hardcoded atmospheric & surface variables + 
 ! 1 * n_Absorbers
 ! 2 * n_Clouds (mass content and effective radius)
 ! if sss is in ObsOptions + sss
 nvars_in = size(varin_default) + self%conf%n_Absorbers + 2 * self%conf%n_Clouds
 
 if (f_confOpts%has("SalinityOption")) then
    nvars_in =  nvars_in + 1
    allocate(self%varin(nvars_in))
    call f_confOpts%get_or_die("SalinityOption",str)
    if (str =="sss") then 
       self%varin(nvars_in) = var_sfc_sss
    endif
 else
    allocate(self%varin(nvars_in))
 endif

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
   self%varin(ind) = self%conf%Clouds(jspec,2)
   ind = ind + 1
 end do

 ! save channels
 allocate(self%channels(size(channels)))
 self%channels(:) = channels(:)

end subroutine ufo_radiancecrtm_setup

! ------------------------------------------------------------------------------

subroutine ufo_radiancecrtm_delete(self)

implicit none
class(ufo_radiancecrtm), intent(inout) :: self

 call crtm_conf_delete(self%conf)

end subroutine ufo_radiancecrtm_delete

! ------------------------------------------------------------------------------

subroutine ufo_radiancecrtm_simobs(self, geovals, obss, nvars, nlocs, hofx, hofxdiags)

implicit none
class(ufo_radiancecrtm),  intent(in) :: self         !Radiance object
type(ufo_geovals),        intent(in) :: geovals      !Inputs from the model
integer,                  intent(in) :: nvars, nlocs
real(c_double),        intent(inout) :: hofx(nvars, nlocs) !h(x) to return
type(ufo_geovals),     intent(inout) :: hofxdiags    !non-h(x) diagnostics
type(c_ptr), value,       intent(in) :: obss         !ObsSpace

! Local Variables
character(*), parameter :: PROGRAM_NAME = 'ufo_radiancecrtm_mod.F90'
character(255) :: message, version
integer        :: err_stat, alloc_stat
integer        :: l, m, n
type(ufo_geoval), pointer :: temp

integer :: n_Profiles
integer :: n_Layers
integer :: n_Channels

! Define the "non-demoninational" arguments
type(CRTM_ChannelInfo_type)             :: chinfo(self%conf%n_Sensors)
type(CRTM_Geometry_type),   allocatable :: geo(:)

! Define the FORWARD variables
type(CRTM_Atmosphere_type), allocatable :: atm(:)
type(CRTM_Surface_type),    allocatable :: sfc(:)
type(CRTM_RTSolution_type), allocatable :: rts(:,:)
type(CRTM_Options_type), allocatable    :: Options(:) 

! Define the K-MATRIX variables for hofxdiags
type(CRTM_Atmosphere_type), allocatable :: atm_K(:,:)
type(CRTM_Surface_type),    allocatable :: sfc_K(:,:)
type(CRTM_RTSolution_type), allocatable :: rts_K(:,:)

! Used to parse hofxdiags
character(len=MAXVARLEN) :: varstr
character(len=MAXVARLEN), dimension(hofxdiags%nvar) :: &
                          ystr_diags, xstr_diags
character(10), parameter :: jacobianstr = "_jacobian_"
integer :: str_pos(4), ch_diags(hofxdiags%nvar)
integer :: jvar, jprofile, jlevel, jchannel, ichannel, jspec

logical :: jacobian_needed
character(max_string) :: err_msg

 ! Get number of profile and layers from geovals
 ! ---------------------------------------------
 n_Profiles = geovals%nlocs
 call ufo_geovals_get_var(geovals, var_ts, temp)
 n_Layers = temp%nval
 nullify(temp)


 ! Program header
 ! --------------
 call CRTM_Version( Version )
 call Program_Message( PROGRAM_NAME, &
                       'Check/example program for the CRTM Forward and K-Matrix functions using '//&
                       trim(self%conf%ENDIAN_type)//' coefficient datafiles', &
                       'CRTM Version: '//TRIM(Version) )


 ! Initialise all the sensors at once
 ! ----------------------------------
 !** NOTE: CRTM_Init points to the various binary files needed for CRTM.  See the
 !**       CRTM_Lifecycle.f90 for more details.

 write( *,'(/5x,"Initializing the CRTM...")' )
 err_stat = CRTM_Init( self%conf%SENSOR_ID, chinfo, &
                       File_Path=trim(self%conf%COEFFICIENT_PATH), &
                       IRwaterCoeff_File=trim(self%conf%IRwaterCoeff_File), &
                       IRlandCoeff_File=trim(self%conf%IRlandCoeff_File), &
                       IRsnowCoeff_File=trim(self%conf%IRsnowCoeff_File), &
                       IRiceCoeff_File=trim(self%conf%IRiceCoeff_File), &
                       VISwaterCoeff_File=trim(self%conf%VISwaterCoeff_File), &
                       VISlandCoeff_File=trim(self%conf%VISlandCoeff_File), &
                       VISsnowCoeff_File=trim(self%conf%VISsnowCoeff_File), &
                       VISiceCoeff_File=trim(self%conf%VISiceCoeff_File), &
                       MWwaterCoeff_File=trim(self%conf%MWwaterCoeff_File), &
                       Quiet=.TRUE.)
 if ( err_stat /= SUCCESS ) THEN
   message = 'Error initializing CRTM'
   call Display_Message( PROGRAM_NAME, message, FAILURE )
   stop
 end if


 ! Loop over all sensors. Not necessary if we're calling CRTM for each sensor
 ! ----------------------------------------------------------------------------
 Sensor_Loop:do n = 1, self%conf%n_Sensors


   ! Pass channel list to CRTM
   ! -------------------------
   err_stat = CRTM_ChannelInfo_Subset(chinfo(n), self%channels, reset=.false.)
   if ( err_stat /= SUCCESS ) THEN
      message = 'Error subsetting channels!'
      call Display_Message( PROGRAM_NAME, message, FAILURE )
      stop
   end if


   ! Determine the number of channels for the current sensor
   ! -------------------------------------------------------
   n_Channels = CRTM_ChannelInfo_n_Channels(chinfo(n))

   ! Allocate the ARRAYS (for CRTM_Forward)
   ! --------------------------------------
   allocate( geo( n_Profiles ),               &
             atm( n_Profiles ),               &
             sfc( n_Profiles ),               &
             rts( n_Channels, n_Profiles ),   &
             Options( n_Profiles ),           &
             STAT = alloc_stat )
   if ( alloc_stat /= 0 ) THEN
      message = 'Error allocating structure arrays'
      call Display_Message( PROGRAM_NAME, message, FAILURE )
      stop
   end if

   if (n_Layers > 0) call CRTM_RTSolution_Create (rts, n_Layers) 

   ! Create the input FORWARD structure (atm)
   ! ----------------------------------------
   call CRTM_Atmosphere_Create( atm, n_Layers, self%conf%n_Absorbers, self%conf%n_Clouds, self%conf%n_Aerosols )
   if ( ANY(.NOT. CRTM_Atmosphere_Associated(atm)) ) THEN
      message = 'Error allocating CRTM Forward Atmosphere structure'
      CALL Display_Message( PROGRAM_NAME, message, FAILURE )
      STOP
   END IF


   ! Create the input FORWARD structure (sfc)
   ! ----------------------------------------
   call CRTM_Surface_Create(sfc, n_Channels)
   IF ( ANY(.NOT. CRTM_Surface_Associated(sfc)) ) THEN
      message = 'Error allocating CRTM Surface structure'
      CALL Display_Message( PROGRAM_NAME, message, FAILURE )
      STOP
   END IF

   CALL CRTM_RTSolution_Create(rts, n_Layers )

   !Assign the data from the GeoVaLs
   !--------------------------------
   call Load_Atm_Data(n_Profiles,n_Layers,geovals,atm,self%conf)
   call Load_Sfc_Data(n_Profiles,n_Channels,self%channels,geovals,sfc,chinfo,obss,self%conf)
   call Load_Geom_Data(obss,geo)

   ! Call THE CRTM inspection
   ! ------------------------
   if (self%conf%inspect > 0) then
     call CRTM_Atmosphere_Inspect(atm(self%conf%inspect))
     call CRTM_Surface_Inspect(sfc(self%conf%inspect))
     call CRTM_Geometry_Inspect(geo(self%conf%inspect))
     call CRTM_ChannelInfo_Inspect(chinfo(n))
   endif

   !! Parse hofxdiags%variables into independent/dependent variables and channel
   !! assumed formats:
   !!   jacobian var -->     <ystr>_jacobian_<xstr>_<chstr>
   !!   non-jacobian var --> <ystr>_<chstr>

   jacobian_needed = .false.
   do jvar = 1, hofxdiags%nvar
      varstr = hofxdiags%variables(jvar)
      str_pos(4) = len_trim(varstr)
      if (str_pos(4) < 1) cycle
      str_pos(3) = index(varstr,"_",back=.true.)        !final "_" before channel
      read(varstr(str_pos(3)+1:str_pos(4)),*) ch_diags(jvar)
      str_pos(1) = index(varstr,jacobianstr) - 1        !position before jacobianstr
      if (str_pos(1) == 0) then
         write(err_msg,*) 'ufo_radiancecrtm_simobs: _jacobian_ must be // &
                           & preceded by dependent variable in config: ', &
                           & hofxdiags%variables(jvar)
         call abor1_ftn(err_msg)
      else if (str_pos(1) > 0) then
         !Diagnostic is a Jacobian member (dy/dx)
         ystr_diags(jvar) = varstr(1:str_pos(1))
         str_pos(2) = str_pos(1) + len(jacobianstr) + 1 !begin xstr_diags
         jacobian_needed = .true.
         str_pos(4) = str_pos(3) - str_pos(2)
         xstr_diags(jvar)(1:str_pos(4)) = varstr(str_pos(2):str_pos(3)-1)
         xstr_diags(jvar)(str_pos(4)+1:) = ""
      else !null
         !Diagnostic is a dependent variable (y)
         xstr_diags(jvar) = ""
         ystr_diags(jvar)(1:str_pos(3)-1) = varstr(1:str_pos(3)-1)
         ystr_diags(jvar)(str_pos(3):) = ""
      end if 
   end do

   if (jacobian_needed) then
      ! Allocate the ARRAYS (for CRTM_K_Matrix)
      ! --------------------------------------
      allocate( atm_K( n_Channels, n_Profiles ),               &
                sfc_K( n_Channels, n_Profiles ),   &
                rts_K( n_Channels, n_Profiles ),   &
                STAT = alloc_stat )

      if ( alloc_stat /= 0 ) THEN
         message = 'Error allocating K structure arrays'
         call Display_Message( PROGRAM_NAME, message, FAILURE )
         stop
      end if

      ! Create output K-MATRIX structure (atm)
      ! --------------------------------------
      call CRTM_Atmosphere_Create( atm_K, n_Layers, self%conf%n_Absorbers, self%conf%n_Clouds, self%conf%n_Aerosols )
      if ( ANY(.NOT. CRTM_Atmosphere_Associated(atm_K)) ) THEN
         message = 'Error allocating CRTM K-matrix Atmosphere structure (setTraj)'
         CALL Display_Message( PROGRAM_NAME, message, FAILURE )
         STOP
      END IF

      ! Create output K-MATRIX structure (sfc)
      ! --------------------------------------
      call CRTM_Surface_Create( sfc_K, n_Channels)
      IF ( ANY(.NOT. CRTM_Surface_Associated(sfc_K)) ) THEN
         message = 'Error allocating CRTM K-matrix Surface structure (setTraj)'
         CALL Display_Message( PROGRAM_NAME, message, FAILURE )
         STOP
      END IF

      ! Zero the K-matrix OUTPUT structures
      ! -----------------------------------
      call CRTM_Atmosphere_Zero( atm_K )
      call CRTM_Surface_Zero( sfc_K )

      ! Inintialize the K-matrix INPUT so that the results are dTb/dx
      ! -------------------------------------------------------------
      rts_K%Radiance               = ZERO
      rts_K%Brightness_Temperature = ONE


      ! Call the K-matrix model
      ! -----------------------
      err_stat = CRTM_K_Matrix( atm         , &  ! FORWARD  Input
                                sfc         , &  ! FORWARD  Input
                                rts_K       , &  ! K-MATRIX Input
                                geo         , &  ! Input
                                chinfo(n:n) , &  ! Input
                                atm_K  , &  ! K-MATRIX Output
                                sfc_K  , &  ! K-MATRIX Output
                                rts           )  ! FORWARD  Output
      if ( err_stat /= SUCCESS ) THEN
         message = 'Error calling CRTM (setTraj) K-Matrix Model for '//TRIM(self%conf%SENSOR_ID(n))
         call Display_Message( PROGRAM_NAME, message, FAILURE )
         stop
      end if

   else

      ! Call the forward model call for each sensor
      ! -------------------------------------------
      if (self%conf%salinity_option == "sss") THEN
      	 Options%Use_Old_MWSSEM = .TRUE.
	 err_stat = CRTM_Forward( atm        , &  ! Input
                                  sfc        , &  ! Input
                                  geo        , &  ! Input
                                  chinfo(n:n), &  ! Input
                                  rts        , &  ! Output
                                  Options     )   ! Optional input
      else
         err_stat = CRTM_Forward( atm        , &  ! Input
                                  sfc        , &  ! Input
                                  geo        , &  ! Input
                                  chinfo(n:n), &  ! Input
                                  rts         )   ! Output
      end if

      if ( err_stat /= SUCCESS ) THEN
         message = 'Error calling CRTM Forward Model for '//TRIM(self%conf%SENSOR_ID(n))
         call Display_Message( PROGRAM_NAME, message, FAILURE )
         stop
      end if

   end if ! jacobian_needed

   !call CRTM_RTSolution_Inspect(rts)

   ! Put simulated brightness temperature into hofx
   ! ----------------------------------------------

   !Set to zero and initialize counter
   hofx(:,:) = 0.0_kind_real

   do m = 1, n_Profiles
     do l = 1, size(self%channels)
       hofx(l,m) = rts(l,m)%Brightness_Temperature
     end do
   end do

   ! Put simulated diagnostics into hofxdiags
   ! ----------------------------------------------
   do jvar = 1, hofxdiags%nvar
      if (len(trim(hofxdiags%variables(jvar))) < 1) cycle

      if (size(pack(self%channels,self%channels==ch_diags(jvar))) /= 1) then
         write(err_msg,*) 'ufo_radiancecrtm_simobs: mismatch between// &
                           & h(x) channels(', self%channels,') and// &
                           & ch_diags(jvar) = ', ch_diags(jvar)
         call abor1_ftn(err_msg)
      end if

      jchannel = -1
      do ichannel = 1, size(self%channels)
         if (ch_diags(jvar) == self%channels(ichannel)) then
            jchannel = ichannel
            exit
         end if
      end do

      if (allocated(hofxdiags%geovals(jvar)%vals)) &
         deallocate(hofxdiags%geovals(jvar)%vals)

      !=========================
      ! Diagnostics used for QC
      !=========================
      if (trim(xstr_diags(jvar)) == "") then
         ! forward h(x) diags
         select case(ystr_diags(jvar))
            ! variable: optical_thickness_of_atmosphere_layer_CH
            case (var_opt_depth)
               hofxdiags%geovals(jvar)%nval = n_Layers
               allocate(hofxdiags%geovals(jvar)%vals(hofxdiags%geovals(jvar)%nval,n_Profiles))
               do jprofile = 1, n_Profiles
                  do jlevel = 1, hofxdiags%geovals(jvar)%nval
                     hofxdiags%geovals(jvar)%vals(jlevel,jprofile) = &
                        rts(jchannel,jprofile) % layer_optical_depth(jlevel)
                  end do
               end do

            ! variable: toa_outgoing_radiance_per_unit_wavenumber_CH [mW / (m^2 sr cm^-1)] (nval=1)
            case (var_radiance)
               hofxdiags%geovals(jvar)%nval = 1
               allocate(hofxdiags%geovals(jvar)%vals(hofxdiags%geovals(jvar)%nval,n_Profiles))
               do jprofile = 1, n_Profiles
                  hofxdiags%geovals(jvar)%vals(1,jprofile) = &
                     rts(jchannel,jprofile) % Radiance
               end do

            ! variable: brightness_temperature_assuming_clear_sky_CH
            case (var_tb_clr)
               hofxdiags%geovals(jvar)%nval = 1
               allocate(hofxdiags%geovals(jvar)%vals(hofxdiags%geovals(jvar)%nval,n_Profiles))
               do jprofile = 1, n_Profiles
                  ! Note: Using Tb_Clear requires CRTM_Atmosphere_IsFractional(cloud_coverage_flag) 
                  ! to be true. For CRTM v2.3.0, that happens when 
                  ! atm(jprofile)%Cloud_Fraction > MIN_COVERAGE_THRESHOLD (1e.-6)
                  hofxdiags%geovals(jvar)%vals(1,jprofile) = &
                     rts(jchannel,jprofile) % Tb_Clear 
               end do
            ! variable: brightness_temperature_CH
            case (var_tb)
               hofxdiags%geovals(jvar)%nval = 1
               allocate(hofxdiags%geovals(jvar)%vals(hofxdiags%geovals(jvar)%nval,n_Profiles))
               do jprofile = 1, n_Profiles
                  hofxdiags%geovals(jvar)%vals(1,jprofile) = &
                     rts(jchannel,jprofile) % Brightness_Temperature 
               end do
            case default
               write(err_msg,*) 'ufo_radiancecrtm_simobs: //&
                                 & ObsDiagnostic is unsupported, ', &
                                 & hofxdiags%variables(jvar)
               call abor1_ftn(err_msg)
         end select
      else if (ystr_diags(jvar) == var_tb) then
         ! var_tb jacobians
         select case (xstr_diags(jvar))
            ! variable: brightness_temperature_jacobian_air_temperature_CH
            case (var_ts)
               hofxdiags%geovals(jvar)%nval = n_Layers
               allocate(hofxdiags%geovals(jvar)%vals(hofxdiags%geovals(jvar)%nval,n_Profiles))
               do jprofile = 1, n_Profiles
                  do jlevel = 1, hofxdiags%geovals(jvar)%nval
                     hofxdiags%geovals(jvar)%vals(jlevel,jprofile) = &
                        atm_K(jchannel,jprofile) % Temperature(jlevel)
                  end do
               end do
            ! variable: brightness_temperature_jacobian_humidity_mixing_ratio_CH (nval==n_Layers) --> requires MAXVARLEN=58
            case (var_mixr)
               hofxdiags%geovals(jvar)%nval = n_Layers
               allocate(hofxdiags%geovals(jvar)%vals(hofxdiags%geovals(jvar)%nval,n_Profiles))
               jspec = ufo_vars_getindex(self%conf%Absorbers, var_mixr)
               do jprofile = 1, n_Profiles
                  do jlevel = 1, hofxdiags%geovals(jvar)%nval
                     hofxdiags%geovals(jvar)%vals(jlevel,jprofile) = &
                        atm_K(jchannel,jprofile) % Absorber(jlevel,jspec)
                  end do
               end do

            ! variable: brightness_temperature_jacobian_surface_temperature_CH (nval=1)
            case (var_sfc_t)
               hofxdiags%geovals(jvar)%nval = 1
               allocate(hofxdiags%geovals(jvar)%vals(hofxdiags%geovals(jvar)%nval,n_Profiles))
               do jprofile = 1, n_Profiles
                  hofxdiags%geovals(jvar)%vals(1,jprofile) = &
                     sfc_K(jchannel,jprofile) % water_temperature &
                   + sfc_K(jchannel,jprofile) % land_temperature &
                   + sfc_K(jchannel,jprofile) % ice_temperature &
                   + sfc_K(jchannel,jprofile) % snow_temperature
               end do

            ! variable: brightness_temperature_jacobian_surface_emissivity_CH (nval=1)
            case (var_sfc_emiss)
               hofxdiags%geovals(jvar)%nval = 1
               allocate(hofxdiags%geovals(jvar)%vals(hofxdiags%geovals(jvar)%nval,n_Profiles))
               do jprofile = 1, n_Profiles
                  hofxdiags%geovals(jvar)%vals(1,jprofile) = &
                     rts_K(jchannel,jprofile) % surface_emissivity
               end do
            case default
               write(err_msg,*) 'ufo_radiancecrtm_simobs: //&
                                 & ObsDiagnostic is unsupported, ', &
                                 & hofxdiags%variables(jvar)
               call abor1_ftn(err_msg)
         end select
      else
         write(err_msg,*) 'ufo_radiancecrtm_simobs: //&
                           & ObsDiagnostic is unsupported, ', &
                           & hofxdiags%variables(jvar)
         call abor1_ftn(err_msg)
      end if
   end do


   ! Deallocate the structures
   ! -------------------------
   call CRTM_Geometry_Destroy(geo)
   call CRTM_Atmosphere_Destroy(atm)
   call CRTM_Surface_Destroy(sfc)
   call CRTM_RTSolution_Destroy(rts)

   ! Deallocate all arrays
   ! ---------------------
   deallocate(geo, atm, sfc, rts, STAT = alloc_stat)
   if ( alloc_stat /= 0 ) THEN
      message = 'Error deallocating structure arrays'
      call Display_Message( PROGRAM_NAME, message, FAILURE )
      stop
   end if

   if (jacobian_needed) then
      ! Deallocate the K structures
      ! ---------------------------
      call CRTM_Atmosphere_Destroy(atm_K)
      call CRTM_Surface_Destroy(sfc_K)
      call CRTM_RTSolution_Destroy(rts_K)

      ! Deallocate all K arrays
      ! -----------------------
      deallocate(atm_K, sfc_K, rts_K, STAT = alloc_stat)
      if ( alloc_stat /= 0 ) THEN
         message = 'Error deallocating K structure arrays'
         call Display_Message( PROGRAM_NAME, message, FAILURE )
         stop
      end if
   end if

 end do Sensor_Loop


 ! Destroy CRTM instance
 ! ---------------------
 write( *, '( /5x, "Destroying the CRTM..." )' )
 err_stat = CRTM_Destroy( chinfo )
 if ( err_stat /= SUCCESS ) THEN
    message = 'Error destroying CRTM'
    call Display_Message( PROGRAM_NAME, message, FAILURE )
    stop
 end if

end subroutine ufo_radiancecrtm_simobs

! ------------------------------------------------------------------------------

end module ufo_radiancecrtm_mod
