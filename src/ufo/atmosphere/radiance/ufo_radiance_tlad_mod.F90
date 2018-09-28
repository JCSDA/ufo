! (C) Copyright 2017-2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Fortran module to handle tl/ad for radiance observations

module ufo_radiance_tlad_mod
  
 use iso_c_binding
 use config_mod
 use kinds

 use ioda_obsdb_mod, only: ioda_obsdb
 use ioda_obs_vectors, only: obs_vector

 use ufo_geovals_mod, only: ufo_geovals
 use ufo_basis_tlad_mod, only: ufo_basis_tlad
 use ufo_vars_mod

 !YT hack needs
 use ufo_geovals_mod, only: ufo_geovals_init, ufo_geovals_zero, ufo_geovals_copy, ufo_geovals_setup

 use ufo_radiance_utils_mod

 use crtm_module

 implicit none
 private

 real(kind=kind_real), parameter :: rmissing = -9.9e10_kind_real

 !> Fortran derived type for radiance trajectory
 type, extends(ufo_basis_tlad), public :: ufo_radiance_tlad
 private
  type(rad_conf) :: rc
  integer :: n_Profiles
  integer :: n_Channels
  type(ufo_geovals) :: crtm_K
  type(ufo_geovals) :: geohack
 contains
  procedure :: setup  => ufo_radiance_tlad_setup
  procedure :: delete  => ufo_radiance_tlad_delete
  procedure :: settraj => ufo_radiance_tlad_settraj 
  procedure :: simobs_tl  => ufo_radiance_simobs_tl
  procedure :: simobs_ad  => ufo_radiance_simobs_ad
 end type ufo_radiance_tlad

contains

! ------------------------------------------------------------------------------

subroutine ufo_radiance_tlad_setup(self, c_conf)

implicit none
class(ufo_radiance_tlad), intent(inout) :: self
type(c_ptr),              intent(in)    :: c_conf

 call rad_conf_setup(self%rc,c_conf)

end subroutine ufo_radiance_tlad_setup

! ------------------------------------------------------------------------------

subroutine ufo_radiance_tlad_delete(self)

implicit none
class(ufo_radiance_tlad), intent(inout) :: self

 self%ltraj = .false.
 call rad_conf_delete(self%rc)

end subroutine ufo_radiance_tlad_delete
  
! ------------------------------------------------------------------------------
  
subroutine ufo_radiance_tlad_settraj(self, geovals, obss)

implicit none

class(ufo_radiance_tlad), intent(inout) :: self
type(ufo_geovals),        intent(in)    :: geovals
type(ioda_obsdb),         intent(in)    :: obss

! Local Variables
character(*), parameter :: PROGRAM_NAME = 'ufo_radiance_mod.F90'
character(255) :: message, version
integer        :: err_stat, alloc_stat
integer        :: n, k1

! Define the "non-demoninational" arguments
type(CRTM_ChannelInfo_type)             :: chinfo(self%rc%n_Sensors)
type(CRTM_Geometry_type),   allocatable :: geo(:)

! Define the FORWARD variables
type(CRTM_Atmosphere_type), allocatable :: atm(:)
type(CRTM_Surface_type),    allocatable :: sfc(:)
type(CRTM_RTSolution_type), allocatable :: rts(:,:)

! Define the K-MATRIX variables
type(CRTM_Atmosphere_type), allocatable :: atm_K(:,:)
type(CRTM_Surface_type)   , allocatable :: sfc_K(:,:)
type(CRTM_RTSolution_type), allocatable :: rts_K(:,:)


 ! Get number of profile and layers from geovals
 ! ---------------------------------------------
 self%n_Profiles = geovals%nobs


 ! Program header
 ! --------------
 call CRTM_Version( Version )
 call Program_Message( PROGRAM_NAME, &
                       'Check/example program for the CRTM Forward and K-Matrix (setTraj) functions using '//&
                       trim(self%rc%ENDIAN_type)//' coefficient datafiles', &
                       'CRTM Version: '//TRIM(Version) )    
    

 ! Initialise all the sensors at once
 ! ----------------------------------
 !** NOTE: CRTM_Init points to the various binary files needed for CRTM.  See the
 !**       CRTM_Lifecycle.f90 for more details. 

 write( *,'(/5x,"Initializing the CRTM (setTraj) ...")' )
 err_stat = CRTM_Init( self%rc%SENSOR_ID, &
            chinfo, &
            File_Path=trim(self%rc%COEFFICIENT_PATH), &
            Quiet=.TRUE.)
 if ( err_stat /= SUCCESS ) THEN
   message = 'Error initializing CRTM (setTraj)'
   call Display_Message( PROGRAM_NAME, message, FAILURE )
   stop
 end if

 ! Loop over all sensors. Not necessary if we're calling CRTM for each sensor
 ! ----------------------------------------------------------------------------
 Sensor_Loop:do n = 1, self%rc%n_Sensors


   ! Determine the number of channels for the current sensor
   ! -------------------------------------------------------
   self%N_Channels = CRTM_ChannelInfo_n_Channels(chinfo(n))       


   ! Allocate the ARRAYS
   ! -------------------
   allocate( geo( self%n_Profiles ),               &
             atm( self%n_Profiles ),               &
             sfc( self%n_Profiles ),               &
             rts( self%N_Channels, self%n_Profiles ),   &
             atm_K( self%N_Channels, self%n_Profiles ), &
             sfc_K( self%N_Channels, self%n_Profiles ), &
             rts_K( self%N_Channels, self%n_Profiles ), &
             STAT = alloc_stat )
   if ( alloc_stat /= 0 ) THEN
      message = 'Error allocating structure arrays (setTraj)'
      call Display_Message( PROGRAM_NAME, message, FAILURE )
      stop
   end if


   ! Create the input FORWARD structure (atm)
   ! ----------------------------------------
   call CRTM_Atmosphere_Create( atm, self%rc%n_Layers, self%rc%n_Absorbers, self%rc%n_Clouds, self%rc%n_Aerosols )
   if ( ANY(.NOT. CRTM_Atmosphere_Associated(atm)) ) THEN
      message = 'Error allocating CRTM Forward Atmosphere structure (setTraj)'
      CALL Display_Message( PROGRAM_NAME, message, FAILURE )
      STOP
   END IF

       
   ! Create the input FORWARD structure (sfc)
   ! ----------------------------------------
   call CRTM_Surface_Create(sfc, self%N_Channels)
   IF ( ANY(.NOT. CRTM_Surface_Associated(sfc)) ) THEN
      message = 'Error allocating CRTM Surface structure (setTraj)'
      CALL Display_Message( PROGRAM_NAME, message, FAILURE )
      STOP
   END IF


   ! Create output K-MATRIX structure (atm)
   ! --------------------------------------
   call CRTM_Atmosphere_Create( atm_K, self%rc%n_Layers, self%rc%n_Absorbers, self%rc%n_Clouds, self%rc%n_Aerosols )
   if ( ANY(.NOT. CRTM_Atmosphere_Associated(atm_K)) ) THEN
      message = 'Error allocating CRTM K-matrix Atmosphere structure (setTraj)'
      CALL Display_Message( PROGRAM_NAME, message, FAILURE )
      STOP
   END IF


   ! Create output K-MATRIX structure (sfc)
   ! --------------------------------------
   call CRTM_Surface_Create(sfc_K, self%N_Channels)
   IF ( ANY(.NOT. CRTM_Surface_Associated(sfc_K)) ) THEN
      message = 'Error allocating CRTM K-matrix Surface structure (setTraj)'
      CALL Display_Message( PROGRAM_NAME, message, FAILURE )
      STOP
   END IF


   !Assign the data from the GeoVaLs
   !--------------------------------
   call Load_Atm_Data(self%N_PROFILES,self%rc%N_LAYERS,geovals,atm)
   call Load_Sfc_Data(self%N_PROFILES,self%rc%N_LAYERS,self%N_Channels,geovals,sfc,chinfo,obss)
   call Load_Geom_Data(obss,geo)

   !Hack absorbers and clouds
   do k1 = 1,self%N_PROFILES
     atm(k1)%Absorber = 0.0
     atm(k1)%Cloud(1)%Water_Content = 0.0
     atm(k1)%Cloud(2)%Water_Content = 0.0
   enddo

   ! Zero the K-matrix OUTPUT structures
   ! -----------------------------------
   call CRTM_Atmosphere_Zero( atm_K )
   call CRTM_Surface_Zero( sfc_K )


   ! Inintialize the K-matrix INPUT so that the results are dTb/dx
   ! -------------------------------------------------------------
   rts_K%Radiance               = ZERO
   rts_K%Brightness_Temperature = ONE


   ! Call the K-matrix model (really need this here or just in settraj?)
   ! ----------------------
   err_stat = CRTM_K_Matrix( atm, &  ! FORWARD  Input
        sfc                     , &  ! FORWARD  Input
        rts_K                   , &  ! K-MATRIX Input
        geo                     , &  ! Input
        chinfo(n:n)             , &  ! Input
        atm_K                   , &  ! K-MATRIX Output
        sfc_K                   , &  ! K-MATRIX Output
        rts          )               ! FORWARD  Output
   if ( err_stat /= SUCCESS ) THEN
      message = 'Error calling CRTM (setTraj) K-Matrix Model for '//TRIM(self%rc%SENSOR_ID(n))
      call Display_Message( PROGRAM_NAME, message, FAILURE )
      stop
   end if


   ! Populate the CRTM K matrix
   ! --------------------------
   call populate_crtm_K(self,geovals,atm_K,sfc_K)


   ! Deallocate the structures
   ! -------------------------
   call CRTM_Geometry_Destroy(geo)
   call CRTM_Atmosphere_Destroy(atm_K)
   call CRTM_Atmosphere_Destroy(atm)
   call CRTM_RTSolution_Destroy(rts_K)
   call CRTM_RTSolution_Destroy(rts)
   call CRTM_Surface_Destroy(sfc)
   call CRTM_Surface_Destroy(sfc_K)


   ! Deallocate all arrays
   ! ---------------------
   deallocate(geo, atm, sfc, rts, rts_K, sfc_K, atm_K, STAT = alloc_stat)
   if ( alloc_stat /= 0 ) THEN
      message = 'Error deallocating structure arrays (setTraj)'
      call Display_Message( PROGRAM_NAME, message, FAILURE )
      stop
   end if

 end do Sensor_Loop


 ! Destroy CRTM instance
 ! ---------------------
 write( *, '( /5x, "Destroying the CRTM (setTraj)..." )' )
 err_stat = CRTM_Destroy( chinfo )
 if ( err_stat /= SUCCESS ) THEN
    message = 'Error destroying CRTM (setTraj)'
    call Display_Message( PROGRAM_NAME, message, FAILURE )
    stop
 end if

 
 ! Hack | Save a copy of the geovals for the tl/ad 
 call ufo_geovals_init(self%geohack)
 call ufo_geovals_copy(geovals, self%geohack)
 call ufo_geovals_zero(self%geohack)


end subroutine ufo_radiance_tlad_settraj

! ------------------------------------------------------------------------------

subroutine populate_crtm_K(self,geovals,atm_k,sfc_k)

implicit none

class(ufo_radiance_tlad),   intent(inout) :: self
type(ufo_geovals),          intent(in)    :: geovals
type(CRTM_Atmosphere_type), intent(inout) :: atm_K(:,:)
type(CRTM_Surface_type),    intent(inout) :: sfc_K(:,:)

! Local variables
integer :: k1, k2, k3, ivar

 !** NOTES: this is to populate the jacobian structures, using exactly the same structures as geovals, including variable names.
 !**        this is the laziest possible way to do this, and will likely need to be changed.
 !**        use public subroutine ufo_geovals_assign to initialize crtm_K with the exact same values geovals?

 call ufo_geovals_init(self%crtm_K)
 call ufo_geovals_setup(self%crtm_K, geovals%variables, self%n_Profiles*self%n_Channels) !** setup jacobian structure using geovals structure.

 !** atmosphere
 ivar = ufo_vars_getindex(geovals%variables, var_tv  )
 self%crtm_K%geovals(ivar)%nval = self%rc%n_Layers
 ivar = ufo_vars_getindex(geovals%variables, var_prs )
 self%crtm_K%geovals(ivar)%nval = self%rc%n_Layers
 ivar = ufo_vars_getindex(geovals%variables, var_prsi)
 self%crtm_K%geovals(ivar)%nval = self%rc%n_Layers+1
 ivar = ufo_vars_getindex(geovals%variables, var_mixr)
 self%crtm_K%geovals(ivar)%nval = self%rc%n_Layers
 ivar = ufo_vars_getindex(geovals%variables, var_oz  )
 self%crtm_K%geovals(ivar)%nval = self%rc%n_Layers
 ivar = ufo_vars_getindex(geovals%variables, var_co2 )
 self%crtm_K%geovals(ivar)%nval = self%rc%n_Layers

 !** cloud 1 
 ivar = ufo_vars_getindex(geovals%variables, var_clw )
 self%crtm_K%geovals(ivar)%nval = self%rc%n_Layers
 ivar = ufo_vars_getindex(geovals%variables, var_clwefr)
 self%crtm_K%geovals(ivar)%nval = self%rc%n_Layers

 !** cloud 2
 ivar = ufo_vars_getindex(geovals%variables, var_cli )
 self%crtm_K%geovals(ivar)%nval = self%rc%n_Layers
 ivar = ufo_vars_getindex(geovals%variables, var_cliefr)
 self%crtm_K%geovals(ivar)%nval = self%rc%n_Layers

 !** surface
 ivar = ufo_vars_getindex(geovals%variables, var_sfc_wspeed )
 self%crtm_K%geovals(ivar)%nval = 1
 ivar = ufo_vars_getindex(geovals%variables, var_sfc_wdir   )
 self%crtm_K%geovals(ivar)%nval = 1
 ivar = ufo_vars_getindex(geovals%variables, var_sfc_wfrac  )
 self%crtm_K%geovals(ivar)%nval = 1
 ivar = ufo_vars_getindex(geovals%variables, var_sfc_wtmp   )
 self%crtm_K%geovals(ivar)%nval = 1
 ivar = ufo_vars_getindex(geovals%variables, var_sfc_ifrac  )
 self%crtm_K%geovals(ivar)%nval = 1
 ivar = ufo_vars_getindex(geovals%variables, var_sfc_itmp   )
 self%crtm_K%geovals(ivar)%nval = 1
 ivar = ufo_vars_getindex(geovals%variables, var_sfc_sfrac  )
 self%crtm_K%geovals(ivar)%nval = 1
 ivar = ufo_vars_getindex(geovals%variables, var_sfc_stmp   )
 self%crtm_K%geovals(ivar)%nval = 1
 ivar = ufo_vars_getindex(geovals%variables, var_sfc_sdepth )
 self%crtm_K%geovals(ivar)%nval = 1
 ivar = ufo_vars_getindex(geovals%variables, var_sfc_landtyp)
 self%crtm_K%geovals(ivar)%nval = 1
 ivar = ufo_vars_getindex(geovals%variables, var_sfc_lfrac  )
 self%crtm_K%geovals(ivar)%nval = 1
 ivar = ufo_vars_getindex(geovals%variables, var_sfc_ltmp   )
 self%crtm_K%geovals(ivar)%nval = 1
 ivar = ufo_vars_getindex(geovals%variables, var_sfc_lai    )
 self%crtm_K%geovals(ivar)%nval = 1
 ivar = ufo_vars_getindex(geovals%variables, var_sfc_vegfrac)
 self%crtm_K%geovals(ivar)%nval = 1
 ivar = ufo_vars_getindex(geovals%variables, var_sfc_vegtyp )
 self%crtm_K%geovals(ivar)%nval = 1
 ivar = ufo_vars_getindex(geovals%variables, var_sfc_soiltyp)
 self%crtm_K%geovals(ivar)%nval = 1
 ivar = ufo_vars_getindex(geovals%variables, var_sfc_soilm  )
 self%crtm_K%geovals(ivar)%nval = 1
 ivar = ufo_vars_getindex(geovals%variables, var_sfc_soilt  )
 self%crtm_K%geovals(ivar)%nval = 1

 do k3 = 1,geovals%nvar
    allocate(self%crtm_K%geovals(k3)%vals(self%crtm_K%geovals(k3)%nval,self%n_Profiles*self%n_Channels))
    self%crtm_K%geovals(k3)%vals(:,:) = 0.0_fp
 end do
 self%crtm_K%linit = .true.

 !** atm_K and sfc_K contain the jacobian structures, and are populated prior to this subroutine being called.
 !** next step is to copy the values into the self%crtm_K structure, with appropriate unit conversions.
 !** the purpose is to have an exact mapping between geovals and the jacobians.  I will eventually need an automated
 !** method to copy the values from atm_K and sfc_K into self%crtm_K.

 !** populate the atmosphere K-matrix (jacobian) structures for CRTM (self%crtm_K(k1), for the k1-th profile)
 !!$      var_tv     atm%Temperature
 !!$      var_prs    atm%Pressure
 !!$      var_prsi   atm%Level_Pressure
 !!$      var_mixr   atm%Absorber(1) !** verify
 !!$      var_oz     atm%Absorber(2) !** verify
 !!$      var_co2    atm%Absorber(3) !** verify
 !!$      var_clw    atm%Cloud(1)%Water_Content !** water cloud content
 !!$      var_clwefr atm%Cloud(1)%Effective_Radius !** water cloud effective radius 
 !!$      var_cli    atm%Cloud(2)%Water_Content !** ice cloud content
 !!$      var_cliefr atm%Cloud(2)%Effective_Radius !** ice cloud effective radius 
 
 !!$      var_sfc_wspeed     sfc%Wind_Speed
 !!$      var_sfc_wdir       sfc%Wind_Direction
 !!$      var_sfc_wfrac      sfc%Water_Coverage
 !!$      var_sfc_wtmp       sfc%Water_Temperature
 !!$      var_sfc_ifrac      sfc%Ice_Coverage 
 !!$      var_sfc_itmp       sfc%Ice_Temperature
 !!$      var_sfc_sfrac      sfc%Snow_Coverage
 !!$      var_sfc_stmp       sfc%Snow_Temperature
 !!$      var_sfc_sdepth     sfc%Snow_Depth
 !!$      var_sfc_landtyp    sfc%Land_Type
 !!$      var_sfc_lfrac      sfc%Land_Coverage 
 !!$      var_sfc_ltmp       sfc%Land_Temperature
 !!$      var_sfc_lai        sfc%Lai       
 !!$      var_sfc_vegfrac    sfc%Vegetation_Fraction
 !!$      var_sfc_vegtyp     sfc%Vegetation_Type
 !!$      var_sfc_soiltyp    sfc%Soil_Type
 !!$      var_sfc_soilm      sfc%Soil_Moisture
 !!$      var_sfc_soilt      sfc%Soil_Temperature

 k3 = 0

 do k1 = 1,self%n_Profiles
   do k2 = 1,self%n_Channels
     k3 = k3 + 1  !** jacobian is self%n_profiles, self%n_channels, geoval is self%rc%n_layers, self%n_obs.  
                  !** k3 flattens self%n_profiles,self%n_channels.
     !** atmosphere
     ivar = ufo_vars_getindex(geovals%variables, var_tv  )
     self%crtm_K%geovals(ivar)%vals(1:self%rc%n_Layers,k3)   = atm_K(k2,k1)%Temperature(1:self%rc%n_Layers)
     ivar = ufo_vars_getindex(geovals%variables, var_prs )
     self%crtm_K%geovals(ivar)%vals(1:self%rc%n_Layers,k3)   = 0.0 !atm_K(k2,k1)%Pressure(1:self%rc%n_Layers) 
     ivar = ufo_vars_getindex(geovals%variables, var_prsi)
     self%crtm_K%geovals(ivar)%vals(1:self%rc%n_Layers+1,k3)   = 0.0 !atm_K(k2,k1)%Level_Pressure(0:self%rc%n_Layers) 
     ivar = ufo_vars_getindex(geovals%variables, var_mixr)
     self%crtm_K%geovals(ivar)%vals(1:self%rc%n_Layers,k3)   = 0.0 !atm_K(k2,k1)%Absorber(1:self%rc%n_Layers,1)
     ivar = ufo_vars_getindex(geovals%variables, var_oz  )
     self%crtm_K%geovals(ivar)%vals(1:self%rc%n_Layers,k3)   = 0.0 !atm_K(k2,k1)%Absorber(1:self%rc%n_Layers,2)
     ivar = ufo_vars_getindex(geovals%variables, var_co2 )
     self%crtm_K%geovals(ivar)%vals(1:self%rc%n_Layers,k3)   = 0.0 !atm_K(k2,k1)%Absorber(1:self%rc%n_Layers,3)

     !** cloud 1 
     ivar = ufo_vars_getindex(geovals%variables, var_clw )
     self%crtm_K%geovals(ivar)%vals(1:self%rc%n_Layers,k3)   = 0.0 !atm_K(k2,k1)%Cloud(1)%Water_Content(1:self%rc%n_Layers)
     ivar = ufo_vars_getindex(geovals%variables, var_clwefr)
     self%crtm_K%geovals(ivar)%vals(1:self%rc%n_Layers,k3)   = 0.0 !atm_K(k2,k1)%Cloud(1)%Effective_Radius(1:self%rc%n_Layers)

     !** cloud 2
     ivar = ufo_vars_getindex(geovals%variables, var_cli )
     self%crtm_K%geovals(ivar)%vals(1:self%rc%n_Layers,k3)   = 0.0 !atm_K(k2,k1)%Cloud(2)%Water_Content(1:self%rc%n_Layers)
     ivar = ufo_vars_getindex(geovals%variables, var_cliefr)
     self%crtm_K%geovals(ivar)%vals(1:self%rc%n_Layers,k3)   = 0.0 !atm_K(k2,k1)%Cloud(2)%Effective_Radius(1:self%rc%n_Layers) 

     !** surface
     ivar = ufo_vars_getindex(geovals%variables, var_sfc_wspeed )
     self%crtm_K%geovals(ivar)%vals(1,k3)      = 0.0 !sfc_K(k2,k1)%Wind_Speed        
     ivar = ufo_vars_getindex(geovals%variables, var_sfc_wdir   )
     self%crtm_K%geovals(ivar)%vals(1,k3)      = 0.0 !sfc_K(k2,k1)%Wind_Direction    
     ivar = ufo_vars_getindex(geovals%variables, var_sfc_wfrac  )
     self%crtm_K%geovals(ivar)%vals(1,k3)      = 0.0 !sfc_K(k2,k1)%Water_Coverage    
     ivar = ufo_vars_getindex(geovals%variables, var_sfc_wtmp   )
     self%crtm_K%geovals(ivar)%vals(1,k3)      = 0.0 !sfc_K(k2,k1)%Water_Temperature 
     ivar = ufo_vars_getindex(geovals%variables, var_sfc_ifrac  )
     self%crtm_K%geovals(ivar)%vals(1,k3)      = 0.0 !sfc_K(k2,k1)%Ice_Coverage      
     ivar = ufo_vars_getindex(geovals%variables, var_sfc_itmp   )
     self%crtm_K%geovals(ivar)%vals(1,k3)      = 0.0 !sfc_K(k2,k1)%Ice_Temperature   
     ivar = ufo_vars_getindex(geovals%variables, var_sfc_sfrac  )
     self%crtm_K%geovals(ivar)%vals(1,k3)      = 0.0 !sfc_K(k2,k1)%Snow_Coverage     
     ivar = ufo_vars_getindex(geovals%variables, var_sfc_stmp   )
     self%crtm_K%geovals(ivar)%vals(1,k3)      = 0.0 !sfc_K(k2,k1)%Snow_Temperature 
     ivar = ufo_vars_getindex(geovals%variables, var_sfc_sdepth )
     self%crtm_K%geovals(ivar)%vals(1,k3)      = 0.0 !sfc_K(k2,k1)%Snow_Depth       
     ivar = ufo_vars_getindex(geovals%variables, var_sfc_landtyp)
     self%crtm_K%geovals(ivar)%vals(1,k3)      = 0.0 !sfc_K(k2,k1)%Land_Type        
     ivar = ufo_vars_getindex(geovals%variables, var_sfc_lfrac  )
     self%crtm_K%geovals(ivar)%vals(1,k3)      = 0.0 !sfc_K(k2,k1)%Land_Coverage    
     ivar = ufo_vars_getindex(geovals%variables, var_sfc_ltmp   )
     self%crtm_K%geovals(ivar)%vals(1,k3)      = 0.0 !sfc_K(k2,k1)%Land_Temperature  
     ivar = ufo_vars_getindex(geovals%variables, var_sfc_lai    )
     self%crtm_K%geovals(ivar)%vals(1,k3)      = 0.0 !sfc_K(k2,k1)%Lai               
     ivar = ufo_vars_getindex(geovals%variables, var_sfc_vegfrac)
     self%crtm_K%geovals(ivar)%vals(1,k3)      = 0.0 !sfc_K(k2,k1)%Vegetation_Fraction
     ivar = ufo_vars_getindex(geovals%variables, var_sfc_vegtyp )
     self%crtm_K%geovals(ivar)%vals(1,k3)      = 0.0 !sfc_K(k2,k1)%Vegetation_Type   
     ivar = ufo_vars_getindex(geovals%variables, var_sfc_soiltyp)
     self%crtm_K%geovals(ivar)%vals(1,k3)      = 0.0 !sfc_K(k2,k1)%Soil_Type        
     ivar = ufo_vars_getindex(geovals%variables, var_sfc_soilm  )
     self%crtm_K%geovals(ivar)%vals(1,k3)      = 0.0 !sfc_K(k2,k1)%Soil_Moisture_Content
     ivar = ufo_vars_getindex(geovals%variables, var_sfc_soilt  )
     self%crtm_K%geovals(ivar)%vals(1,k3)      = 0.0 !sfc_K(k2,k1)%Soil_Temperature
   end do
 end do

 self%ltraj = .true.  !** set trajectory is true
 !print '(A,I5,2G12.4)', 'xjacobian1:', k3, maxval(abs(self%crtm_K%geovals(1)%vals(:,k3))), maxval(abs(geovals%geovals(1)%vals(:,k3)))

end subroutine populate_crtm_k
    
! ------------------------------------------------------------------------------
  
subroutine ufo_radiance_simobs_tl(self, geovals, hofx, obss)

implicit none
class(ufo_radiance_tlad), intent(in)  :: self
type(ufo_geovals),     intent(in)     :: geovals
type(obs_vector),      intent(inout)  :: hofx
type(ioda_obsdb),      intent(in)     :: obss

character(len=*), parameter :: myname_="ufo_radiance_simobs_tl"
character(len=MAXVARLEN) :: fldname
integer jvar, jobs, jprofile, jchannel, jlev, ivar

 if (geovals%nvar /= self%crtm_K%nvar) call abor1_ftn("radiance_tl: error nvar")
 if (hofx%nobs /= self%n_Profiles*self%n_Channels) call abor1_ftn("radiance_tl: obsvector wrong size")

 hofx%values(:) = 0.0

 do jvar = 1, geovals%nvar
 fldname=geovals%variables%fldnames(jvar)
 if (trim(fldname)=="temperature" .or. trim(fldname)=="virtual_temperature" .or. &
   & trim(fldname)=="humidity_mixing_ratio") then
   ivar = ufo_vars_getindex(self%crtm_K%variables, geovals%variables%fldnames(jvar))
   write(*,*)'radiance_tl geovals nvar fldname ', geovals%nvar, geovals%variables%fldnames(jvar)
   write(*,*)'radiance_tl crtm_K ivar fldname ', ivar, self%crtm_K%variables%fldnames(ivar)
   write(*,*)'radiance_tl nval ', geovals%geovals(jvar)%nval, self%crtm_K%geovals(ivar)%nval
   if (geovals%geovals(jvar)%nval /= self%crtm_K%geovals(ivar)%nval) call abor1_ftn("radiance_tl: error nval")
   write(*,*)'radiance_tl geovals ',geovals%variables%fldnames(jvar),' min max = ', &
           & minval(geovals%geovals(jvar)%vals(:,:)),maxval(geovals%geovals(jvar)%vals(:,:))
   jobs = 0
   do jprofile = 1, self%n_Profiles
   do jchannel = 1, self%n_Channels
     jobs = jobs + 1
     do jlev = 1, geovals%geovals(jvar)%nval
       hofx%values(jobs) = hofx%values(jobs) + self%crtm_K%geovals(ivar)%vals(jlev,jobs) * geovals%geovals(jvar)%vals(jlev,jprofile)
     enddo
   enddo
   enddo
 endif
 enddo
 write(*,*)'radiance_tl hofx min max = ',minval(hofx%values(:)),maxval(hofx%values(:))

end subroutine ufo_radiance_simobs_tl
  
! ------------------------------------------------------------------------------
  
subroutine ufo_radiance_simobs_ad(self, geovals, hofx, obss)

implicit none
class(ufo_radiance_tlad), intent(in) :: self
type(ufo_geovals),     intent(inout) :: geovals
type(obs_vector),      intent(in)    :: hofx
type(ioda_obsdb),      intent(in)    :: obss

character(len=*), parameter :: myname_="ufo_radiance_simobs_ad"
integer jvar, jobs, jprofile, jchannel, jlev, ivar, jj, iavg
real(kind_real) :: zmin, zmax, zavg
character(len=MAXVARLEN) :: fldname

 call ufo_geovals_copy(self%geohack, geovals)

 if (geovals%nvar /= self%crtm_K%nvar) call abor1_ftn("radiance_ad: error nvar")
 if (hofx%nobs /= self%n_Profiles*self%n_Channels) call abor1_ftn("radiance_ad: obsvector wrong size")

 write(*,*)'radiance_ad starting'
 zmin = huge(zmin)
 zmax = -huge(zmax)
 zavg = 0.0_kind_real
 iavg = 0
 do jj=1,hofx%nobs
   if (hofx%values(jj)>rmissing) then
     zmin = min(hofx%values(jj), zmin)
     zmax = max(hofx%values(jj), zmax)
     zavg = zavg + hofx%values(jj)
     iavg = iavg + 1
   endif
 enddo
 zavg = zavg/real(iavg)
 write(*,*)'radiance_ad hofx min max avg = ',zmin,zmax,zavg

 do jvar = 1, geovals%nvar
 fldname=geovals%variables%fldnames(jvar)
 if (trim(fldname)=="temperature" .or. trim(fldname)=="virtual_temperature" .or. &
   & trim(fldname)=="humidity_mixing_ratio") then
   ivar = ufo_vars_getindex(self%crtm_K%variables, geovals%variables%fldnames(jvar))
   write(*,*)'radiance_ad geovals nvar fldname ', geovals%nvar, geovals%variables%fldnames(jvar)
   write(*,*)'radiance_ad crtm_K ivar fldname ', ivar, self%crtm_K%variables%fldnames(ivar)
   write(*,*)'radiance_ad nval ', geovals%geovals(jvar)%nval, self%crtm_K%geovals(ivar)%nval
   if (geovals%geovals(jvar)%nval /= self%crtm_K%geovals(ivar)%nval) call abor1_ftn("radiance_ad: error nval")
   jobs = 0
   do jprofile = 1, self%n_Profiles
   do jchannel = 1, self%n_Channels
     jobs = jobs + 1
     if (hofx%values(jobs)>rmissing) then
       do jlev = 1, geovals%geovals(jvar)%nval
         geovals%geovals(jvar)%vals(jlev,jprofile) = geovals%geovals(jvar)%vals(jlev,jprofile) &
           & + self%crtm_K%geovals(ivar)%vals(jlev,jobs) * hofx%values(jobs)
       enddo
     endif
   enddo
   enddo
   write(*,*)'radiance_ad geovals ',geovals%variables%fldnames(jvar),' min max = ', &
           & minval(geovals%geovals(jvar)%vals(:,:)),maxval(geovals%geovals(jvar)%vals(:,:))
 endif
 enddo
    
 write(*,*)'radiance_ad finished'
    
end subroutine ufo_radiance_simobs_ad
  
! ------------------------------------------------------------------------------
  
end module ufo_radiance_tlad_mod
