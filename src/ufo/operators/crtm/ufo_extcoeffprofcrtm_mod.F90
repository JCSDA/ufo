! (C) Copyright 2025-2026 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module to handle extinction coefficient profile observations

module ufo_extcoeffprofcrtm_mod

 use fckit_configuration_module, only: fckit_configuration
 use fckit_mpi_module,   only: fckit_mpi_comm
 use , intrinsic :: iso_c_binding
 use kinds
 use missing_values_mod

 use ufo_geovals_mod, only: ufo_geovals, ufo_geoval, ufo_geovals_get_var
 use ufo_vars_mod
 use ufo_crtm_utils_mod
 use crtm_module
 use ODPS_CoordinateMapping, only: Geopotential_Height
 use obsspace_mod

 use vert_interp_mod, only: vert_interp_weights
 use ufo_constants_mod, only : thousand

 implicit none
 private

 !> Fortran derived type for aod trajectory
 type, public :: ufo_extcoeffprofcrtm
 private
  character(len=MAXVARLEN), public, allocatable :: varin(:), varin_aero(:) ! variablesrequested from the model
  integer, allocatable                          :: channels(:)
  type(crtm_conf) :: conf
 contains
   procedure :: setup  => ufo_extcoeffprofcrtm_setup
   procedure :: delete => ufo_extcoeffprofcrtm_delete
   procedure :: simobs => ufo_extcoeffprofcrtm_simobs
 end type ufo_extcoeffprofcrtm
 CHARACTER(len=maxvarlen), DIMENSION(6), PARAMETER :: varin_default = [var_ts, var_mixr, var_rh, var_prs, var_prsi, var_sfc_z]

contains

! Remaps a per-model-layer quantity (profile_model, one value per
! CRTM model layer) onto the obs layers' own height grid (profile_obslayers)
subroutine interpolate_to_obslayers(nlayers_obs, nlayers_model, iheight_obs, &
                                    iheight_model, profile_model, profile_obslayers)
  implicit none
  integer, intent(in   ) :: nlayers_obs, nlayers_model
  real(kind_real), intent(in   ), dimension(nlayers_obs+1)   :: iheight_obs
  real(kind_real), intent(in   ), dimension(nlayers_model+1) :: iheight_model
  real(kind_real), intent(in   ), dimension(nlayers_model)   :: profile_model
  real(kind_real) :: wf_a, wf_b
  real(kind_real), intent(inout), dimension(nlayers_obs)     :: profile_obslayers
  integer :: j, k, wi_a, wi_b

  do k=1, nlayers_obs

     call vert_interp_weights(nlayers_model+1, iheight_obs(k),   iheight_model, wi_a, wf_a)
     call vert_interp_weights(nlayers_model+1, iheight_obs(k+1), iheight_model, wi_b, wf_b)

     profile_obslayers(k) = zero
     ! when multiple geovals levels are in a obs layer
     if ( wi_a < wi_b ) then
        profile_obslayers(k) = profile_obslayers(k) + profile_model(wi_a) * wf_a
        do j = wi_a+1, wi_b-1
           profile_obslayers(k) = profile_obslayers(k) + profile_model(j)
        end do
        profile_obslayers(k) = profile_obslayers(k) + profile_model(wi_b) * (one-wf_b)

     ! when multiple obs layers are in a geovals level
     else if ( wi_a == wi_b ) then
        profile_obslayers(k) = profile_obslayers(k) + profile_model(wi_a) * (wf_a-wf_b)
     end if
  end do
end subroutine interpolate_to_obslayers

! ------------------------------------------------------------------------------

subroutine ufo_extcoeffprofcrtm_setup(self, f_confOper, channels, midPointJulday, comm)

implicit none
class(ufo_extcoeffprofcrtm), intent(inout) :: self
type(fckit_configuration), intent(in)      :: f_confOper

! List of CRTM/wavelength channels.
integer(c_int),            intent(in)    :: channels(:)
integer(c_int64_t),        intent(in)    :: midPointJulday
type(fckit_mpi_comm),      intent(in)    :: comm

integer                   :: nvars_in
character(len=max_string) :: err_msg
type(fckit_configuration) :: f_confOpts

CHARACTER(len=MAXVARLEN), ALLOCATABLE :: var_aerosols(:)

 call f_confOper%get_or_die("obs options",f_confOpts)

 call crtm_conf_setup(self%conf, f_confOpts, f_confOper, midPointJulday, comm)
 if ( ufo_vars_getindex(self%conf%Absorbers, var_mixr) /= 1 ) then
   write(err_msg,*) "ufo_extcoeffprofcrtm_setup error: H2O must be first in CRTM Absorbers for CRTM_AOD"
   call abor1_ftn(err_msg)
 end if
 if ( ufo_vars_getindex(self%conf%Absorbers, var_oz) < 2 ) then
   write(err_msg,*) "ufo_extcoeffprofcrtm_setup error: O3 must be included in CRTM Absorbers"
   call abor1_ftn(err_msg)
 end if

 CALL assign_aerosol_names(self%conf%aerosol_option, var_aerosols)

 nvars_in = SIZE(varin_default)+SIZE(var_aerosols)
 allocate(self%varin(nvars_in))
 self%varin(1:size(varin_default)) = varin_default
 self%varin(SIZE(varin_default)+1:) = var_aerosols

 allocate(self%varin_aero(SIZE(var_aerosols)))
 self%varin_aero(:) = var_aerosols(:)

 ! save channels
 allocate(self%channels(size(channels)))
 self%channels(:) = channels(:)

end subroutine ufo_extcoeffprofcrtm_setup

! ------------------------------------------------------------------------------

subroutine ufo_extcoeffprofcrtm_delete(self)

implicit none
class(ufo_extcoeffprofcrtm), intent(inout) :: self

 call crtm_conf_delete(self%conf)

end subroutine ufo_extcoeffprofcrtm_delete

! ------------------------------------------------------------------------------

SUBROUTINE ufo_extcoeffprofcrtm_simobs(self, geovals, obss, nvars, nlocs, nobslayers, hofx)
use fckit_mpi_module,   only : fckit_mpi_comm
use obsspace_mod
use string_f_c_mod

implicit none
class(ufo_extcoeffprofcrtm),     intent(in)    :: self
type(ufo_geovals),               intent(in)    :: geovals
integer,                         intent(in)    :: nvars, nlocs, nobslayers
real(c_double),                  intent(inout) :: hofx(nvars, nlocs)
type(c_ptr), value,              intent(in)    :: obss

! Local Variables
! Locations are one-per-(profile,layer), grouped into records (one record per
! physical profile) via "obsgrouping" on the obs space's "sequenceNumber" MetaData.
! nobslayers is the *maximum* layers one profile can have,
! n_profileLayerCnt is each profile's actual layer count.
! heightTop/heightBottom are per-location (each layer's own interface bounds);
! obs_iheight is gathered from these into a per-profile view,
! and obs_thick is just the difference between each layer's bounds.

character(*), parameter :: PROGRAM_NAME = "ufo_extcoeffprofcrtm_mod.F90"
character(len=MAXVARLEN) :: def_aero_mod
character(255) :: message
integer        :: err_stat, alloc_stat
integer        :: j, k, l, m, n, kobs, iloc, n_profileLayerCnt
integer        :: nrecs, icount, iobs
type(ufo_geoval), pointer :: temp, sfc_z
real(c_double) :: missing

integer :: n_Profiles
integer :: n_Layers
integer :: n_Channels
integer :: n_validlayers
integer :: h2o_idx

logical :: jacobian_needed

integer(c_size_t), allocatable                :: recnum(:)
integer, allocatable                          :: nlocs_begin(:), nlocs_end(:), row_indices(:)
real(kind_real), dimension(nlocs)             :: heightTop, heightBottom
real(kind_real), dimension(nobslayers)        :: obs_thick
real(kind_real), dimension(nobslayers+1)      :: obs_iheight
real(kind_real), allocatable, dimension(:, :) :: profile_model, profile_obslayers

! Define the Channel Info and Geometry  arguments
type(CRTM_ChannelInfo_type)             :: chinfo(self%conf%n_Sensors)
type(CRTM_Geometry_type),   allocatable :: geo(:)

! Define the FORWARD variables
type(CRTM_Atmosphere_type), allocatable :: atm(:)
type(CRTM_Surface_type),    allocatable :: sfc(:)
type(CRTM_RTSolution_type), allocatable :: rts(:,:)

! Define the K-MATRIX variables - necessary for CRTM_AOD call
! ------------------------------------------------------------
type(CRTM_Atmosphere_type), allocatable :: atm_K(:,:)
type(CRTM_RTSolution_type), allocatable :: rts_K(:,:)

type(fckit_mpi_comm)  :: f_comm

 call obsspace_get_comm(obss, f_comm)

 ! Group flattened (profile,layer) locations back into physical profiles
 ! (records), Requires "obsgrouping on [sequenceNumber]"
 nrecs = obsspace_get_nrecs(obss)

 ! for MPI parallel, if there are no records on this rank, return early to avoid
 ! allocating arrays of size 0, which is not allowed in Fortran.
 if (nrecs==0) return

 allocate(recnum(nlocs), nlocs_begin(nrecs), nlocs_end(nrecs))
 call obsspace_get_recnum(obss, recnum)

 nlocs_begin = 1
 nlocs_end = 1
 icount = 1
 do iobs = 1, nlocs-1
   if (recnum(iobs+1) /= recnum(iobs)) then
     icount = icount + 1
     nlocs_end(icount-1) = iobs
     nlocs_begin(icount) = iobs+1
   end if
 end do
 nlocs_end(nrecs) = nlocs
 if (icount /= nrecs) then
   write(message,*) "ufo_extcoeffprofcrtm_simobs error: record number is not consistent: ", icount, nrecs
   call abor1_ftn(message)
 end if
 n_Profiles = nrecs

 ! The obs operator's locations() builds one interpolation path
 ! per profile via recidx_begin()/recidx_end() [ascending record-number]
 ! And geovals column m may not be profile m's data;
 ! so row_indices maps profile m to its actual geovals column.
 allocate(row_indices(n_Profiles))
 do m = 1, n_Profiles
   row_indices(m) = count(recnum(nlocs_begin) < recnum(nlocs_begin(m))) + 1
 end do

 call ufo_geovals_get_var(geovals, var_ts, temp)
 n_Layers = temp%nval
 nullify(temp)

 call obsspace_get_db(obss, "MetaData", "heightTop", heightTop)
 call obsspace_get_db(obss, "MetaData", "heightBottom", heightBottom)

 call ufo_geovals_get_var(geovals, var_sfc_z, sfc_z)

 ! Program header
 ! --------------
 ! call CRTM_Version( Version )
 ! call Program_Message( PROGRAM_NAME, &
 !                       "Check/example program for the CRTM Forward and K-Matrix functions using "//&
 !                       trim(self%conf%ENDIAN_type)//" coefficient datafiles", &
 !                       "CRTM Version: "//TRIM(Version) )


 ! Initialise all the sensors at once
 ! ----------------------------------
 !** NOTE: CRTM_Init points to the various binary files needed for CRTM.  See the
 !**       CRTM_Lifecycle.f90 for more details.

 ! write( *,"(/5x,"Initializing the CRTM...")" )
 call define_aerosol_model(self%conf%AerosolCoeff_File, def_aero_mod)
 err_stat = CRTM_Init( self%conf%SENSOR_ID, &
            chinfo, &
            File_Path           = trim(self%conf%COEFFICIENT_PATH), &
            NC_File_Path        = trim(self%conf%NC_COEFFICIENT_PATH), &
            SpcCoeff_Format     = trim(self%conf%SpcCoeff_Format), &
            TauCoeff_Format     = trim(self%conf%TauCoeff_Format), &
            Aerosol_Model       = trim(def_aero_mod), &
            AerosolCoeff_Format = trim(self%conf%AerosolCoeff_Format), &
            AerosolCoeff_File   = trim(self%conf%AerosolCoeff_File), &
            Quiet=.TRUE.)
 message = "Error initializing CRTM"
 call crtm_comm_stat_check(err_stat, PROGRAM_NAME, message, f_comm)



 ! Loop over all sensors. Not necessary if we"re calling CRTM for each sensor
 ! ----------------------------------------------------------------------------
 Sensor_Loop:DO n = 1, self%conf%n_Sensors

   ! Determine the number of channels for the current sensor
   ! -------------------------------------------------------
   n_Channels = CRTM_ChannelInfo_n_Channels(chinfo(n))


   ! Allocate the ARRAYS
   ! -------------------
   allocate( geo( n_Profiles ),               &
             atm( n_Profiles ),               &
             sfc( n_Profiles ),               &
             rts( n_Channels, n_Profiles ),   &
             profile_model(n_Channels, n_Layers),     &
             profile_obslayers(size(self%channels), nobslayers), &
             STAT = alloc_stat )
   message = "Error allocating structure arrays"
   call crtm_comm_stat_check(alloc_stat, PROGRAM_NAME, message, f_comm)

   ! Create the input FORWARD structure (atm)
   ! ----------------------------------------
   call CRTM_Atmosphere_Create( atm, n_Layers, self%conf%n_Absorbers, self%conf%n_Clouds, self%conf%n_Aerosols )
   if ( any(.not. CRTM_Atmosphere_Associated(atm)) ) then
      message = "Error allocating CRTM Forward Atmosphere structure"
      call Display_Message( PROGRAM_NAME, message, FAILURE )
      stop
   end if

   call CRTM_RTSolution_Create(rts, n_Layers )

   ! Assign the data from the GeoVaLs
   ! ---------------------------------
   call Load_Atm_Data(n_Profiles,n_Layers,geovals,atm,self%conf,row_indices=row_indices)

   if (trim(self%conf%aerosol_option) /= "") then
     call load_aerosol_data(n_Profiles, n_layers, geovals, &
       & self%conf, self%varin_aero, trim(def_aero_mod), atm, row_indices=row_indices)
   end if

   ! Call THE CRTM inspection
   ! ------------------------
   if (self%conf%inspect > 0) then
     call CRTM_Atmosphere_Inspect(atm(self%conf%inspect))
     call CRTM_ChannelInfo_Inspect(chinfo(n))
   end if

   ! Start processing the CRTM_AOD
   ! |-> Jacobian K-Matrix .OR.
   ! |-> Forward Operator
   ! ------------------------------------------------
   jacobian_needed = .false.
   kmatrix : select case ( jacobian_needed )

     case(.true.)
       !
       ! Description:
       ! ===========
       ! For the time being this branch is NEVER executed because
       ! jacobian_needed = .false. . The Jacobian may be needed in the future
       ! for QC and bias correction but running the Forward operator alone for
       ! ufo simobs requires less resources. The CRTM_AOD Jacobian code will
       ! remain unexecuted until the necessary configuration code is added.
       !

       ! The output K-MATRIX structure:

       ! Allocate the k-matrix arrays
       allocate( atm_K( n_channels, N_PROFILES ), &
                 rts_K( n_channels, N_PROFILES ), &
                 STAT = alloc_stat )
       if ( alloc_stat /= 0 ) then
         message = "Error allocating structure arrays"
         call Display_Message( PROGRAM_NAME, message, FAILURE )
         stop
       end if

       ! Call the constructor for the k-matrix output structure
       call CRTM_Atmosphere_Create( atm_K, n_layers, self%conf%n_Absorbers, self%conf%n_Clouds, self%conf%n_Aerosols)
       if ( any(.not. CRTM_Atmosphere_Associated(atm_K)) ) then
         message = "Error allocating CRTM K-matrix Atmosphere structure"
         call Display_Message( PROGRAM_NAME, message, FAILURE )
         stop
       end if

       ! Call the constructor for the k-matrix input structure
       call CRTM_RTSolution_Create(rts_k, n_Layers )

       ! ============================================================================
       !    **** INITIALIZE THE K-MATRIX ARGUMENTS ****
       !
       !     Zero the K-matrix OUTPUT structures
       ! ---------------------------------------
       call CRTM_Atmosphere_Zero( atm_k )

       !     Initialize the K-matrix INPUT
       ! ----------------------------------
       do m = 1, n_Profiles
         do l = 1, n_Channels
           rts_K(l, m)%Layer_Optical_Depth = ONE
         end do
       end do
       ! ============================================================================


       !     The K-matrix model y = K.x for AOD
       ! -----------------------------------------
       err_stat = CRTM_AOD_K( atm,    &  ! FORWARD  Input
            rts_K                   , &  ! K-MATRIX Input
            chinfo(n:n)             , &  ! Input
            rts                     , &  ! FORWARD  Output
            atm_k        )               ! K-MATRIX Output

       message = "Error calling CRTM K-Matrix Model for "//TRIM(self%conf%SENSOR_ID(n))
       call crtm_comm_stat_check(err_stat, PROGRAM_NAME, message, f_comm)

       call CRTM_Atmosphere_Destroy(atm_k)
       call CRTM_RTSolution_Destroy(rts_k)

       deallocate( atm_K, rts_K, STAT = alloc_stat )

       message = "Error deallocating Jacobian structure arrays"
       call crtm_comm_stat_check(alloc_stat, PROGRAM_NAME, message, f_comm)

     case default

       !     The Forward Operator y = H(x) for AOD
       ! ------------------------------------------
       err_stat = CRTM_AOD( atm          , &  ! FORWARD  Input
                            chinfo(n:n)  , &  ! Input
                            rts          )    ! FORWARD  Output

       message = "Error calling CRTM Forward Model for "//TRIM(self%conf%SENSOR_ID(n))
       call crtm_comm_stat_check(err_stat, PROGRAM_NAME, message, f_comm)

   end select kmatrix

   ! Put extinction coefficient into hofx
   ! ----------------------------------------------
   ! Set missing value
   missing = missing_value(missing)
   hofx = missing

   ! Identify h2o_idx
   h2o_idx = 1
   DO j = 1,atm(1)%n_Absorbers
        IF (atm(1)%Absorber_ID(j) == H2O_ID) h2o_idx = j
   END DO

   do m = 1, n_Profiles
      profile_obslayers = missing

      ! This profile's own layer count -- its record's location range.
      n_profileLayerCnt = nlocs_end(m) - nlocs_begin(m) + 1
      if (n_profileLayerCnt > nobslayers) then
        write(message,*) "ufo_extcoeffprofcrtm_simobs error: profile has more layers "// &
                          "than maxlayer allows: ", n_profileLayerCnt, nobslayers
        call abor1_ftn(message)
      end if

      ! Gather this profile's per-layer interface heights out of per-location arrays
      ! layer k of profile m lives at location nlocs_begin(m)+k-1.
      do k = 1, n_profileLayerCnt
         iloc = nlocs_begin(m) + k - 1
         obs_iheight(k) = heightTop(iloc)
      end do
      obs_iheight(n_profileLayerCnt+1) = heightBottom(nlocs_end(m))

      ! Layer thickness is just the difference between its own top and bottom bounds.
      do k = 1, n_profileLayerCnt
         obs_thick(k) = obs_iheight(k) - obs_iheight(k+1)
      end do

      ! Calculate geopotential height at layer interface
      call Geopotential_Height(atm(m)%Level_Pressure, &
                               atm(m)%Temperature, &
                               atm(m)%Absorber(:, h2o_idx), &
                               sfc_z%vals(1,row_indices(m))/thousand, &   ! geopotential_height_at_surface in meter
                               atm(m)%Height)
      atm(m)%Height = atm(m)%Height * thousand  ! km to m

      ! Remove obs layers that fall below the model's lowest computed layer height
      n_validlayers = n_profileLayerCnt
      do kobs = 1, n_profileLayerCnt+1
         if (obs_iheight(kobs) < atm(m)%Height(n_Layers) ) then
            n_validlayers = kobs - 1
            exit
         end if
      end do
      ! in case "channels" ever configures a subset of sensor's complete channel list
      do l = 1, size(self%channels)
         ! AOD conserved interpolation from extinction profile (unitless) to observation layers
         call interpolate_to_obslayers(n_validlayers, n_Layers, obs_iheight(1:n_validlayers+1), atm(m)%Height, &
                                       rts(self%channels(l), m)%Layer_Optical_Depth, profile_obslayers(l,1:n_validlayers))

         ! Convert unitless optical depth to an extinction coefficient: dividing
         ! by thickness (m) gives m-1, then x1000 converts m-1 to km-1
         do k = 1, n_validlayers
            profile_obslayers(l, k) = profile_obslayers(l, k) / obs_thick(k) * thousand
         end do
      end do

      ! Save this profile's interpolated layers into hofx(nvars, nlocs): layer k
      ! of profile m goes to flattened location nlocs_begin(m)+k-1.
      do k = 1, n_profileLayerCnt
         iloc = nlocs_begin(m) + k - 1
         do l = 1, nvars
            hofx(l, iloc) = profile_obslayers(l, k)
         end do
      end do

   end do

   ! Deallocate the structures
   ! -------------------------
   call CRTM_Atmosphere_Destroy(atm)
   call CRTM_RTSolution_Destroy(rts)

   ! Deallocate all arrays
   ! ---------------------
   deallocate(geo, atm, sfc, rts, STAT = alloc_stat)
   message = "Error deallocating structure arrays"
   call crtm_comm_stat_check(alloc_stat, PROGRAM_NAME, message, f_comm)

end do Sensor_Loop

 ! Destroy CRTM instance
 ! ---------------------
 ! write( *, "( /5x, "Destroying the CRTM..." )" )
 err_stat = CRTM_Destroy( chinfo )
 message = "Error destroying CRTM"
 call crtm_comm_stat_check(err_stat, PROGRAM_NAME, message, f_comm)

 call f_comm%final()

 deallocate(recnum, nlocs_begin, nlocs_end, row_indices)

end subroutine ufo_extcoeffprofcrtm_simobs

! ------------------------------------------------------------------------------

end module ufo_extcoeffprofcrtm_mod
