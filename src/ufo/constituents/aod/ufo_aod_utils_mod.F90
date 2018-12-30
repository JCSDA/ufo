MODULE ufo_aod_utils_mod

  USE iso_c_binding
  USE config_mod
  USE kinds
  USE CRTM_module
  USE ufo_vars_mod
  USE ufo_geovals_mod, ONLY: ufo_geovals, ufo_geoval, ufo_geovals_get_var
  USE ufo_basis_mod, ONLY: ufo_basis
  USE obsspace_mod
  
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: load_atm_data, load_aerosol_data
  PUBLIC :: max_string, aerosol_concentration_minvalue, aod_conf

  TYPE aod_conf
     INTEGER :: n_sensors
     INTEGER :: n_absorbers
     INTEGER :: n_clouds
     INTEGER :: n_aerosols
     INTEGER, ALLOCATABLE :: skiplist(:)
     CHARACTER(len=255), ALLOCATABLE :: sensor_id(:)
     CHARACTER(len=255) :: endian_type
     CHARACTER(len=255) :: coefficient_path
  END TYPE aod_conf

  REAL, PARAMETER :: aerosol_concentration_minvalue=1.e-16

  INTEGER, PARAMETER :: max_string=800
  LOGICAL, PARAMETER :: ice4qsat=.TRUE.


!later from aod.yaml file
       INTEGER, PARAMETER :: n_absorbers=2,n_clouds=0,n_aerosols=14

  REAL(kind_real), PARAMETER:: &
       &ttp = 2.7316e+2_kind_real, &
       &psat = 6.1078e+2_kind_real,&
       &rd = 2.8705e+2_kind_real,&
       &rv = 4.6150e+2_kind_real,&
       &cv = 7.1760e+2_kind_real,&
       &cliq = 4.1855e+3_kind_real,&
       &csol = 2.1060e+3_kind_real,&
       &cvap = 1.8460e+3_kind_real,&
       &hvap = 2.5000e+6_kind_real,&
       &hfus = 3.3358e+5_kind_real,&
       &grav = 9.81_kind_real

  REAL(kind_real), PARAMETER ::  &
       &tmix = ttp-20_kind_real,&
       &hsub = hvap+hfus,&
       &eps = rd/rv,&
       &eps_p1= one+eps,&
       &omeps=one-eps,&
       &dldt =cvap-cliq,&
       &dldti = cvap-csol,&
       &xa = -(dldt/rv),&
       &xai = -(dldti/rv),&
       &xb = xa+hvap/(rv*ttp),&
       &xbi = xai+hsub/(rv*ttp)

  CHARACTER(len=MAXVARLEN) :: varname


CONTAINS


  SUBROUTINE aod_conf_setup(rc, c_conf)

    IMPLICIT NONE
    TYPE(aod_conf), INTENT(inout) :: rc
    TYPE(c_ptr),    INTENT(in)    :: c_conf
    
    CHARACTER(len=1023) :: skipchannels
    INTEGER :: nskip, i
    CHARACTER(len=100), ALLOCATABLE :: skiplist_str(:)
    
!some config needs to come from user
!-----------------------------------
    
!number of sensors, each call to CRTM will be for a single sensor
!type (zenith/scan angle will be different)
    rc%n_sensors = 1
    
!number of absorbers, clouds and aerosols (should match what model will provide)

!@mzp begin
!    rc%n_absorbers = config_get_int(c_conf,"n_absorbers")
!    rc%n_clouds    = config_get_int(c_conf,"n_clouds"   )
!    rc%n_aerosols  = config_get_int(c_conf,"n_aerosols" )
!@mzp end    

!allocate sensor_id
    ALLOCATE(rc%sensor_id(rc%n_sensors))
    
!get sensor id from config

!@mzp begin
!    rc%sensor_id(rc%n_sensors) = config_get_string(c_conf,LEN(rc%sensor_id(rc%n_sensors)),"sensor_id")
!@mzp end
    
!endian type
!@mzp begin
!    rc%endian_type = config_get_string(c_conf,LEN(rc%endian_type),"endiantype")
!@mzp end 
   
!path to coefficient files
!@mzp begin
!    rc%coefficient_path = config_get_string(c_conf,LEN(rc%coefficient_path),"coefficientpath")
!@mzp end 
   
!channels to skip
    IF (config_element_exists(c_conf,"skipchannels")) THEN
       skipchannels = config_get_string(c_conf,LEN(skipchannels),"skipchannels")
       nskip = 1 + COUNT(TRANSFER(skipchannels, 'a', LEN(skipchannels)) == ",")
       ALLOCATE(skiplist_str(nskip))
       READ(skipchannels,*) skiplist_str
    ELSE
       nskip = 0
    ENDIF
    ALLOCATE(rc%skiplist(nskip))
    DO i = 1,nskip
       READ(skiplist_str(i),*)  rc%skiplist(i)
    ENDDO
    
  END SUBROUTINE aod_conf_setup
  

  SUBROUTINE aod_conf_delete(rc)
    
    IMPLICIT NONE
    TYPE(aod_conf), INTENT(inout) :: rc
    
    DEALLOCATE(rc%sensor_id)
    DEALLOCATE(rc%skiplist)
    
  END SUBROUTINE aod_conf_delete


  SUBROUTINE get_var_name(varname_tmplate,n,varname)
    
    CHARACTER(len=*), INTENT(in) :: varname_tmplate
    INTEGER, INTENT(in) :: n
    CHARACTER(len=*), INTENT(out) :: varname
    
    CHARACTER(len=3) :: chan
    
! pass in varname_tmplate = "brigtness_temperature"
    WRITE(chan, '(i0)') n
    varname = TRIM(varname_tmplate) // '_' // TRIM(chan) // '_'
    
  END SUBROUTINE get_var_name
  

  SUBROUTINE load_atm_data(n_profiles,n_layers,geovals,atm)

    IMPLICIT NONE

    INTEGER, INTENT(in) :: n_profiles,n_layers
    TYPE(ufo_geovals), INTENT(in) :: geovals
    TYPE(CRTM_atmosphere_type), INTENT(inout) :: atm(:)

    TYPE(ufo_geoval), POINTER :: geoval
    INTEGER :: nc, nl
    INTEGER :: k1, k2

!!$ 1   Temperature
!!$ 2   Water vapor
!!$ 3   Pressure
!!$ 4   Level pressure

!** populate the atmosphere structures for CRTM (atm(k1), for the k1-th profile)

    DO k1 = 1,n_profiles

       varname=var_t
       CALL ufo_geovals_get_var(geovals, varname,geoval)
       atm(k1)%temperature(1:n_layers) = geoval%vals(:,k1) 

!         print *, 'Temperature:', atm(k1)%Temperature(1:2), geoval%vals(1:2,k1)

       varname=var_prs
       CALL ufo_geovals_get_var(geovals, varname, geoval)

       atm(k1)%pressure(1:n_layers) = geoval%vals(:,k1) 

!         print *, 'Pressure:', atm(k1)%Pressure(1:2), geoval%vals(1:2,k1)

       varname=var_prsi
       CALL ufo_geovals_get_var(geovals, varname, geoval)
       atm(k1)%Level_Pressure(0:n_layers) = geoval%vals(:,k1)

!         print *, 'level_pressure:', atm(k1)%Level_Pressure(0:1), geoval%vals(1:2,k1)

       atm(k1)%climatology         = US_STANDARD_ATMOSPHERE

       atm(k1)%absorber_id(1:1)    = (/ H2O_ID /)
       atm(k1)%absorber_units(1:1) = (/ MASS_MIXING_RATIO_UNITS /)
       varname=var_mixr
       CALL ufo_geovals_get_var(geovals, varname, geoval)

       atm(k1)%absorber(1:n_layers,1)       = geoval%vals(:,k1) 

!         print *, 'water vapor:', atm(k1)%absorber(1:2,1), geoval%vals(1:2,k1)

       atm(k1)%absorber_id(2:2)    = (/ O3_ID /)
       atm(k1)%absorber_units(2:2) = (/ VOLUME_MIXING_RATIO_UNITS /)
       atm(k1)%absorber(1:n_layers,2)=1.e-10

    ENDDO

    NULLIFY(geoval)

  END SUBROUTINE load_atm_data


  SUBROUTINE load_aerosol_data(n_profiles,n_layers,geovals,atm)

    USE CRTM_aerosolcoeff, ONLY: aeroc

    INTEGER, INTENT(in) :: n_profiles,n_layers
    TYPE(ufo_geovals), INTENT(in) :: geovals
    TYPE(CRTM_atmosphere_type), INTENT(inout) :: atm(:)

    REAL(kind_real), DIMENSION(n_layers) :: ugkg_kgm2,qsat,rh,prsl,tsen

    TYPE(ufo_geoval), POINTER :: geoval
    INTEGER :: nc, nl
    INTEGER :: k1, k2

    INTEGER :: i,k,m

    DO m=1,n_profiles

!ug2kg && hPa2Pa
       DO k=1,n_layers
!correct for mixing ratio factor ugkg_kgm2 
!being calculated from dry pressure, cotton eq. (2.4)
!p_dry=p_total/(1+1.61*mixing_ratio)
          ugkg_kgm2(k)=1.0e-9_kind_real*(atm(m)%level_pressure(k)-&
               &atm(m)%level_pressure(k-1))*100_kind_real/grav/&
               &(one+eps_p1*atm(m)%absorber(k,1)*1.e-3_kind_real)
          prsl(k)=atm(m)%pressure(n_layers-k+1)*0.1_kind_real ! must be in cb for genqsat
          tsen(k)=atm(m)%temperature(n_layers-k+1)
       ENDDO

       CALL genqsat(qsat,tsen,prsl,n_layers,ice4qsat)

!relative humidity is ratio of specific humidities not mixing ratios
       DO k=1,n_layers
          rh(k)=(atm(m)%absorber(k,1)/(one+atm(m)%absorber(k,1)))*1.e-3_kind_real/&
               &qsat(n_layers-k+1)
       ENDDO

       DO i=1,n_aerosols_gocart_default
          varname=var_aerosols_gocart_default(i)
          CALL ufo_geovals_get_var(geovals,varname, geoval)

          atm(m)%aerosol(i)%Concentration(1:n_layers)=&
               &MAX(geoval%vals(:,m)*ugkg_kgm2,aerosol_concentration_minvalue)

          SELECT CASE ( TRIM(varname))
          CASE ('sulf')
             atm(m)%aerosol(i)%type  = SULFATE_AEROSOL
!rh needs to be from top to bottom
             DO k=1,n_layers
                atm(m)%aerosol(i)%effective_radius(k)=&
                     &gocart_aerosol_size(atm(m)%aerosol(i)%type, &
                     &rh(k))
             ENDDO

          CASE ('bc1')
             atm(m)%aerosol(i)%type  = BLACK_CARBON_AEROSOL
             atm(m)%aerosol(i)%effective_radius(:)=&
                  &AeroC%Reff(1,atm(m)%aerosol(i)%type)
          CASE ('bc2')
             atm(m)%aerosol(i)%type  = BLACK_CARBON_AEROSOL
             DO k=1,n_layers
                atm(m)%aerosol(i)%effective_radius(k)=&
                     &gocart_aerosol_size(atm(m)%aerosol(i)%type, &
                     &rh(k))
             ENDDO

          CASE ('oc1')
             atm(m)%aerosol(i)%type  = ORGANIC_CARBON_AEROSOL
             atm(m)%aerosol(i)%effective_radius(:)=&
                  &AeroC%Reff(1,atm(m)%aerosol(i)%type)
          CASE ('oc2')
             atm(m)%aerosol(i)%type  = ORGANIC_CARBON_AEROSOL
             DO k=1,n_layers
                atm(m)%aerosol(i)%effective_radius(k)=&
                     &gocart_aerosol_size(atm(m)%aerosol(i)%type, &
                     &rh(k))
             ENDDO

          CASE ('dust1')
             atm(m)%aerosol(i)%type  = DUST_AEROSOL
             atm(m)%aerosol(i)%effective_radius(:)=0.55_kind_real
          CASE ('dust2')
             atm(m)%aerosol(i)%type  = DUST_AEROSOL
             atm(m)%aerosol(i)%effective_radius(:)=1.4_kind_real
          CASE ('dust3')
             atm(m)%aerosol(i)%type  = DUST_AEROSOL
             atm(m)%aerosol(i)%effective_radius(:)=2.4_kind_real
          CASE ('dust4')
             atm(m)%aerosol(i)%type  = DUST_AEROSOL
             atm(m)%aerosol(i)%effective_radius(:)=4.5_kind_real
          CASE ('dust5')
             atm(m)%aerosol(i)%type  = DUST_AEROSOL
             atm(m)%aerosol(i)%effective_radius(:)=8.0_kind_real

          CASE ('seas1')
             atm(m)%aerosol(i)%type  = SEASALT_SSAM_AEROSOL
             DO k=1,n_layers
                atm(m)%aerosol(i)%effective_radius(k)=&
                     &gocart_aerosol_size(atm(m)%aerosol(i)%type, &
                     &rh(k))
             ENDDO
          CASE ('seas2')
             atm(m)%aerosol(i)%type  = SEASALT_SSCM1_AEROSOL
             DO k=1,n_layers
                atm(m)%aerosol(i)%effective_radius(k)=&
                     &gocart_aerosol_size(atm(m)%aerosol(i)%type, &
                     &rh(k))
             ENDDO
          CASE ('seas3')
             atm(m)%aerosol(i)%type  = SEASALT_SSCM2_AEROSOL
             DO k=1,n_layers
                atm(m)%aerosol(i)%effective_radius(k)=&
                     &gocart_aerosol_size(atm(m)%aerosol(i)%type, &
                     &rh(k))
             ENDDO
          CASE ('seas4')
             atm(m)%aerosol(i)%type  = SEASALT_SSCM3_AEROSOL
             DO k=1,n_layers
                atm(m)%aerosol(i)%effective_radius(k)=&
                     &gocart_aerosol_size(atm(m)%aerosol(i)%type, &
                     &rh(k))
             ENDDO
          END SELECT

       ENDDO

    ENDDO

    NULLIFY(geoval)

  END SUBROUTINE load_aerosol_data


  FUNCTION gocart_aerosol_size( itype, eh ) & ! eh input in 0-1
       &RESULT(r_eff)   ! in micrometer

    USE CRTM_aerosolcoeff, ONLY: aeroc
    IMPLICIT NONE

!
!   modified from a function provided by quanhua liu
!
    INTEGER ,INTENT(in) :: itype
    REAL(kind_real)    ,INTENT(in) :: eh

    INTEGER :: j1,j2,k
    REAL(kind_real)    :: h1
    REAL(kind_real)    :: r_eff

    j2 = 0
    IF ( eh <= aeroc%rh(1) ) THEN
       j1 = 1
    ELSE IF ( eh >= aeroc%rh(aeroc%n_rh) ) THEN
       j1 = aeroc%n_rh
    ELSE
       DO k = 1, aeroc%n_rh-1
          IF ( eh < aeroc%rh(k+1) .AND. eh > aeroc%rh(k) ) THEN
             j1 = k
             j2 = k+1
             h1 = (eh-aeroc%rh(k))/(aeroc%rh(k+1)-aeroc%rh(k))
             EXIT
          ENDIF
       ENDDO
    ENDIF

    IF ( j2 == 0 ) THEN
       r_eff = aeroc%reff(j1,itype )
    ELSE
       r_eff = (one-h1)*aeroc%reff(j1,itype ) + h1*aeroc%reff(j2,itype )
    ENDIF

    RETURN

  END FUNCTION gocart_aerosol_size


  SUBROUTINE genqsat(qsat,tsen,prsl,nsig,ice)

!   input argument list:
!     tsen      - input sensibile temperature field (nlat,nlon,nsig)
!     prsl      - input layer mean pressure field (nlat,nlon,nsig)
!     nsig      - number of levels                              
!     ice       - logical flag:  t=include ice and ice-water effects,
!                 depending on t, in qsat calcuations.
!                 otherwise, compute qsat with respect to water surface
!
!   output argument list:
!     qsat      - saturation specific humidity (output)
!
! remarks: see modules used
!
! attributes:
!   language: f90
!   machine:  ibm rs/6000 sp
!

    IMPLICIT NONE

    LOGICAL                               ,INTENT(in   ) :: ice
    REAL(kind_real),DIMENSION(nsig), INTENT(  out) :: qsat
    REAL(kind_real),DIMENSION(nsig),INTENT(in   ) :: tsen,prsl
    INTEGER                       ,INTENT(in   ) :: nsig

    INTEGER k
    REAL(kind_real) pw,tdry,tr,es,es2
    REAL(kind_real) w,onep3,esmax
    REAL(kind_real) desidt,deswdt,dwdt,desdt,esi,esw
    REAL(kind_real) :: mint,estmax
    INTEGER :: lmint

    onep3 = 1.e3_kind_real

    mint=340_kind_real
    lmint=1

    DO k=1,nsig
       IF((prsl(k) < 30_kind_real .AND.  &
            prsl(k) > 2_kind_real) .AND.  &
            tsen(k) < mint)THEN
          lmint=k
          mint=tsen(k)
       END IF
    END DO

    tdry = mint
    tr = ttp/tdry

    IF (tdry >= ttp .OR. .NOT. ice) THEN
       estmax = psat * (tr**xa) * EXP(xb*(one-tr))
    ELSEIF (tdry < tmix) THEN
       estmax = psat * (tr**xai) * EXP(xbi*(one-tr))
    ELSE
       w  = (tdry - tmix) / (ttp - tmix)
       estmax =  w * psat * (tr**xa) * EXP(xb*(one-tr)) &
            + (one-w) * psat * (tr**xai) * EXP(xbi*(one-tr))
    ENDIF

    DO k = 1,nsig

       tdry = tsen(k)
       tr = ttp/tdry
       IF (tdry >= ttp .OR. .NOT. ice) THEN
          es = psat * (tr**xa) * EXP(xb*(one-tr))
       ELSEIF (tdry < tmix) THEN
          es = psat * (tr**xai) * EXP(xbi*(one-tr))
       ELSE
          esw = psat * (tr**xa) * EXP(xb*(one-tr)) 
          esi = psat * (tr**xai) * EXP(xbi*(one-tr)) 
          w  = (tdry - tmix) / (ttp - tmix)
          es =  w * psat * (tr**xa) * EXP(xb*(one-tr)) &
               + (one-w) * psat * (tr**xai) * EXP(xbi*(one-tr))
       ENDIF

       pw = onep3*prsl(k)
       esmax = es
       IF(lmint < k)THEN
          esmax=0.1_kind_real*pw
          esmax=MIN(esmax,estmax)
       END IF
       es2=MIN(es,esmax)
       qsat(k) = eps * es2 / (pw - omeps * es2)

    END DO

    RETURN

  END SUBROUTINE genqsat


END MODULE ufo_aod_utils_mod

