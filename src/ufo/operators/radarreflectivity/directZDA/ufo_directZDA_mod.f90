! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module for directZDA observation operator
!
module ufo_directZDA_mod

 use kinds
 use ufo_vars_mod
 use oops_variables_mod
 use obs_variables_mod
 
 implicit none
 private

!> Fortran derived type for the observation type

 type, public :: ufo_directZDA
 private
   type(obs_variables), public :: obsvars
   type(oops_variables), public :: geovars 

   character(len=MAXVARLEN), public :: v_coord ! GeoVaL to use to interpolate in vertical 
   character(len=MAXVARLEN), public :: micro_option    ! Choice (enum) of microphysics option
   integer, public :: mphyopt                  ! Integer value for microphysics option
   logical, public :: use_variational = .false.

 contains 
   procedure :: setup  => ufo_directZDA_setup
   procedure :: simobs => ufo_directZDA_simobs
 end type ufo_directZDA
 
 integer :: n_geovars
 character(len=maxvarlen), dimension(:), allocatable :: geovars_list

contains

! ------------------------------------------------------------------------------
subroutine ufo_directZDA_setup(self, yaml_conf)
  use fckit_configuration_module, only: fckit_configuration
  use iso_c_binding
  use radarz_iface

  implicit none
  class(ufo_directZDA), intent(inout)     :: self
  type(fckit_configuration), intent(in) :: yaml_conf

  character(kind=c_char,len=:), allocatable :: coord_name
  character(kind=c_char,len=:), allocatable :: micro_option
  character(kind=c_char,len=:), allocatable :: this_varname
  character(len=maxvarlen) :: var_string
  integer :: iret_init_mphyopt = -1

  ! YAML option for microphysics scheme
  call yaml_conf%get_or_die("microphysics option", micro_option)
  if(trim(micro_option) .eq. "Thompson") then
    !self%micro_option = THOMPSON
    self%mphyopt = 108
    mp_option = 108
  else if(trim(micro_option) .eq. "NSSL") then
    !self%micro_option = NSSL
    self%mphyopt = 14
    mp_option = 14
  else if(trim(micro_option) .eq. "Lin") then
    !self%micro_option = LIN
    self%mphyopt = 2
    mp_option = 2
  else if(trim(micro_option) .eq. "GFDL") then
    !self%micro_option = GFDL
    self%mphyopt = 5
    mp_option = 5
  else
    print*, ' microphysics picked is: ', trim(micro_option)
    call abor1_ftn("microphysics option not set or unsupported, aborting")
  endif

  ! If YAML option indicates use of variational method then update mp_option
  self%use_variational=.false.
  call yaml_conf%get_or_die("use variational method", self%use_variational)

  if (self%use_variational) then
    if (trim(micro_option) .ne. "Thompson" .and. trim(micro_option) .ne. "Lin" ) then ! Jun: Or we can use if ( NSSL) 
      call abor1_ftn("variational method not available for requested microphysics option")
    endif
  endif

  ! YAML option for vertical coordinate name
  self%v_coord = var_z
  if( yaml_conf%has("vertical coordinate") ) then
      call yaml_conf%get_or_die("vertical coordinate",coord_name)
      self%v_coord = coord_name
      if( trim(self%v_coord) .ne. var_z ) then
        call abor1_ftn("ufo_directZDA: incorrect vertical coordinate specified")
      endif
  endif

  ! Set up the atmospheric state variables and microphysics species (into geovars_list).

  if ( self%mphyopt .eq. 14 ) then         !  NSSL EnKF OP: 14
    n_geovars=13 
    if (.not.allocated(geovars_list) ) allocate(geovars_list(n_geovars))
    var_string="var_rain_mixing_ratio"
    if( yaml_conf%has(trim(var_string)) ) then
       !call yaml_conf%get_or_die(trim(var_string), geovars_list(1))
       call yaml_conf%get_or_die("var_rain_mixing_ratio", this_varname)
       geovars_list(1) = this_varname
    endif
    var_string="var_snow_mixing_ratio"
    if( yaml_conf%has(trim(var_string)) ) then
       call yaml_conf%get_or_die(trim(var_string), this_varname)
       geovars_list(2) = this_varname
    endif
    var_string="var_graupel_mixing_ratio"
    if( yaml_conf%has(trim(var_string)) ) then
       call yaml_conf%get_or_die(trim(var_string), this_varname)
       geovars_list(3) = this_varname
    endif
    var_string="var_hail_mixing_ratio"
    if( yaml_conf%has(trim(var_string)) ) then
       call yaml_conf%get_or_die(trim(var_string), this_varname)
       geovars_list(4) = this_varname
    endif
    var_string="var_rain_number_concentration"
    if( yaml_conf%has(trim(var_string)) ) then
       call yaml_conf%get_or_die(trim(var_string), this_varname)
       geovars_list(5) = this_varname
    endif
    var_string="var_snow_number_concentration"
    if( yaml_conf%has(trim(var_string)) ) then
       call yaml_conf%get_or_die(trim(var_string), this_varname)
       geovars_list(6) = this_varname
    endif
    var_string="var_graupel_number_concentration"
    if( yaml_conf%has(trim(var_string)) ) then
       call yaml_conf%get_or_die(trim(var_string), this_varname)
       geovars_list(7) = this_varname
    endif
    var_string="var_hail_number_concentration"
    if( yaml_conf%has(trim(var_string)) ) then
       call yaml_conf%get_or_die(trim(var_string), this_varname)
       geovars_list(8) = this_varname
    endif
    var_string="var_graupel_vol_mixing_ratio"
    if( yaml_conf%has(trim(var_string)) ) then
       call yaml_conf%get_or_die(trim(var_string), this_varname)
       geovars_list(9) = this_varname
    endif
    var_string="var_hail_vol_mixing_ratio"
    if( yaml_conf%has(trim(var_string)) ) then
       call yaml_conf%get_or_die(trim(var_string), this_varname)
       geovars_list(10) = this_varname
    endif
    geovars_list(11) = var_prs
    geovars_list(12) = var_ts
    geovars_list(13) = var_q
  else if ( self%mphyopt .eq. 108 ) then   ! TM OP : 108
    n_geovars=7
    if (.not.allocated(geovars_list) ) allocate(geovars_list(n_geovars))
    var_string="var_rain_mixing_ratio"
    if( yaml_conf%has(trim(var_string)) ) then
       call yaml_conf%get_or_die(trim(var_string), this_varname)
       geovars_list(1) = this_varname
    endif
    var_string="var_snow_mixing_ratio"
    if( yaml_conf%has(trim(var_string)) ) then
       call yaml_conf%get_or_die(trim(var_string), this_varname)
       geovars_list(2) = this_varname
    endif
    var_string="var_graupel_mixing_ratio"
    if( yaml_conf%has(trim(var_string)) ) then
       call yaml_conf%get_or_die(trim(var_string), this_varname)
       geovars_list(3) = this_varname
    endif
    var_string="var_rain_number_concentration"
    if( yaml_conf%has(trim(var_string)) ) then
       call yaml_conf%get_or_die(trim(var_string), this_varname)
       geovars_list(4) = this_varname
    endif
    geovars_list(5) = var_prs
    geovars_list(6) = var_ts
    geovars_list(7) = var_q
  else if ( self%mphyopt .eq. 2 .or. self%mphyopt .eq. 5 ) then   ! LIN=2, GFDL=5
    n_geovars=6
    if (.not.allocated(geovars_list) ) allocate(geovars_list(n_geovars))
    var_string="var_rain_mixing_ratio"
    if( yaml_conf%has(trim(var_string)) ) then
       call yaml_conf%get_or_die(trim(var_string), this_varname)
       geovars_list(1) = this_varname
    endif
    var_string="var_snow_mixing_ratio"
    if( yaml_conf%has(trim(var_string)) ) then
       call yaml_conf%get_or_die(trim(var_string), this_varname)
       geovars_list(2) = this_varname
    endif
    var_string="var_graupel_mixing_ratio"
    if( yaml_conf%has(trim(var_string)) ) then
       call yaml_conf%get_or_die(trim(var_string), this_varname)
       geovars_list(3) = this_varname
    endif
    geovars_list(4) = var_prs
    geovars_list(5) = var_ts
    geovars_list(6) = var_q
  end if

  call self%geovars%push_back(geovars_list)
  call self%geovars%push_back(self%v_coord)

  iret_init_mphyopt = init_mphyopt(self%mphyopt)
  if ( iret_init_mphyopt .ne. 0 ) then
     call abor1_ftn("microphysics option not set properly, aborting")
  end if


end subroutine ufo_directZDA_setup

! ------------------------------------------------------------------------------
! Code in this routine is for radarreflectivity only
subroutine ufo_directZDA_simobs(self, geovals, obss, nvars, nlocs, hofx)
  use vert_interp_mod
  use ufo_geovals_mod, only: ufo_geovals, ufo_geoval, ufo_geovals_get_var
  use obsspace_mod

  use radarz_iface, only: P_qr, P_qs, P_qg, P_qh,   &
                          P_nr, P_ns, P_ng, P_nh,   & 
                          P_vg, P_vh, qgh_opt, ta,  &
                          t_obs_dual, t_para_dsd, nscalar
  use radarz_module, only: init_refl, init_para_dsd, calcMDR, calcMu, &
                           set_dsd_para, rdr_obs
  use ufo_directZDA_util_mod, only : coef4dbzfwrd, Cr, Pr, &
                                     Cs_dry, Cs_wet, &
                                     Cg_dry, Cg_wet, &
                                     Ps_dry, Ps_wet, Pg_dry, Pg_wet, &
                                     calc_coeffs_dry_snow_tm

  implicit none
  class(ufo_directZDA), intent(in)  :: self
  integer, intent(in)               :: nvars, nlocs
  type(ufo_geovals),  intent(in)    :: geovals
  real(c_double),     intent(inout) :: hofx(nvars, nlocs)
  type(c_ptr), value, intent(in)    :: obss

  ! Local variables
  integer :: iobs, ivar, nvars_geovars
  real(kind_real),  dimension(:), allocatable :: obsvcoord
  type(ufo_geoval), pointer :: vcoordprofile, profile
  real(kind_real),  allocatable :: wf(:)
  integer,          allocatable :: wi(:)

  character(len=MAXVARLEN) :: geovar

  real(kind_real), allocatable :: tmp(:) 
  real(kind_real) :: tmp2 
  real(kind_real), allocatable :: vfields(:,:)  ! background fields interplated vertically to the observation height

  integer :: i_melt_snow, i_melt_graupel
  real(kind_real) :: qrges, qsges, qgges, qhges
  real(kind_real) :: qnrges, qnsges, qngges, qnhges
  real(kind_real) :: qvgges, qvhges
  real(kind_real) :: rdBZ, rdBZr, rdBZs, rdBZg, rdBZh
  real(kind_real) :: P1D, Q1D, T1D, RHO
  real(kind_real) :: Ze, Zer, Zes, Zeg, Zeh
  real(kind_real) :: wgt_dry, wgt_wet
  real(kind_real) :: Zeg_dry, Zeg_wet

  real(kind_real), parameter :: qx_min = 1.0E-8_kind_real
  real(kind_real), parameter :: qn_min = 1.0_kind_real
  real(kind_real), parameter :: qb_min = 1.0E-8_kind_real

  real(kind_real), parameter :: rd=287.04_kind_real
  real(kind_real), parameter :: one=1._kind_real
  real(kind_real), parameter :: D608=0.608_kind_real
  real(kind_real), parameter :: zero=0.0_kind_real
  real(kind_real), parameter :: ten=10.0_kind_real
  real(kind_real), parameter :: T_melt=273.15_kind_real

  real(kind_real), parameter :: pi = 3.141592_kind_real  ! pi
  real(kind_real), parameter :: rhor=1000._kind_real     ! Density of rain (kg m**-3)
  real(kind_real), parameter :: rhoh=913._kind_real      ! Density of hail (kg m**-3)
  real(kind_real), parameter :: rhos=100._kind_real      ! Density of snow (kg m**-3)
  real(kind_real), parameter :: rhog=400._kind_real      ! Density of graupel (kg m**-3)
  real(kind_real),    parameter :: am_s = 0.069_kind_real ! from WRF
  real(kind_real) :: a_dry_snow_tm, b_dry_snow_tm ! dry snow coeffs for TM
  real(kind_real) :: smoz, alpha_const_tm_dry_snow
  real(kind_real) :: oams

  integer :: idx_p, idx_t, idx_q
  integer :: iret_coef4dbzfwrd

! variables added for radarZ
  real(kind_real),dimension(nscalar) :: qscalar
  type(t_obs_dual) :: obs_dual
  type(t_para_dsd) :: var_dsd

  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
  ! Get height profiles from geovals
  call ufo_geovals_get_var(geovals, self%v_coord, vcoordprofile)

  ! Get the observation vertical coordinates
  allocate(obsvcoord(nlocs))
  call obsspace_get_db(obss, "MetaData", "height", obsvcoord)

  ! Allocate arrays for interpolation weights

  allocate(wi(nlocs))
  allocate(wf(nlocs))

  ! Calculate the interpolation weights
  allocate(tmp(vcoordprofile%nval))
  do iobs = 1, nlocs
    tmp = vcoordprofile%vals(:,iobs)
    tmp2 = obsvcoord(iobs)
    call vert_interp_weights(vcoordprofile%nval, tmp2, tmp, wi(iobs), wf(iobs))
  enddo

! Number of variables in geovars (without the vertical coordinate)
  nvars_geovars = self%geovars%nvars() - 1
  allocate(vfields(nvars_geovars,nlocs))

  do ivar = 1, nvars_geovars
    ! Get the name of input variable in geovals
    geovar = self%geovars%variable(ivar)

    ! Get profile for this variable from geovals
    call ufo_geovals_get_var(geovals, geovar, profile)

    ! Interpolate from geovals to observational location into hofx
    do iobs = 1, nlocs
      call vert_interp_apply(profile%nval, profile%vals(:,iobs), &
                             & vfields(ivar,iobs), wi(iobs), wf(iobs))
    enddo
  enddo

  if ( self%mphyopt .eq. 2 .or. self%mphyopt .eq. 5 ) then
! Initialize coefficients and power index numers used in the dbz obs operator
! for single moment scheme (LIN)
    i_melt_snow=-1
    i_melt_graupel=100
    call coef4dbzfwrd(self%mphyopt,iret_coef4dbzfwrd)
     if ( iret_coef4dbzfwrd .ne. 0 ) then
        call abor1_ftn("iret_coef4dbzfwrd is not 0, stop")
     end if
  else if ( self%mphyopt .eq. 108 ) then
  ! Define Variables for dry Zs in TM
     oams = one/am_s
     alpha_const_tm_dry_snow =  (0.176_kind_real/0.93_kind_real) * (6.0_kind_real/pi)*(6.0_kind_real/pi)     &
                              *(am_s/900.0_kind_real)*(am_s/900.0_kind_real)*1.e18_kind_real
  end if

  do ivar = 1,1
    do iobs=1,nlocs 

      ! calculate RHO first
      if ( self%mphyopt .eq. 14 ) then
         idx_p = 11
         idx_t = 12
         idx_q = 13
      else if ( self%mphyopt .eq. 108 ) then
         idx_p = 5
         idx_t = 6
         idx_q = 7
      else if ( self%mphyopt .eq. 2 .or. self%mphyopt .eq. 5 ) then
         idx_p = 4
         idx_t = 5
         idx_q = 6
      end if
      P1D = vfields(idx_p,iobs) ! pressure (Pa)
      T1D = vfields(idx_t,iobs) ! temperature (K)
      Q1D = vfields(idx_q,iobs) ! specific humidity
      Q1D = Q1D/(one-Q1D)   ! convert from specific humidity to mixing ratio
      RHO = P1D/(rd*T1D*(one+D608*Q1D))

      ! Change MP variables from kg/kg to kg/m3, assume lower bound minimum values.
      qrges = max(vfields(1,iobs)*RHO, qx_min)  ! qr
      qsges = max(vfields(2,iobs)*RHO, qx_min)  ! qs
      qgges = max(vfields(3,iobs)*RHO, qx_min)  ! qg

      ! Some microphysics options have more species and double-moment.
      if ( self%mphyopt .eq. 108 ) then
         qnrges = max(vfields(4,iobs)*RHO, qn_min)  ! qnr
      elseif ( self%mphyopt .eq. 14 ) then
         qhges = max(vfields(4,iobs)*RHO, qx_min)   ! qh
         qnrges = max(vfields(5,iobs)*RHO, qn_min)  ! qnr
         qnsges = max(vfields(6,iobs)*RHO, qn_min)  ! qns
         qngges = max(vfields(7,iobs)*RHO, qn_min)  ! qng
         qnhges = max(vfields(8,iobs)*RHO, qn_min)  ! qnh
         qvgges = max(vfields(9,iobs)*RHO, qb_min)  ! qvg
         qvhges = max(vfields(10,iobs)*RHO, qb_min) ! qvh
      end if

      ! intialization
      Zer = zero
      Zes = zero
      Zeg = zero
      rdBZ = zero
      rdBZr = zero
      rdBZs = zero
      rdBZg = zero
      rdBZh = zero
      wgt_dry = zero
      wgt_wet = zero ;
      Zeg_dry = zero
      Zeg_wet = zero ;

! --------------- LIN operator code from GSI 'setupdbz.f90' --------------
      if ( self%mphyopt .eq. 2 .or. self%mphyopt .eq. 5 ) then ! LIN operator
         ! rain
         if (qrges > qx_min) then
            Zer = Cr * qrges**Pr
         end if

         ! snow
         if (qsges > qx_min) then
            if ( i_melt_snow < 0 ) then
               ! no melting: dry snow at any temperature
               Zes = Cs_dry * qsges**Ps_dry
            else if ( i_melt_snow  .eq. 100 ) then
               ! melting: wet snow at any temperature
               Zes = Cs_wet * qsges**Ps_wet
            else
               ! melting: depending on temperature
               if (T1D < T_melt) then
                  Zes = Cs_dry * qsges**Ps_dry
               else
                  Zes = Cs_wet * qsges**Ps_wet
               end if
            end if
         end if

         ! graupel/hail
         if (qgges > qx_min) then
            if ( i_melt_graupel < 0 ) then
               ! no melting: dry grauple/hail at any temperature
               Zeg = Cg_dry * qgges**Pg_dry
            else if ( i_melt_graupel  .eq. 100 ) then
               ! melting: wet graupel at any temperature
               Zeg = Cg_wet * qgges**Pg_wet
            else
               ! melting: depending on the temperature
               if (T1D < (T_melt - 2.5_kind_real)) then
                  Zeg = Cg_dry * qgges**Pg_dry
               else if (T1D > (T_melt + 2.5_kind_real)) then
                  Zeg = Cg_wet * qgges**Pg_wet
               else
                  wgt_dry = abs(T1D - (T_melt + 2.5_kind_real))/5.0_kind_real
                  wgt_wet = abs(T1D - (T_melt - 2.5_kind_real))/5.0_kind_real
                  Zeg_dry = Cg_dry * qgges**Pg_dry
                  Zeg_wet = Cg_wet * qgges**Pg_wet
                  Zeg     = wgt_dry*Zeg_dry + wgt_wet*Zeg_wet
               end if
            end if
         end if

         Ze=Zer+Zes+Zeg

!        Zelim treatment
         if (Ze < 1.0_kind_real) then
            Ze = Ze + 1.0_kind_real
         end if

!        Convert to simulated radar reflectivity in units of dBZ
         if (Ze > 0) rdBZ = ten * log10(Ze) ! hofx
         if (Zer > 0) rdBZr = ten * log10(Zer)
         if (Zes > 0) rdBZs = ten * log10(Zes)
         if (Zeg > 0) rdBZg = ten * log10(Zeg)

      end if

! --------------- TM operator code from GSI 'setupdbz.f90' --------------
      if ( self%mphyopt .eq. 108 .and. self%use_variational ) then ! TM operator
!----TM-------------
         Cg_wet= 5.54914E+12_kind_real
         Cg_dry= 1.90409E+12_kind_real
         Pg_wet= 2.5_kind_real

       ! rain
         if (qrges .gt. qx_min*100. .and. qnrges .ge. qn_min*100.) then
            Zer = 720._kind_real * qrges**2 / (pi**2*rhor**2*qnrges) * 1.E18
         endif

       ! snow
         if (qsges > qx_min*100.) then
             if ( T1D > T_melt ) then
                Zes = (1.47E+05_kind_real)*(qsges*1000._kind_real)**2.67_kind_real
             else
             ! calc coeffs for dry snow based on TM MP code 
                call calc_coeffs_dry_snow_tm(T1D,a_dry_snow_tm,b_dry_snow_tm)
                smoz = a_dry_snow_tm * (qsges*oams)**b_dry_snow_tm
                Zes  = alpha_const_tm_dry_snow*smoz
             end if
         endif

       ! graupel/hail
         if (qgges > qx_min*100.) then
             if ( T1D > T_melt ) then
                Zeg = Cg_wet*RHO**1.75_kind_real * qgges**Pg_wet
             else
                Zeg = Cg_dry*RHO**1.75_kind_real * qgges**Pg_wet
             end if
         endif

         Ze=Zer+Zes+Zeg

!      Zelim treatment
         if(Ze < 1.0_kind_real) then
              Ze=Ze + 1.0_kind_real
         end if

!      Convert to simulated radar reflectivity in units of dBZ
         rdBZ = ten * log10(Ze)
         if (Zer .gt. zero) rdBZr = ten * log10(Zer)
         if (Zes .gt. zero) rdBZs = ten * log10(Zes)
         if (Zeg .gt. zero) rdBZg = ten * log10(Zeg)

      end if

     !write(*,'(a,i5,3(a,f10.8,a,f6.2,a),a,f6.2,a)') 'iobs = ', iobs,  &
     !              ' Rain = ', qrges, ", ", rdBZr, " dBZ",            &
     !              ' Snow = ', qsges, ", ", rdBZs, " dBZ",            &
     !              ' Graupel = ', qgges, ", ", rdBZg, " dBZ",         &
     !              ' Total: ', rdBZ, " dBZ "

      if ( (self%mphyopt .eq. 108 .or. self%mphyopt .eq. 14 ) .and. &
           (.not.(self%use_variational)) ) then

! NSSL operator for EnKF / use radarZ
          if ( self%mphyopt .eq. 108 ) then
             if ( P_qr > 0 ) qscalar(P_qr) = qrges
             if ( P_qs > 0 ) qscalar(P_qs) = qsges
             if ( P_qg > 0 ) qscalar(P_qg) = qgges
             if ( P_nr > 0 ) qscalar(P_nr) = qnrges

          else if ( self%mphyopt .eq. 14 ) then
! NSSL operator for EnKF / use radarZ
             if ( P_qr > 0 ) qscalar(P_qr) = qrges
             if ( P_qs > 0 ) qscalar(P_qs) = qsges
             if ( P_qg > 0 ) qscalar(P_qg) = qgges
             if ( P_qh > 0 ) qscalar(P_qh) = qhges
             if ( P_nr > 0 ) qscalar(P_nr) = qnrges
             if ( P_ns > 0 ) qscalar(P_ns) = qnsges
             if ( P_ng > 0 ) qscalar(P_ng) = qngges
             if ( P_nh > 0 ) qscalar(P_nh) = qnhges
             if ( P_vg > 0 ) qscalar(P_vg) = qvgges
             if ( P_vh > 0 ) qscalar(P_vh) = qvhges
          end if

          call set_dsd_para()
          call calcMDR()
          call calcMu()

          obs_dual = init_refl()
          var_dsd = init_para_dsd()
          ta = T1D ! pass air temperature to radarZ
          call rdr_obs(real(RHO, kind=kind_real), real(qscalar, kind=kind_real),&
                            obs_dual, var_dsd, 1, 1)
          rdBZ = real(obs_dual%T_log_ref, kind=kind_real)
          Ze = 10._kind_real ** (rdBZ / 10._kind_real)
      end if

      hofx(ivar,iobs) = rdBZ ! assign hofx

    end do
  end do

  ! Cleanup memory
  deallocate(obsvcoord)
  deallocate(wi)
  deallocate(wf)
  
  deallocate(vfields)

end subroutine ufo_directZDA_simobs

end module ufo_directZDA_mod
