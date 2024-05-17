! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module for directZDA tl/ad observation operator

module ufo_directZDA_tlad_mod

 use oops_variables_mod
 use obs_variables_mod
 use ufo_vars_mod
 use ufo_geovals_mod
 use kinds
 use vert_interp_mod


 use ufo_directZDA_util_mod, only : coef4dbzfwrd, Cr, Pr, &
                                    Cs_dry, Cs_wet, &
                                    Cg_dry, Cg_wet, &
                                    Ps_dry, Ps_wet, Pg_dry, Pg_wet, &
                                    calc_coeffs_dry_snow_tm

 
 !> Fortran derived type for the tl/ad observation operator
 type, public  :: ufo_directZDA_tlad
 private
  type(obs_variables), public :: obsvars
  type(oops_variables), public :: geovars
  character(len=MAXVARLEN), public :: v_coord ! GeoVaL to use to interpolate in vertical
  character(len=MAXVARLEN), public :: micro_option    ! Choice (enum) of microphysics option
  integer, public :: mphyopt                  ! Integer value for microphysics option
  logical :: use_variational  = .false.
  integer :: nval, nlocs
  real(kind_real), allocatable  :: jac_qr(:), jac_qs(:), jac_qg(:), jac_qnr(:)
  real(kind_real),  allocatable :: wf(:)
  integer,          allocatable :: wi(:)
 contains
  procedure :: setup  => ufo_directZDA_tlad_setup_
  procedure :: cleanup => ufo_directZDA_tlad_cleanup_
  procedure :: settraj => ufo_directZDA_tlad_settraj_
  procedure :: simobs_tl  => ufo_directZDA_simobs_tl_
  procedure :: simobs_ad  => ufo_directZDA_simobs_ad_
  final :: destructor
 end type ufo_directZDA_tlad

 integer :: n_geovars
 character(len=maxvarlen), dimension(:), allocatable:: geovars_list
 character(len=maxvarlen):: varname_qr, varname_qs, varname_qg, varname_qnr

  integer, parameter         :: max_string=800

contains

! ------------------------------------------------------------------------------
subroutine ufo_directZDA_tlad_setup_(self, yaml_conf)
  use fckit_configuration_module, only: fckit_configuration
  use iso_c_binding
  use radarz_iface

  implicit none
  class(ufo_directZDA_tlad), intent(inout) :: self
  type(fckit_configuration), intent(in)  :: yaml_conf

  character(kind=c_char,len=:), allocatable :: coord_name
  character(kind=c_char,len=:), allocatable :: micro_option
  character(kind=c_char,len=:), allocatable :: this_varname
  character(len=maxvarlen):: var_string

  ! YAML option for microphysics scheme
  call yaml_conf%get_or_die("microphysics option", micro_option)
  if(trim(micro_option) .eq. "Thompson") then
    !self%micro_option = THOMPSON
    self%mphyopt = 108
    mp_option = 108
  else if(trim(micro_option) == 'Lin') then
    !self%micro_option = LIN
    self%mphyopt = 2
    mp_option = 2
  else if(trim(micro_option) == 'GFDL') then
    !self%micro_option = GFDL
    self%mphyopt = 5
    mp_option = 5
  else
    print*, ' TLAD micro_option set to ', trim(micro_option)
    call abor1_ftn("microphysics option not set or unsupported, aborting")
  endif

  ! If YAML option indicates use of variational method then update mp_option
  self%use_variational=.false.
  call yaml_conf%get_or_die("use variational method", self%use_variational)
  if (self%use_variational) then
    if (trim(micro_option) .ne. "Thompson" .and. trim(micro_option) .ne. "Lin" ) then
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

  if ( self%mphyopt .eq. 108 ) then      ! TM VarOP: 108
    n_geovars=7
    if (.not.allocated(geovars_list) ) allocate(geovars_list(n_geovars))
    var_string="var_rain_mixing_ratio"
    if( yaml_conf%has(trim(var_string)) ) then
       call yaml_conf%get_or_die(trim(var_string), this_varname)
       geovars_list(1) = this_varname
       varname_qr = this_varname
    endif
    var_string="var_snow_mixing_ratio"
    if( yaml_conf%has(trim(var_string)) ) then
       call yaml_conf%get_or_die(trim(var_string), this_varname)
       geovars_list(2) = this_varname
       varname_qs = this_varname
    endif
    var_string="var_graupel_mixing_ratio"
    if( yaml_conf%has(trim(var_string)) ) then
       call yaml_conf%get_or_die(trim(var_string), this_varname)
       geovars_list(3) = this_varname
       varname_qg = this_varname
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
       varname_qr = this_varname
    endif
    var_string="var_snow_mixing_ratio"
    if( yaml_conf%has(trim(var_string)) ) then
       call yaml_conf%get_or_die(trim(var_string), this_varname)
       geovars_list(2) = this_varname
       varname_qs = this_varname
    endif
    var_string="var_graupel_mixing_ratio"
    if( yaml_conf%has(trim(var_string)) ) then
       call yaml_conf%get_or_die(trim(var_string), this_varname)
       geovars_list(3) = this_varname
       varname_qg = this_varname
    endif
    geovars_list(4) = var_prs
    geovars_list(5) = var_ts
    geovars_list(6) = var_q
  end if

  call self%geovars%push_back(geovars_list)
  call self%geovars%push_back(self%v_coord)

end subroutine ufo_directZDA_tlad_setup_

! ------------------------------------------------------------------------------
subroutine ufo_directZDA_tlad_settraj_(self, geovals, obss)
  use obsspace_mod

  implicit none
  class(ufo_directZDA_tlad), intent(inout) :: self
  type(ufo_geovals),       intent(in)    :: geovals
  type(c_ptr), value,      intent(in)    :: obss

 ! Local variables
  integer :: iobs, ivar, nvars_geovars
  real(kind_real),  dimension(:), allocatable :: obsvcoord
  type(ufo_geoval), pointer :: vcoordprofile, profile

  real(kind_real), allocatable :: tmp(:)
  real(kind_real) :: tmp2
  real(kind_real), allocatable :: vfields(:,:)  ! background fields interplated vertically to the observation height

  integer :: i_melt_snow, i_melt_graupel
  real(kind_real) :: denom
  real(kind_real) :: RHO
  real(kind_real) :: qrexp, qsexp, qgexp, qnr1D
  real(kind_real) :: P1D, Q1D, T1D
  real(kind_real) :: Ze, Zer, Zes, Zeg
  real(kind_real) :: wgt_dry, wgt_wet
  real(kind_real) :: Zeg_dry, Zeg_wet
  real(kind_real) :: Cs_tmp, Cg_tmp
  real(kind_real) :: jqr_num,jqs_num,jqg_num, jqnr_num
  real(kind_real) :: jqg_num_dry, jqg_num_wet 
  real(kind_real) :: jqr, jqs, jqg, jqnr
  real(kind_real) :: qrges, qsges, qgges, qnrges
  real(kind_real) :: rdBZ, rdBZr, rdBZs, rdBZg

  real(kind_real), parameter :: qr_min = 1.0E-8_kind_real
  real(kind_real), parameter :: qs_min = 1.0E-8_kind_real
  real(kind_real), parameter :: qg_min = 1.0E-8_kind_real
  real(kind_real), parameter :: qh_min = 1.0E-8_kind_real
  real(kind_real), parameter :: qnr_min = 1.0E-8_kind_real

  real(kind_real), parameter :: rd=287.04_kind_real
  real(kind_real), parameter :: one=1._kind_real
  real(kind_real), parameter :: D608=0.608_kind_real
  real(kind_real), parameter :: zero=0.0_kind_real
  real(kind_real), parameter :: ten=10.0_kind_real
  real(kind_real), parameter :: T_melt=273.15_kind_real

  real(kind_real), parameter :: pi = 3.141592_kind_real   ! pi
  real(kind_real), parameter :: rhor=1000._kind_real      ! Density of rain (kg m**-3)
  real(kind_real), parameter :: rhoh=913._kind_real       ! Density of hail (kg m**-3)
  real(kind_real), parameter :: rhos=100._kind_real       ! Density of snow (kg m**-3)
  real(kind_real), parameter :: rhog=500._kind_real       ! Density of graupel (kg m**-3)
  real(kind_real), parameter :: am_s = 0.069_kind_real ! from WRF
  real(kind_real) :: a_dry_snow_tm, b_dry_snow_tm ! dry snow coeffs for TM
  real(kind_real) :: smoz, alpha_const_tm_dry_snow
  real(kind_real) :: oams

  real(kind_real) :: Ze_orig
  integer :: iret_coef4dbzfwrd

  type(ufo_geoval),    pointer :: qr, qs, qg, qnr, prs, t, qv

  ! Make sure nothing already allocated
  call self%cleanup()

  ! SetTraj - calc jacobian
  ! TL/AD - only load GEOVALS for QR/QS/QG and QNR and interpolate to VFIELD
  ! Let's see if this eliminate the use of analysis/state variable..

! ---  get vcoord
  call ufo_geovals_get_var(geovals, self%v_coord, vcoordprofile)

! Keep copy of dimensions
  self%nval  = vcoordprofile%nval
  self%nlocs = obsspace_get_nlocs(obss)

! alllocate arrays for jacobian/interpolation weights at TLAD
  allocate(self%jac_qr(self%nlocs))
  allocate(self%jac_qs(self%nlocs))
  allocate(self%jac_qg(self%nlocs))
  if ( self%mphyopt .eq. 108 ) then ! TM operator only
      allocate(self%jac_qnr(self%nlocs))
  end if
  allocate(self%wi(self%nlocs))
  allocate(self%wf(self%nlocs))

! get geoval column
  call ufo_geovals_get_var(geovals, varname_qr, qr)
  call ufo_geovals_get_var(geovals, varname_qs, qs)
  call ufo_geovals_get_var(geovals, varname_qg, qg)
  call ufo_geovals_get_var(geovals, var_prs, prs)
  call ufo_geovals_get_var(geovals, var_ts, t)
  call ufo_geovals_get_var(geovals, var_q, qv)
  if ( self%mphyopt .eq. 108 ) then ! TM operator only
     call ufo_geovals_get_var(geovals, varname_qnr, qnr)
  end if

  allocate(obsvcoord(self%nlocs))
  call obsspace_get_db(obss, "MetaData", "height", obsvcoord)

! Calculate the interpolation weights (wi, wf)
  allocate(tmp(vcoordprofile%nval))
  do iobs = 1, self%nlocs
    tmp = vcoordprofile%vals(:,iobs)
    tmp2 = obsvcoord(iobs)
    call vert_interp_weights(vcoordprofile%nval, tmp2, tmp, self%wi(iobs), self%wf(iobs))
  end do

  if ( self%mphyopt .eq. 2 .or. self%mphyopt .eq. 5 ) then
! Initialize coefficients and power index numers used in the dbz obs operator
! for single moment scheme
    i_melt_snow=-1
    i_melt_graupel=100
    call coef4dbzfwrd(self%mphyopt,iret_coef4dbzfwrd)
  else if ( self%mphyopt .eq. 108 ) then
  ! Define Variables for dry Zs in TM
    oams = one/am_s
    alpha_const_tm_dry_snow =  (0.176_kind_real/0.93_kind_real) * (6.0_kind_real/pi)*(6.0_kind_real/pi)     &
                                *(am_s/900.0_kind_real)*(am_s/900.0_kind_real)*1.e18_kind_real
  end if

  do iobs = 1, self%nlocs
    call vert_interp_apply(qr%nval, qr%vals(:,iobs), qrges, self%wi(iobs),self%wf(iobs))
    call vert_interp_apply(qs%nval, qs%vals(:,iobs), qsges, self%wi(iobs),self%wf(iobs))
    call vert_interp_apply(qg%nval, qg%vals(:,iobs), qgges, self%wi(iobs),self%wf(iobs))
    call vert_interp_apply(prs%nval, prs%vals(:,iobs), P1D, self%wi(iobs),self%wf(iobs))
    call vert_interp_apply(t%nval, t%vals(:,iobs), T1D, self%wi(iobs),self%wf(iobs))
    call vert_interp_apply(qv%nval, qv%vals(:,iobs), Q1D, self%wi(iobs),self%wf(iobs))
    if ( self%mphyopt .eq. 108 ) then ! TM operator only
       call vert_interp_apply(qnr%nval, qnr%vals(:,iobs), qnrges, self%wi(iobs),self%wf(iobs))
    end if

! Treat minimum Qx in case of variational DA
    if ( self%use_variational ) then
       qrexp = max(qrges, qr_min)
       qsexp = max(qsges, qs_min)
       qgexp = max(qgges, qg_min)
       if ( self%mphyopt .eq. 108 ) qnr1D=max(qnrges, qnr_min)
    else
       qrexp=qrges
       qsexp=qsges
       qgexp=qgges
       if ( self%mphyopt .eq. 108 ) qnr1D=qnrges
    end if
!
    Q1D=Q1D/(one-Q1D)   ! convert to mixing ratio
    RHO=P1D/(rd*T1D*(one+D608*Q1D))

    if ( self%mphyopt  .eq. 108 ) then ! TM operator
        qnr1D=qnr1D*RHO ! convert qnr to unit demanded by operators
    end if

    ! intialization
    Zer = zero ;         Zes = zero ;          Zeg = zero ;
    wgt_dry = zero ;     wgt_wet = zero ;
    Zeg_dry = zero ;     Zeg_wet = zero ;

! --------------- LIN operator code from GSI 'setupdbz.f90' --------------
    if ( self%mphyopt .eq. 2 .or. self%mphyopt .eq. 5 ) then ! LIN operator
      ! rain
        Zer = Cr  * (RHO * qrexp)**(Pr)

      ! snow
        if ( i_melt_snow < 0 ) then
!       no melting: dry snow at any temperature
           Zes = Cs_dry * (RHO * qsexp)**(Ps_dry)
           Cs_tmp  = Cs_dry
        else if ( i_melt_snow  .eq. 100 ) then
!       melting: wet snow at any temperature
           Zes = Cs_wet * (RHO * qsexp)**(Ps_wet)
           Cs_tmp  = Cs_wet
        else
!       melting: depending on temperature
           if (T1D < T_melt) then
              Zes = Cs_dry * (RHO * qsexp)**(Ps_dry)
              Cs_tmp  = Cs_dry
           else
              Zes = Cs_wet * (RHO * qsexp)**(Ps_wet)
              Cs_tmp  = Cs_wet
           end if
        end if

      ! graupel/hail
        if ( i_melt_graupel < 0 ) then
!       no melting: dry grauple/hail at any temperature
           Zeg = Cg_dry * (RHO * qgexp)**(Pg_dry)
           Cg_tmp  = Cg_dry
        else if ( i_melt_graupel  .eq. 100 ) then
!       melting: wet graupel at any temperature
           Zeg = Cg_wet * (RHO * qgexp)**(Pg_wet)
           Cg_tmp  = Cg_wet
        else
!       melting: depending on the temperature
           if (T1D < (T_melt - 2.5_kind_real)) then
              Zeg = Cg_dry * (RHO * qgexp)**(Pg_dry)
              Cg_tmp  = Cg_dry
           else if (T1D > (T_melt + 2.5_kind_real)) then
              Zeg = Cg_wet * (RHO * qgexp)**(Pg_wet)
              Cg_tmp  = Cg_wet
           else
              wgt_dry = abs(T1D - (T_melt + 2.5_kind_real))/5.0_kind_real
              wgt_wet = abs(T1D - (T_melt - 2.5_kind_real))/5.0_kind_real
              Zeg_dry = Cg_dry * (RHO * qgexp)**(Pg_dry)
              Zeg_wet = Cg_wet * (RHO * qgexp)**(Pg_wet)
              Zeg     = wgt_dry*Zeg_dry + wgt_wet*Zeg_wet
              Cg_tmp  = wgt_dry*Cg_dry  + wgt_wet*Cg_wet
           end if
        end if

        Ze=Zer+Zes+Zeg

!       Zelim treatment
        if (Ze < 1.0_kind_real) then
           Ze=Ze + 1.0_kind_real
        end if

      ! Convert to simulated radar reflectivity in units of dBZ
        Ze_orig = Ze
        rdBZ = ten * log10(Ze)
        rdBZr = ten * log10(Zer)
        rdBZs = ten * log10(Zes)
        rdBZg = ten * log10(Zeg)

! ----- CALCULATE JACOBIAN ----
!   find dqr/ddBZ, dqs/ddBZ, dqg/ddBZ (used in inner loop routine)
!   Jacobian used for TLM and ADM
!       rain
        jqr_num = ten*Cr*((RHO)**Pr)*Pr*((qrexp)**(Pr - one))

!       snow
        if ( i_melt_snow < 0 ) then
!          no melting: dry snow at any temperature
           jqs_num = ten*Cs_dry*((RHO)**Ps_dry)*Ps_dry*((qsexp)**(Ps_dry - one))
        else if ( i_melt_snow  .eq. 100 ) then
!          melting: wet snow at any temperature
           jqs_num = ten*Cs_wet*((RHO)**Ps_wet)*Ps_wet*((qsexp)**(Ps_wet - one))
        else
!          melting: depending on temperature
           if (T1D < T_melt) then
              jqs_num = ten*Cs_dry*((RHO)**Ps_dry)*Ps_dry*((qsexp)**(Ps_dry - one))
           else
              jqs_num = ten*Cs_wet*((RHO)**Ps_wet)*Ps_wet*((qsexp)**(Ps_wet - one))
           end if
        end if

!       graupel/hail
        if ( i_melt_graupel < 0 ) then
!          no melting: dry grauple/hail at any temperature
           jqg_num = ten*Cg_dry*((RHO)**(Pg_dry))*(Pg_dry)*((qgexp**(Pg_dry - one)))
        else if ( i_melt_graupel  .eq. 100 ) then
!          melting: wet graupel at any temperature
           jqg_num = ten*Cg_wet*((RHO)**(Pg_wet))*(Pg_wet)*((qgexp**(Pg_wet - one)))
        else
!          melting: depending on the temperature
           if (T1D < (T_melt - 2.5_kind_real)) then
              jqg_num = ten*Cg_dry*((RHO)**(Pg_dry))*(Pg_dry)*((qgexp**(Pg_dry - one)))
           else if (T1D > (T_melt + 2.5_kind_real)) then
              jqg_num = ten*Cg_wet*((RHO)**(Pg_wet))*(Pg_wet)*((qgexp**(Pg_wet - one)))
           else
              wgt_dry = abs(T1D - (T_melt + 2.5_kind_real))/5.0_kind_real
              wgt_wet = abs(T1D - (T_melt - 2.5_kind_real))/5.0_kind_real
              jqg_num_dry = ten*Cg_dry*((RHO)**(Pg_dry))*(Pg_dry)*((qgexp**(Pg_dry - one)))
              jqg_num_wet = ten*Cg_wet*((RHO)**(Pg_wet))*(Pg_wet)*((qgexp**(Pg_wet - one)))
              jqg_num = wgt_dry*jqg_num_dry + wgt_wet*jqg_num_wet
           end if
        end if

        denom=(log(ten))*Ze

        jqr  = jqr_num/denom
        jqs  = jqs_num/denom
        jqg  = jqg_num/denom

        self%jac_qr(iobs)=jqr
        self%jac_qs(iobs)=jqs
        self%jac_qg(iobs)=jqg

    else if ( self%mphyopt  .eq. 108 ) then ! TM operator
!----TM-------------
         Cg_wet= 5.54914E+12_kind_real ! wet Zeg
         Cg_dry= 1.90409E+12_kind_real ! dry Zeg
         Pg_wet= 2.5_kind_real         ! used for both dry and wet

       ! rain
         if (  qrexp > 10E-7_kind_real .and. qnr1D > 20._kind_real  ) then
            Zer = 720._kind_real *(RHO*qrexp)**2*ten**18/(pi**2*rhor**2*qnr1D)
         else
            Zer = 10E-8_kind_real
         endif

       ! snow
         if (qsexp > 10E-7_kind_real) then
            if ( T1D > T_melt ) then
                Zes = (1.47E+05_kind_real)*(qsexp*1000._kind_real)**2.67_kind_real
            else
             ! calc coeffs for dry snow based on TM MP code 
                call calc_coeffs_dry_snow_tm(T1D,a_dry_snow_tm,b_dry_snow_tm)
                smoz = a_dry_snow_tm * (qsexp*RHO*oams)**b_dry_snow_tm
                Zes  = alpha_const_tm_dry_snow*smoz
            end if
         else
            Zes = 10E-8_kind_real
         endif

       ! graupel/hail
         if ( qgexp > 10E-7_kind_real ) then
             if ( T1D > T_melt ) then
                Zeg = Cg_wet*RHO**1.75_kind_real*qgexp**Pg_wet
             else
                Zeg = Cg_dry*RHO**1.75_kind_real*qgexp**Pg_wet
             end if
         else
             Zeg=10E-8_kind_real
         endif

         Ze=Zer+Zes+Zeg

    !    Zelim treatment
         if (Ze < 1.0_kind_real) then
            Ze=Ze + 1.0_kind_real
         end if

!        Convert to simulated radar reflectivity in units of dBZ
         rdBZ = ten * log10(Ze)
         rdBZr = ten * log10(Zer)
         rdBZs = ten * log10(Zes)
         rdBZg = ten * log10(Zeg)

!       find dqr/ddBZ, dqs/ddBZ, dqg/ddBZ (used in inner loop routine)
!       Jacobian used for TLM and ADM
!        rain
         if ( qrexp > 10E-7_kind_real .and. qnr1D > 1._kind_real ) then
             jqr_num = (2_kind_real*720_kind_real*RHO**2_kind_real*qrexp/&
                       (pi**2_kind_real*rhor**2_kind_real*qnr1D))*ten**12_kind_real*ten**6_kind_real
             jqnr_num = (-720_kind_real*RHO**2_kind_real*qrexp**2_kind_real/&
                       (pi**2_kind_real*rhor**2_kind_real*qnr1D**2_kind_real))*ten**12_kind_real*ten &
                       **6_kind_real
         else
             jqr_num = 10E-8_kind_real
             jqnr_num = 10E-8_kind_real
         endif

!        snow
         if (qsexp > 10E-7_kind_real) then
            if ( T1D > T_melt ) then
               jqs_num= (1.47E+05_kind_real)*(1000._kind_real)**2.67_kind_real*(2.67_kind_real) &
                        *(qsexp)**(2.67_kind_real-1.0_kind_real)
            else
               jqs_num= alpha_const_tm_dry_snow*a_dry_snow_tm*b_dry_snow_tm &
                        *((RHO)**b_dry_snow_tm)*((oams)**b_dry_snow_tm)*((qsexp) &
                        **(b_dry_snow_tm-1.0_kind_real))
            end if
         else
            jqs_num = 10E-8_kind_real
         endif

!        graupel/hail
         if ( qgexp > 10E-7_kind_real ) then
            if ( T1D > T_melt ) then
               jqg_num=Cg_wet*((RHO)**1.75_kind_real)*Pg_wet*((qgexp)**(Pg_wet-1.0_kind_real))
            else
               jqg_num=Cg_dry*((RHO)**1.75_kind_real)*Pg_wet*((qgexp)**(Pg_wet-1.0_kind_real))
            end if

         else
            jqg_num=10E-8_kind_real
         endif

         denom=(log(ten))*Ze

       ! multiply ten for numerators
         jqr  = ten*jqr_num/denom
         jqs  = ten*jqs_num/denom
         jqg  = ten*jqg_num/denom
         jqnr  = ten*jqnr_num/denom

         self%jac_qr(iobs)=jqr
         self%jac_qs(iobs)=jqs
         self%jac_qg(iobs)=jqg
         self%jac_qnr(iobs)=jqnr

    end if

  end do

! cleanup
  deallocate(obsvcoord)
  deallocate(tmp)

end subroutine ufo_directZDA_tlad_settraj_

! ------------------------------------------------------------------------------
! Input geovals parameter represents dx for tangent linear model
subroutine ufo_directZDA_simobs_tl_(self, geovals, obss, nvars, nlocs, hofx)

  implicit none
  class(ufo_directZDA_tlad), intent(in)    :: self
  type(ufo_geovals),       intent(in)    :: geovals
  integer,                   intent(in) :: nvars, nlocs
  real(c_double),         intent(inout) :: hofx(nvars, nlocs)
  type(c_ptr), value,        intent(in) :: obss

  character(len=*), parameter :: myname="ufo_directZDA_tlad_tl"
  character(max_string)       :: err_msg

  integer :: iobs, ivar, nvars_geovars
  type(ufo_geoval), pointer :: profile
  character(len=MAXVARLEN) :: geovar
  character(len=maxvarlen), dimension(:), allocatable :: geovars_list

  real(kind_real), allocatable :: vfields(:,:)  ! background fields interplated vertically to the observation height

! Just to use related geovals here
! Number of variables in geovars (without the vertical coordinate)
  if ( self%mphyopt .eq. 108 ) then ! TM
     nvars_geovars = 4 ! Qr, Qs, Qg, Qnr
     allocate(geovars_list(nvars_geovars))
     geovars_list=(/varname_qr, varname_qs, varname_qg, varname_qnr/)
  else                              ! LIN 
     nvars_geovars = 3 ! Qr, Qs, Qg
     allocate(geovars_list(nvars_geovars))
     geovars_list=(/varname_qr, varname_qs, varname_qg/)
  end if
  allocate(vfields(nvars_geovars,nlocs))
  vfields=0.0_kind_real

  do ivar = 1, nvars_geovars
    geovar = geovars_list(ivar)

! Get profile for this variable from geovals
    call ufo_geovals_get_var(geovals, geovar, profile)

! Interpolate from geovals to observational location into hofx
    do iobs = 1, nlocs
      call vert_interp_apply_tl(profile%nval, profile%vals(:,iobs), &
                             & vfields(ivar,iobs), self%wi(iobs), self%wf(iobs))
    enddo
  enddo

  do ivar = 1, nvars
    do iobs=1,nlocs
      hofx(ivar,iobs) = vfields(1,iobs)*self%jac_qr(iobs) &
                      + vfields(2,iobs)*self%jac_qs(iobs) &
                      + vfields(3,iobs)*self%jac_qg(iobs)
      if ( self%mphyopt .eq. 108 ) then ! add qnr as well
         hofx(ivar,iobs) = hofx(ivar,iobs) +  vfields(4,iobs)*self%jac_qnr(iobs) 
      end if
    enddo
  end do

  deallocate(vfields)

end subroutine ufo_directZDA_simobs_tl_

! ------------------------------------------------------------------------------
subroutine ufo_directZDA_simobs_ad_(self, geovals, obss, nvars, nlocs, hofx)

  implicit none
  class(ufo_directZDA_tlad), intent(in)    :: self
  type(ufo_geovals),       intent(inout) :: geovals
  integer,                   intent(in)    :: nvars, nlocs
  real(c_double),            intent(in)    :: hofx(nvars, nlocs)
  type(c_ptr), value,        intent(in)    :: obss

  character(len=*), parameter :: myname="ufo_directZDA_tlad_ad"
  character(max_string)       :: err_msg

  integer :: iobs, ivar, nvars_geovars
  type(ufo_geoval), pointer :: profile
  character(len=MAXVARLEN) :: geovar
  character(len=maxvarlen), dimension(:), allocatable :: geovars_list

  real(c_double) :: missing
  real(kind_real), allocatable :: vfields(:,:)  ! background fields interplated vertically to the observation height

    missing = missing_value(missing)

! Just to use related geovals only here
! Number of variables in geovars (without the vertical coordinate)
  if ( self%mphyopt .eq. 108 ) then ! TM
     nvars_geovars = 4 ! Qr, Qs, Qg, Qnr
     allocate(geovars_list(nvars_geovars))
     geovars_list=(/varname_qr, varname_qs, varname_qg, varname_qnr/)
  else                              ! LIN 
     nvars_geovars = 3 ! Qr, Qs, Qg
     allocate(geovars_list(nvars_geovars))
     geovars_list=(/varname_qr, varname_qs, varname_qg/)
  end if

  allocate(vfields(nvars_geovars,nlocs))

  vfields=0.0_kind_real

  do ivar = 1, nvars
    do iobs=1,nlocs
     if (hofx(ivar,iobs) .ne. missing) then
       vfields(1,iobs) = vfields(1,iobs) + hofx(ivar,iobs)*self%jac_qr(iobs)
       vfields(2,iobs) = vfields(2,iobs) + hofx(ivar,iobs)*self%jac_qs(iobs)
       vfields(3,iobs) = vfields(3,iobs) + hofx(ivar,iobs)*self%jac_qg(iobs)
       if ( self%mphyopt .eq. 108 ) then ! TM
          vfields(4,iobs) = vfields(4,iobs) + hofx(ivar,iobs)*self%jac_qnr(iobs)
        end if
     end if
    enddo
  end do

  do ivar = 1, nvars_geovars
! Get the name of input variable in geovals
    geovar = geovars_list(ivar)

! Get profile for this variable from geovals
    call ufo_geovals_get_var(geovals, geovar, profile)

! Interpolate from geovals to observational location into hofx
    do iobs = 1, nlocs
      call vert_interp_apply_ad(profile%nval, profile%vals(:,iobs), &
                             & vfields(ivar,iobs), self%wi(iobs), self%wf(iobs))
    enddo
  enddo

  deallocate(vfields)

end subroutine ufo_directZDA_simobs_ad_

! ------------------------------------------------------------------------------
subroutine ufo_directZDA_tlad_cleanup_(self)
  implicit none
  class(ufo_directZDA_tlad), intent(inout) :: self

  self%nval  = 0
  self%nlocs = 0
  if (allocated(self%jac_qr)) deallocate(self%jac_qr)
  if (allocated(self%jac_qs)) deallocate(self%jac_qs)
  if (allocated(self%jac_qg)) deallocate(self%jac_qg)
  if (allocated(self%jac_qnr)) deallocate(self%jac_qnr)
  if (allocated(self%wf)) deallocate(self%wf)
  if (allocated(self%wi)) deallocate(self%wi)

end subroutine ufo_directZDA_tlad_cleanup_


! ------------------------------------------------------------------------------

subroutine  destructor(self)
  type(ufo_directZDA_tlad), intent(inout)  :: self

  call self%cleanup()

end subroutine destructor

! ------------------------------------------------------------------------------


end module ufo_directZDA_tlad_mod

