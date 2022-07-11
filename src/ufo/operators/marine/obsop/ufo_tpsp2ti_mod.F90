module ufo_tpsp2ti_mod
  implicit none
  private
  public :: insitu_t_nl, insitu_t_tl, insitu_t_tlad, insitu_t_jac
contains
  subroutine insitu_t_nl(temp_i, temp_p, salt_p, lono, lato, deptho)
    !==========================================================================
    !
    ! Calculates insitu temperature from potential temp and practical salinity
    !
    ! Input:
    ! ------
    ! temp_p : Potential teperature                            [deg C]
    ! salt_p : Practical salinity                              [psu]
    ! layer  : Layer thickness                                 [m]
    ! lono   : Ref. Conservative Temperature                   longitude of obs
    ! lato   : Reference sea-surface height                    latitude    "
    ! deptho : Depth of t,s                                    depth       "
    !
    ! Output:
    ! -------
    ! temp_i   : Insitu temperature interpolated at obs depth location  [C]
    ! ti       : Insitu temperature at model levels                     [C]    
    !--------------------------------------------------------------------------
    !
    ! z1=0 ----------- Surface
    !         t1,s1      h1 (thickness)  
    !  z2  -----------
    !         t2,s2      h2
    !  z3  -----------
    !         t3,s3
    !           .
    !           .    
    !           .
    !  zN-1  -----------
    !         tN-1,sN-1  hN-1
    !  zN  -----------
    !         tN,sN      hN
    ! zN+1 ----------- Bottom
    !      ///////////
    
    use gsw_mod_kinds
    use gsw_pot_to_insitu
    use vert_interp_mod
    use kinds
    use ieee_arithmetic
    
    implicit none

    real(kind=kind_real), intent(in)    :: temp_p             !< Potential temperature at observation location [C]
    real(kind=kind_real), intent(in)    :: salt_p             !< Practical Salinity at observation location [ppt]
    real(kind=kind_real), intent(in)    :: lono, lato, deptho !< Observation location    
    real(kind=kind_real), intent(inout) :: temp_i             !< Vertically interpolated Insitu temperature at

    real(kind_real) :: tp, sp, prso


    !< Pressure from depth
    prso = p_from_z(deptho,lato)

    ! Insitu temperature
    temp_i = t_from_pt(temp_p,salt_p,prso,lono,lato)

    if (ieee_is_nan(temp_i)) temp_i = 0.0
    
  end subroutine insitu_t_nl

  subroutine insitu_t_tl(dtemp_i, dtemp_p, dsalt_p, temp_p, salt_p, lono, lato, deptho, Jacobian)
    
    use gsw_mod_kinds
    use gsw_pot_to_insitu
    use vert_interp_mod
    use kinds
    use ieee_arithmetic

    implicit none

    real(kind=kind_real), intent(in)    :: dtemp_p            !< Potential temperature at observation location [C]
    real(kind=kind_real), intent(in)    :: dsalt_p            !< Practical Salinity at observation location [ppt]
    real(kind=kind_real), intent(in)    :: temp_p             !< Bkg potential temperature at observation location [C]
    real(kind=kind_real), intent(in)    :: salt_p             !< Bkg practical Salinity at observation location [ppt]    
    real(kind=kind_real), intent(in)    :: lono, lato, deptho !< Observation location
    real(kind=kind_real), intent(inout) :: dtemp_i            !< Vertically interpolated Insitu temperature at
    real(kind=kind_real), intent(in), optional :: Jacobian(2) !< Precomputed Jacobian    

    real(kind=kind_real) :: prso, jac(2)    

    !< Pressure from depth
    prso = p_from_z(deptho,lato)

    ! Jacobian: dTi/dTp and dTi/dSp
    if (present(Jacobian)) then
       jac=Jacobian
    else
       call insitu_t_jac(jac, temp_p, salt_p, lono, lato, deptho)
    end if
    
    ! Tangent of insitu temperature at (tp,sp)
    dtemp_i = jac(1)*dtemp_p + jac(2)*dsalt_p

    if (ieee_is_nan(dtemp_i)) dtemp_i = 0.0
    
  end subroutine insitu_t_tl

  subroutine insitu_t_tlad(dtemp_i, dtemp_p, dsalt_p, temp_p, salt_p, lono, lato, deptho, Jacobian)
    
    use gsw_mod_kinds
    use gsw_pot_to_insitu
    use vert_interp_mod
    use kinds
    use ieee_arithmetic

    implicit none

    real(kind=kind_real), intent(inout) :: dtemp_p            !< Potential temperature at observation location [C]
    real(kind=kind_real), intent(inout) :: dsalt_p            !< Practical Salinity at observation location [ppt]
    real(kind=kind_real), intent(in)    :: temp_p             !< Bkg potential temperature at observation location [C]
    real(kind=kind_real), intent(in)    :: salt_p             !< Bkg practical Salinity at observation location [ppt]    
    real(kind=kind_real), intent(in)    :: lono, lato, deptho !< Observation location
    real(kind=kind_real), intent(in)    :: dtemp_i            !< Insitu temperature increment at obs loc
    real(kind=kind_real), intent(in), optional :: Jacobian(2) !< Precomputed Jacobian
    
    real(kind=kind_real) :: prso, jac(2)    

    !< Pressure from depth
    prso = p_from_z(deptho,lato)

    ! Jacobian: dTi/dTp and dTi/dSp
    if (present(Jacobian)) then
       jac=Jacobian
    else
       call insitu_t_jac(jac, temp_p, salt_p, lono, lato, deptho)
    end if

    !< Adjoint
    dtemp_p = dtemp_p + jac(1)*dtemp_i
    dsalt_p = dsalt_p + jac(2)*dtemp_i    

    if (ieee_is_nan(dtemp_p+dsalt_p)) then
       dtemp_p = 0.0
       dsalt_p = 0.0
    end if
    
  end subroutine insitu_t_tlad
  
  subroutine insitu_t_jac(jac, temp_p, salt_p, lono, lato, deptho)
    !==========================================================================
    ! return jacobian at model levels
    use gsw_mod_kinds
    use gsw_pot_to_insitu
    use vert_interp_mod
    use kinds
    use ieee_arithmetic

    implicit none

    real(kind=kind_real), intent(in)    :: temp_p             !< Potential temperature at observation location [C]
    real(kind=kind_real), intent(in)    :: salt_p             !< Practical Salinity at observation location [ppt]
    real(kind=kind_real), intent(in)    :: lono, lato, deptho !< Observation location    
    real(kind=kind_real), intent(out)   :: jac(2)             !< Jacobian (dti/dtp, dti,dsp)
    
    real(kind=kind_real) :: pressure
    real(kind=kind_real) :: delta=1.0e-10        
    real(kind=kind_real) :: delta_tp, delta_sp    

    ! Vertical interpolation
    real(kind_real) :: wf
    integer :: wi

    delta_tp=delta
    delta_sp=delta
    
    !Pressure from depth
    pressure=p_from_z(deptho,lato)

    !Jacobian of Insitu temperature
    !dti/dtp
    Jac(1) = ( t_from_pt(temp_p+delta_tp,salt_p,pressure,lono,lato) -&
              &t_from_pt(temp_p,salt_p,pressure,lono,lato) )/delta_tp
    !dti/dsp
    Jac(2) = ( t_from_pt(temp_p,salt_p+delta_sp,pressure,lono,lato) -&
              &t_from_pt(temp_p,salt_p,pressure,lono,lato) )/delta_tp       

    !if (ieee_is_nan(sum(Jac))) Jac = 0.0
    
  end subroutine insitu_t_jac

  
end module ufo_tpsp2ti_mod
