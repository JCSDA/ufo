module steric

  implicit none
  private
  public :: steric_nl, steric_tl, steric_ad
contains
  !==========================================================================
  subroutine steric_nl (etas, t, s, t0, s0, eta0, z)
    !==========================================================================
    !
    ! Calculates the steric height relative to the reference state t0, s0, eta0
    !
    ! Input:
    ! ------
    ! s      : Absolute Salinity                               [g/kg]
    ! t      : Conservative Temperature                        [deg C]
    ! s0     : Ref. Absolute Salinity                          [g/kg]
    ! t0     : Ref. Conservative Temperature                   [deg C]
    ! eta0   : Reference sea-surface height                    [m]
    ! z      : Depth of t,s                                    [m]
    !
    ! Output:
    ! -------
    ! etas   : steric height                                   [m]
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

    use gsw_mod_toolbox, only : gsw_rho
    use gsw_mod_kinds

    implicit none

    real (r8), dimension(:), intent(in) :: t, s, t0, s0, eta0, z
    real (r8), intent(out) :: etas(1)
    real (r8) :: rho0(1), rho(1), h
    real (r8) :: p(1), lat=0.0

    integer :: k, N

    N=size(t,1)
    etas=eta0
    do k=1,N
       h=z(k+1)-z(k)
       p = z(k)+h/2.0 ! assume p~z
       rho0 = gsw_rho(s0(k),t0(k),p)
       rho = gsw_rho(s(k),t(k),p)
       etas = etas + h*(rho-rho0)/rho0
    end do
    
  end subroutine steric_nl

  !==========================================================================
  subroutine steric_tl (detas, dt, ds, t0, s0, z)
    !==========================================================================
    !
    ! Tangent of stericnl relative to the reference state t0, s0
    !
    ! Input:
    ! ------
    ! ds     : Absolute Salinity                               [g/kg]
    ! dt     : Conservative Temperature                        [deg C]
    ! s0     : Ref. Absolute Salinity                          [g/kg]
    ! t0     : Ref. Conservative Temperature                   [deg C]
    ! z      : Depth of t,s                                    [m]
    !
    ! Output:
    ! -------
    ! detas  : steric height                                   [m]
    !
    ! jac    : Jacobian [detas/dt1, ...,detas/dtN;             [m/deg C]  
    !                    detas/ds1, ...,detas/dsN]             [m/(g/kg)]
    !    
    !--------------------------------------------------------------------------

    use gsw_mod_toolbox, only : gsw_rho, gsw_rho_first_derivatives
    use gsw_mod_kinds

    implicit none

    real (r8), dimension(:), intent(in) :: dt, ds, t0, s0, z
    real (r8), intent(out) :: detas
    real (r8), allocatable, dimension(:,:) :: jac 
    real (r8) :: rho0(1), rho(1), h, deriv(3)
    real (r8) :: p(1), drhods, drhodt, drhodp

    integer :: k, N

    detas = 0.0
    N=size(dt,1)
    allocate(jac(2,N)) !jac(1,:)=deta/dt; jac(2,:)=deta/ds at t0, s0
    do k=1,N
       h=z(k+1)-z(k)
       p = z(k)+h/2.0 ! assume p~z
       rho0 = gsw_rho(s0(k),t0(k),p)
       call gsw_rho_first_derivatives(s0(k),t0(k),p(1),drhods, drhodt, drhodp)
       jac(1,k)=h*drhodt/rho0(1)
       jac(2,k)=h*drhods/rho0(1)
       detas = detas + jac(1,k)*dt(k) + jac(2,k)*ds(k)
    end do
    deallocate(jac)
  end subroutine steric_tl

  !==========================================================================
  subroutine steric_ad (detas, dt_ad, ds_ad, t0, s0, z)
    !==========================================================================
    !
    ! Adjoint of sterictl relative to the trajectory t0, s0
    !
    ! Input:
    ! ------
    ! detas  : Steric Height                                   [m]    
    ! s0     : Traj. for Absolute Salinity                     [g/kg]
    ! t0     : Traj. for Conservative Temperature              [deg C]
    ! z      : Depth of t,s                                    [m]
    !
    ! Output:
    ! -------
    !
    ! ds_ad  : Adjoint var for Absolute Salinity               [g/kg]
    ! dt_ad  : Adjoint var for Conservative Temperature        [deg C]
    !
    !--------------------------------------------------------------------------

    use gsw_mod_toolbox, only : gsw_rho, gsw_rho_first_derivatives
    use gsw_mod_kinds

    implicit none

    real (r8), dimension(:), intent(in) :: t0, s0, z
    real (r8), intent(in) :: detas
    real (r8), allocatable, dimension(:,:) :: jac 
    real (r8), dimension(:), intent(out) :: dt_ad, ds_ad

    real (r8) :: rho0(1), rho(1), h, deriv(3)
    real (r8) :: p(1), drhods, drhodt, drhodp

    integer :: k, N

    N=size(dt_ad,1)
    allocate(jac(2,N)) !jac(1,:)=deta/dt; jac(2,:)=deta/ds at t0, s0
    dt_ad = 0.0
    ds_ad = 0.0    
    do k=1,N
       h=z(k+1)-z(k)
       p = z(k)+h/2.0 ! assume p~z
       rho0 = gsw_rho(s0(k),t0(k),p)
       call gsw_rho_first_derivatives(s0(k),t0(k),p(1),drhods, drhodt, drhodp)
       jac(1,k)=h*drhodt/rho0(1)
       jac(2,k)=h*drhods/rho0(1)
       dt_ad(k) = dt_ad(k) + jac(1,k)*detas
       ds_ad(k) = ds_ad(k) + jac(2,k)*detas
    end do
    deallocate(jac)
  end subroutine steric_ad

end module steric
