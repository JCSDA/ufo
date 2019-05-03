module ufo_coolskin_sim_mod
  use kinds
  use ufo_constants_mod 
  implicit none
  public :: ufo_coolskin_sim, ufo_coolskin_jac
  private

contains
  subroutine ufo_coolskin_sim(Ts,dTc,S_ns,H_I,H_s,R_nl,Td,u0)
  use kinds 
  implicit none
 
  real(kind=kind_real),        intent(inout) :: Ts,dTc,u0      !dTc to return
  real(kind=kind_real),        intent(in)    :: S_ns,H_I,H_s,R_nl,Td
  
  ! local variables
  integer:: N_i,i
  real(kind=kind_real) :: delta,fc, u, lamda,Q0,Qb
  

  u     = max(0.0002, u0) !friction velocity over water

  dTc   = 0.0
  
  do i = 1,N_i
      Q0     = H_I + H_s + (eps * sig * (Td - dTc)**4 - R_nl)
      Qb     = Q0 + (S_B*cw/(alpha*L_e))*H_I
      lamda = 6.0*(1.0+((alpha*gr*Qb/(Rou*cw))*(16.0*Rou**2.0 * cw**2.0 * v_w**3.0 /\
      k_t**2.0)* (1/u**4.0))**(3.0/4.0))**(-1.0/3.0)
      delta  =lamda * v_w / u
      fc  = 0.0685 + 11.0 * delta - 3.3E-5 /delta * (1.0-exp(-delta/(8.0E-4)))
      dTc  = (H_I + H_s + (eps * sig * (Td - dTc)**4 - R_nl) - S_ns * fc) * delta/k_t
      dTc = max(dTc,0.0)
  enddo

  Ts = Td - dTc
  
  end subroutine ufo_coolskin_sim
  
  !-----------------------------------------------------------------------
  
  subroutine ufo_coolskin_jac(jac,S_ns,H_I,H_s,R_nl,Td,u0)
  use kinds
  implicit none
  
  real(kind=kind_real),        intent(in)    :: S_ns,H_I,H_s,R_nl,Td,u0
  real(kind=kind_real),        intent(out)   :: jac(6) ! jac calculated for all inputs

  real(kind=kind_real) :: delta ,fc ,u
  real(kind=kind_real) :: lamda ,Q0 ,Qb ,Ts ,dTc ,c0 ,y ,Q
  real(kind=kind_real) :: const ,d_lamda_dQb ,dQb_dTs ,d_lamda_dTs,d_delta_dTs,dfc_d_delta,dQ_dTs

  u     = max(0.0002, u0) !friction velocity over water
  
  call ufo_coolskin_sim(Ts,dTc,S_ns,H_I,H_s,R_nl,Td,u)
  
  Q0      = H_I + H_s + (eps * sig * Ts**4 - R_nl)
   
  !constant apears in net heat equation
  c0      = S_B*cw/(alpha*L_e)   
  Qb      = Q0 + c0*H_I

  !Saunder’s constant
  lamda   = 6.0*(1.0+((alpha*gr*Qb/(Rou*cw))*(16.0*Rou**2.0 * cw**2.0 * v_w**3.0 /\
  k_t**2.0)* (1/u**4.0))**(3.0/4.0))**(-1.0/3.0)
  
  ! cool layer thickness
  delta   = lamda * v_w / u
  
  !solar absorbption profile in cool layer
  fc      = 0.0685 + 11.0 * delta - 3.3E-5 /delta * (1.0-exp(-delta/(8.0E-4)))
  
  ! net heat in cool layer
  Q       = H_I + H_s + (eps * sig * (Td - dTc)**4 - R_nl) - S_ns * fc
  const   = ((alpha*gr/(Rou*cw))*(16.0*Rou**2.0 * cw**2.0 * v_w**3.0 / k_t**2.0))**(3.0/4.0)
  
  ! calculate d(fc)/d(delta)
  dfc_d_delta  = 11 + 3.3E-5 / delta**2 *(1.0-exp(-delta/(8.0E-4))) - 3.3E-5\
  / delta * (1.0/(8.0E-4) * exp(-delta/(8.0E-4)))
  
  ! calculate d(lamda)/d(Qb)
  d_lamda_dQb  = -6.0/3.0 *(1.0+const *(Qb/u**4) **(3.0/4.0))**(-4.0/3.0) *\
  (3.0/4.0) * const * (Qb/u**4) ** (-1.0/4.0)
  
  ! this apears in several of jacobians, I decided to calculate it once (4 eps * sig * Ts**3)
  y            = 4 * eps * sig * Ts**3

  ! d(Ts)/d(S_ns)
  jac(1) = fc * delta /(k_t+y*(delta+(Q - S_ns*delta *dfc_d_delta)*v_w/u*d_lamda_dQb))
  
  ! d(Ts)/d(H_I)
  jac(2) = -((1+c0)*(delta+Q*v_w/u*d_lamda_dQb-S_ns*dfc_d_delta*v_w/u*d_lamda_dQb))/\
  (y*(delta+Q*v_w/u*d_lamda_dQb-S_ns*dfc_d_delta*v_w/u*d_lamda_dQb)+k_t)
  
  ! d(Ts)/d(H_s)
  jac(3) = -((delta+Q*v_w/u*d_lamda_dQb-S_ns*dfc_d_delta*v_w/u*d_lamda_dQb))/\
  (y*(delta+Q*v_w/u*d_lamda_dQb -S_ns*dfc_d_delta*v_w/u*d_lamda_dQb)+k_t)
  
  ! d(Ts)/d(R_nl)
  jac(4) = (delta +(Q-S_ns*dfc_d_delta*delta)*(v_w/u*d_lamda_dQb))/\
  (k_t+(Q-S_ns*dfc_d_delta*delta)*(y*v_w/u*d_lamda_dQb))
  
  ! d(Ts)/d(Td)
  jac(5) = k_t/(k_t + y *(delta - delta*S_ns*dfc_d_delta +Q) * v_w/u * d_lamda_dQb)
  	
  ! d(Ts)/d(u)
  jac(6) = - (Q-S_ns* delta*dfc_d_delta)*(v_w/u**2.0 *lamda + (4.0*u**(-6.0))*d_lamda_dQb*Qb)/\
  (k_t+y *(delta +(Q-delta*S_ns*dfc_d_delta)*(v_w*d_lamda_dQb*u**(-5.0))))
  
  end subroutine ufo_coolskin_jac
  
  
end module ufo_coolskin_sim_mod
