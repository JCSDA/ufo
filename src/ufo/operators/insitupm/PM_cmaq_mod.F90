! (C) Copyright 2021.
!
! This software is developed by NOAA/NWS/EMC under the Apache 2.0 license
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! This module adds up CMAQ aerosol species at observation locations. 
! Scaling factors for three modes (Aitken - At; accumulation - Ac; coarse - Co) can be applied to derive PM2.5 from the total PM
!

    MODULE PM_cmaq_mod

    use aero_kinds_mod
    Implicit None

    PRIVATE
    PUBLIC get_PM_cmaq        
    PUBLIC get_PM_cmaq_tl
    PUBLIC get_PM_cmaq_ad

CONTAINS

!-----------------

 subroutine get_PM_cmaq (km, nv, nobs, nq, modes,       &
                         vf, qm, wi, wf, pm)

! returns surface PM in converted unit (ug/m3) from aerosol Mass Mixing Ratio 

  use vert_interp_mod, only: vert_interp_apply

  implicit NONE

  integer,               intent(in)  :: km               ! number of vertical layers
  integer,               intent(in)  :: nv               ! number of simulated variables
  integer,               intent(in)  :: nobs             ! number of profiles
  integer,               intent(in)  :: nq               ! number of tracers 
  integer,               optional,  intent(in)  :: modes(nq)        ! cmaq modes of tracers 
  integer,               intent(in)  :: wi(nobs)
  real(kind=kind_real),  intent(in)  :: wf(nobs)
  real(kind=kind_real),  intent(in)  :: qm(nq,km,nobs)   ! speciated mass at all levels in ug/m3 
  real(kind=kind_real),  optional,  intent(in)  :: vf(3,km,nobs)     ! At, Ac, Co mode scaling factors
  real(kind=kind_real),  intent(out) :: pm(nv,nobs)      ! PM in ug/m3

!                               ---

  integer :: iq, i, k, iv
  real(kind=kind_real) :: pm_all(km, nobs)      ! PM, all layers, in ug/m3

  pm = 0.0_kind_real

   do iv = 1, nv 

   pm_all = 0.0_kind_real

   if ( present(vf) .and. present(modes)) then
!    Sum over the tracers
   do iq = 1, nq
!    Loop over nobs, km
!    --------------------------
     do i = 1, nobs
           do k =1, km

               pm_all(k,i) = pm_all(k,i) + vf(modes(iq),k,i) * qm(iq,k,i)   
                
           end do  ! end km
     end do  ! end nobs
   end do ! end tracers

   else
!    Sum over the tracers
   do iq = 1, nq
!    Loop over nobs, km
!    --------------------------
     do i = 1, nobs
           do k =1, km

               pm_all(k,i) = pm_all(k,i) + qm(iq,k,i)   
                
           end do  ! end km
     end do  ! end nobs
   end do ! end tracers
  end if

   do i = 1, nobs

    call vert_interp_apply(km, pm_all(:,i), &
                            pm(iv,i), wi(i), wf(i))
   end do  ! end nobs

  end do
 
  end subroutine get_PM_cmaq

!-------------------------------------

  subroutine get_PM_cmaq_tl(km, nv, nobs, nq, modes, vf, &
                            qm_tl, wi, wf, pm_tl)

  use vert_interp_mod, only: vert_interp_apply_tl
  implicit none
  integer, intent(in)    :: km                      ! number of layers
  integer, intent(in)    :: nv                      ! number of simulated variables
  integer, intent(in)    :: nobs                    ! number of profiles
  integer, intent(in)    :: nq                      ! number of tracers 
  integer, optional,  intent(in)    :: modes(nq)               ! cmaq modes of tracers 
  integer, intent(in)    :: wi(nobs)
  real(kind=kind_real),    intent(in)    :: wf(nobs)
  real(kind=kind_real),    optional,  intent(in)    :: vf(3,km,nobs)     ! At, Ac, Co mode scaling factors
  real(kind=kind_real),    intent(in)    :: qm_tl( nq, km, nobs)  
  real(kind=kind_real),    intent(inout) :: pm_tl(nv,nobs)   

  integer :: ob, tr, lv, var  
  real(kind=kind_real) :: pm_all_tl(km, nobs)      ! PM, all layers, in ug/m3

  pm_tl = 0.0_kind_real

  do var = 1, nv

   pm_all_tl = 0.0_kind_real

   if ( present(vf) .and. present(modes)) then
   do ob = 1, nobs
    do lv = 1, km

     do tr = 1,nq
          pm_all_tl(lv, ob) = pm_all_tl(lv, ob) + vf(modes(tr),lv,ob) * qm_tl(tr,lv,ob)
     enddo

    enddo  ! end km
   enddo   ! end nobs
   else
   do ob = 1, nobs
    do lv = 1, km

     do tr = 1,nq
          pm_all_tl(lv, ob) = pm_all_tl(lv, ob) + qm_tl(tr,lv,ob)
     enddo

    enddo  ! end km
   enddo   ! end nobs
   end if

   do ob = 1, nobs

    call vert_interp_apply_tl(km, pm_all_tl(:,ob), &
                               pm_tl(var,ob), wi(ob), wf(ob))
   end do  ! end nobs
  
  end do  
  end subroutine get_PM_cmaq_tl

! -----------------------------------
  subroutine get_PM_cmaq_ad(km, nv, nobs, nq, modes,   &
                              vf, wi, wf, pm_ad, qm_ad)

  use vert_interp_mod, only: vert_interp_apply_ad
  implicit none
  integer, intent(in)    :: km                       ! number of layers
  integer, intent(in)    :: nv                       ! number of simulated variables
  integer, intent(in)    :: nobs                     ! number of profiles
  integer, intent(in)    :: nq                       ! number of tracers
  integer, optional, intent(in)    :: modes(nq)                ! cmaq modes of tracers 
  integer, intent(in)    :: wi(nobs)
  real(kind=kind_real),    intent(in)    :: wf(nobs)
  real(kind=kind_real), optional,    intent(in)    :: vf(3,km,nobs)     ! At, Ac, Co mode scaling factors
  real(kind=kind_real),    intent(in)    :: pm_ad(nv,nobs)     ! PM adjoint
  real(kind=kind_real),    intent(out)   :: qm_ad( nq, km, nobs)      ! aerosol concentration adjoint

  integer :: ob, tr, lv, var
  real(kind=kind_real) :: pm_all_ad(km,nobs)       ! PM adjoint, all layers

  qm_ad = 0.0_kind_real

  do var = 1, nv
   pm_all_ad = 0.0_kind_real

   do ob = nobs,1,-1

    call vert_interp_apply_ad(km, pm_all_ad(:,ob), &
                              pm_ad(var,ob), wi(ob), wf(ob))
   end do  ! end nobs

   if ( present(vf) .and. present(modes)) then
   do ob=nobs,1,-1
    do lv=km,1,-1

      do tr=nq,1,-1
        qm_ad(tr, lv, ob) = qm_ad(tr, lv, ob) + vf(modes(tr), lv, ob) * pm_all_ad(lv, ob)
      end do

     end do
   end do

   else
   do ob=nobs,1,-1
    do lv=km,1,-1

      do tr=nq,1,-1
        qm_ad(tr, lv, ob) = qm_ad(tr, lv, ob) + pm_all_ad(lv, ob)
      end do

    end do
   end do
   end if

  end do

  end subroutine get_PM_cmaq_ad

end module PM_cmaq_mod
