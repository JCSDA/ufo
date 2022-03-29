module satcolumn_mod
use kinds
use ufo_constants_mod, only: zero, one, gas_constant, avogadro, M_dryair, grav
use vert_interp_mod, only: vert_interp_weights
contains

subroutine simulate_column_ob(nlayers_obs, nlayers_model, avgkernel_obs, &
                              prsi_obs, prsi_model, &
                              profile_model, hofx, &
                              troplev_obs, airmass_tot, airmass_trop)
  implicit none
  integer, intent(in   ) :: nlayers_obs, nlayers_model
  real(kind_real), intent(in   ), dimension(nlayers_obs) :: avgkernel_obs
  real(kind_real), intent(in   ), dimension(nlayers_obs+1) :: prsi_obs
  real(kind_real), intent(in   ), dimension(nlayers_model+1) :: prsi_model
  real(kind_real), intent(in   ), dimension(nlayers_model) :: profile_model
  real(kind_real), intent(  out) :: hofx
  real(kind_real), intent(in   ), optional :: airmass_tot, airmass_trop
  integer, intent(in   ), optional :: troplev_obs
  real(kind_real) :: airmass_ratio, wf_a, wf_b
  real(kind_real), dimension(nlayers_obs) :: avgkernel_use, profile_obslayers
  real(kind_real), dimension(nlayers_obs+1) :: pobs
  real(kind_real), dimension(nlayers_model+1) :: pmod
  integer, parameter :: max_string=800
  character(len=max_string) :: err_msg
  integer :: k, j, wi_a, wi_b
  logical :: troposphere
  troposphere = .false.
  ! determine if we are computing a total column or just tropospheric column
  if (present(airmass_tot) .and. present(airmass_trop) .and. present(troplev_obs)) then
     troposphere = .true.
  end if

  if (troposphere) then
    ! compute air mass ratio and apply it to the averaging kernel
    airmass_ratio = airmass_tot / airmass_trop
    avgkernel_use = avgkernel_obs * airmass_ratio

    ! set all layers to zero above tropopause layer
    avgkernel_use(1:troplev_obs) = zero
  else
    avgkernel_use = avgkernel_obs
  end if

  ! make sure we use the highest surface pressure 
  !(do not miss surface values or interpolate below surface)
  pobs = prsi_obs
  pmod = prsi_model
  if (pmod(nlayers_model+1) > pobs(nlayers_obs+1)) then
     pobs(nlayers_obs+1) = pmod(nlayers_model+1)
  else
     pmod(nlayers_model+1) = pobs(nlayers_obs+1)
  end if

  hofx = zero
  profile_obslayers = zero
  do k=1,nlayers_obs  
     ! get obs layer bound model indexes and weights for staggered 
     ! obs and geoval levels
     call vert_interp_weights(nlayers_model+1, pobs(k), pmod, wi_a, wf_a) 
     call vert_interp_weights(nlayers_model+1, pobs(k+1), pmod, wi_b, wf_b)
     ! when multiple geovals levels are in a obs layer
     if ( wi_a < wi_b ) then
        profile_obslayers(k) = profile_obslayers(k) + profile_model(wi_a) * &
             (pmod(wi_a+1)-pmod(wi_a)) * wf_a / (M_dryair*grav)  
        do j=wi_a+1,wi_b-1
           profile_obslayers(k) = profile_obslayers(k) + profile_model(j) * &
                (pmod(j+1)-pmod(j)) / (M_dryair*grav)
        enddo
        profile_obslayers(k) = profile_obslayers(k) + profile_model(wi_b) * &
             (pmod(wi_b+1)-pmod(wi_b)) * (one-wf_b) / (M_dryair*grav)

     ! when multiple obs layers are in a geovals level
     else if ( wi_a == wi_b ) then
        profile_obslayers(k) = profile_obslayers(k) + profile_model(wi_a) * &
             (pmod(wi_a+1)-pmod(wi_a)) * (wf_a-wf_b) / (M_dryair*grav)

     ! if pressures coordinates are inverted return exception
     else if ( wi_a > wi_b) then
        write(err_msg, *) "Error: inverted pressure coordinate in obs, &
                &convention: top->bottom, decreasing pressures"
        call abor1_ftn(err_msg)
  end if
     ! compute A.x
     hofx = hofx + (avgkernel_use(k) * profile_obslayers(k))
  end do
end subroutine simulate_column_ob

subroutine simulate_column_ob_tl(nlayers_obs, nlayers_model, avgkernel_obs, &
                              prsi_obs, prsi_model, &
                              profile_model, hofx, &
                              troplev_obs, airmass_tot, airmass_trop)
  implicit none
  integer, intent(in   ) :: nlayers_obs, nlayers_model
  real(kind_real), intent(in   ), dimension(nlayers_obs) :: avgkernel_obs
  real(kind_real), intent(in   ), dimension(nlayers_obs+1) :: prsi_obs
  real(kind_real), intent(in   ), dimension(nlayers_model+1) :: prsi_model
  real(kind_real), intent(in   ), dimension(nlayers_model) :: profile_model
  real(kind_real), intent(  out) :: hofx
  real(kind_real), intent(in   ), optional :: airmass_tot, airmass_trop
  integer, intent(in   ), optional :: troplev_obs
  real(kind_real) :: airmass_ratio, wf_a, wf_b
  real(kind_real), dimension(nlayers_obs) :: avgkernel_use, profile_obslayers
  integer, parameter :: max_string=800
  real(kind_real), dimension(nlayers_obs+1) :: pobs
  real(kind_real), dimension(nlayers_model+1) :: pmod
  character(len=max_string) :: err_msg
  integer :: k, j, wi_a, wi_b
  logical :: troposphere

  troposphere = .false.
  ! determine if we are computing a total column or just tropospheric column
  if (present(airmass_tot) .and. present(airmass_trop) .and. present(troplev_obs)) then
     troposphere = .true.
  end if

  if (troposphere) then
    ! compute air mass ratio and apply it to the averaging kernel
    airmass_ratio = airmass_tot / airmass_trop
    avgkernel_use = avgkernel_obs * airmass_ratio

    ! set all layers to zero above tropopause layer
    avgkernel_use(troplev_obs+1:nlayers_obs) = zero
  else
    avgkernel_use = avgkernel_obs
  end if

  ! make sure we use the highest surface pressure
  !(do not miss surface values or interpolate below surface)
  pobs = prsi_obs
  pmod = prsi_model
  if (pmod(nlayers_model+1) > pobs(nlayers_obs+1)) then
     pobs(nlayers_obs+1) = pmod(nlayers_model+1)
  else
     pmod(nlayers_model+1) = pobs(nlayers_obs+1)
  end if

  hofx = zero
  profile_obslayers = zero
  do k=1,nlayers_obs
     ! get obs layer bound model indexes and weights for staggered
     ! obs and geoval levels
     call vert_interp_weights(nlayers_model+1, pobs(k), pmod, wi_a, wf_a)
     call vert_interp_weights(nlayers_model+1, pobs(k+1), pmod, wi_b, wf_b)

     ! when multiple geovals levels are in a obs layer
     if ( wi_a < wi_b ) then
        profile_obslayers(k) = profile_obslayers(k) + profile_model(wi_a) * &
             (pmod(wi_a+1)-pmod(wi_a)) * wf_a / (M_dryair*grav)
        do j=wi_a+1,wi_b-1
           profile_obslayers(k) = profile_obslayers(k) + profile_model(j) * &
                (pmod(j+1)-pmod(j)) / (M_dryair*grav)
        enddo
        profile_obslayers(k) = profile_obslayers(k) + profile_model(wi_b) * &
             (pmod(wi_b+1)-pmod(wi_b)) * (one-wf_b) / (M_dryair*grav)

     ! when multiple obs layers are in a geovals level
     else if ( wi_a == wi_b ) then
        profile_obslayers(k) = profile_obslayers(k) + profile_model(wi_a) * &
             (pmod(wi_a+1)-pmod(wi_a)) * (wf_a-wf_b) / (M_dryair*grav)

     ! if pressures coordinates are inverted return exception
     else if ( wi_a > wi_b) then
        write(err_msg, *) "Error: inverted pressure coordinate in obs &
                &convention: top->bottom, decreasing pressures"
        call abor1_ftn(err_msg)
     end if
     ! compute A.x
     hofx = hofx + (avgkernel_use(k) * profile_obslayers(k))
  end do
end subroutine simulate_column_ob_tl

subroutine simulate_column_ob_ad(nlayers_obs, nlayers_model, avgkernel_obs, &
                              prsi_obs, prsi_model, &
                              profile_model_ad, hofx_ad, &
                              troplev_obs, airmass_tot, airmass_trop)
  implicit none
  integer, intent(in   ) :: nlayers_obs, nlayers_model
  real(kind_real), intent(in   ), dimension(nlayers_obs) :: avgkernel_obs
  real(kind_real), intent(in   ), dimension(nlayers_obs+1) :: prsi_obs
  real(kind_real), intent(in   ), dimension(nlayers_model+1) :: prsi_model
  real(kind_real), intent(inout), dimension(nlayers_model) :: profile_model_ad
  real(kind_real), intent(in   ) :: hofx_ad
  real(kind_real), intent(in   ), optional :: airmass_tot, airmass_trop
  integer, intent(in   ), optional :: troplev_obs
  real(kind_real) :: airmass_ratio, wf_a, wf_b
  real(kind_real), dimension(nlayers_obs) :: avgkernel_use, profile_obslayers_ad
  real(kind_real), dimension(nlayers_obs+1) :: pobs
  real(kind_real), dimension(nlayers_model+1) :: pmod
  integer, parameter :: max_string=800
  character(len=max_string) :: err_msg
  integer :: k, j, wi_a, wi_b
  logical :: troposphere

  troposphere = .false.
  ! determine if we are computing a total column or just tropospheric column
  if (present(airmass_tot) .and. present(airmass_trop) .and. present(troplev_obs)) then
     troposphere = .true.
  end if

  if (troposphere) then
    ! compute air mass ratio and apply it to the averaging kernel
    airmass_ratio = airmass_tot / airmass_trop
    avgkernel_use = avgkernel_obs * airmass_ratio

    ! set all layers to zero above tropopause layer
    avgkernel_use(troplev_obs+1:nlayers_obs) = zero
  else
    avgkernel_use = avgkernel_obs
  end if

  ! make sure we use the highest surface pressure
  !(do not miss surface values or interpolate below surface)
  pobs = prsi_obs
  pmod = prsi_model
  if (pmod(nlayers_model+1) > pobs(nlayers_obs+1)) then
     pobs(nlayers_obs+1) = pmod(nlayers_model+1)
  else
     pmod(nlayers_model+1) = pobs(nlayers_obs+1)
  end if

  profile_obslayers_ad = zero
  do k=1,nlayers_obs
     !A.x ad
     profile_obslayers_ad(k) = profile_obslayers_ad(k) + avgkernel_use(k) * hofx_ad

     ! get obs layer bound model indexes and weights for staggered
     ! obs and geoval levels
     call vert_interp_weights(nlayers_model+1, pobs(k), pmod, wi_a, wf_a)
     call vert_interp_weights(nlayers_model+1, pobs(k+1), pmod, wi_b, wf_b)

     ! when multiple geovals levels are in a obs layer
     if ( wi_a < wi_b ) then
           profile_model_ad(wi_a) = profile_model_ad(wi_a) + profile_obslayers_ad(k) * &
             (pmod(wi_a+1)-pmod(wi_a)) * (wf_a) / (M_dryair*grav)
        do j=wi_a+1,wi_b-1
           profile_model_ad(j) = profile_model_ad(j) + profile_obslayers_ad(k) * &
             (pmod(j+1)-pmod(j)) / (M_dryair*grav)
        enddo
        profile_model_ad(wi_b) = profile_model_ad(wi_b) + profile_obslayers_ad(k) * &
             (pmod(wi_b+1)-pmod(wi_b)) * (one-wf_b) / (M_dryair*grav)

     ! when multiple obs layers are in a geovals level
     else if ( wi_a == wi_b ) then
        profile_model_ad(wi_a) = profile_model_ad(wi_a) + profile_obslayers_ad(k) * &
             (pmod(wi_a+1)-pmod(wi_a)) * (wf_a-wf_b) / (M_dryair*grav)

     ! if pressures coordinates are inverted return exception
     else if ( wi_a > wi_b ) then
        write(err_msg, *) "Error: inverted pressure coordinate in obs &
                &convention: top->bottom, decreasing pressures"
        call abor1_ftn(err_msg)
     end if
  end do
end subroutine simulate_column_ob_ad

end module satcolumn_mod
