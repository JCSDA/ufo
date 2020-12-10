module satcolumn_mod
contains

subroutine simulate_column_ob(nlayers_obs, nlayers_model, avgkernel_obs, &
                              prsl_obs, prsl_model, temp_model, zi_model, &
                              profile_model, hofx, &
                              troplev_obs, airmass_tot, airmass_trop)
  use kinds
  use ufo_constants_mod, only: zero, gas_constant, avogadro
  use vert_interp_mod, only: vert_interp_weights, vert_interp_apply
  implicit none
  integer, intent(in   ) :: nlayers_obs, nlayers_model
  real(kind_real), intent(in   ), dimension(nlayers_obs) :: avgkernel_obs
  real(kind_real), intent(in   ), dimension(nlayers_obs) :: prsl_obs
  real(kind_real), intent(in   ), dimension(nlayers_model) :: prsl_model
  real(kind_real), intent(in   ), dimension(nlayers_model) :: profile_model
  real(kind_real), intent(in   ), dimension(nlayers_model) :: temp_model
  real(kind_real), intent(in   ), dimension(nlayers_model+1) :: zi_model
  real(kind_real), intent(  out) :: hofx
  real(kind_real), intent(in   ), optional :: airmass_tot, airmass_trop
  integer, intent(in   ), optional :: troplev_obs
  real(kind_real) :: airmass_ratio, lnp_ob, wf, dz
  real(kind_real), dimension(nlayers_obs) :: avgkernel_use
  real(kind_real), dimension(nlayers_model) :: lnp_model, profile_model_use
  real(kind_real), dimension(nlayers_obs) :: profile_obslayers
  integer :: k, wi
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

  ! need to convert from kg/kg to molec/m3 for each layer's T and P
  profile_model_use = profile_model * ((avogadro*prsl_model)/(temp_model*gas_constant))

  ! need to convert from molec/m3 to molec/cm2 to sum later
  do k=1,nlayers_model
    dz = (zi_model(k+1) - zi_model(k))
    profile_model_use(k) = profile_model_use(k) * dz ! convert to molec/m2
    profile_model_use(k) = profile_model_use(k) / 1.0e4_kind_real ! to molec/cm2
  end do

  ! need to interpolate model profile layers to observation layers
  lnp_model = log(prsl_model)
  do k=1,nlayers_obs
    lnp_ob = log(prsl_obs(k))
    call vert_interp_weights(nlayers_model, lnp_ob, lnp_model, wi, wf)
    call vert_interp_apply(nlayers_model, profile_model_use, &
                            profile_obslayers(k), wi, wf)
  end do

  ! compute hofx as column integral
  hofx = zero
  do k=1,nlayers_obs
    hofx = hofx + (avgkernel_use(k) * profile_obslayers(k))
  end do

end subroutine simulate_column_ob

subroutine simulate_column_ob_tl(nlayers_obs, nlayers_model, avgkernel_obs, &
                              prsl_obs, prsl_model, temp_model, zi_model, &
                              profile_model, hofx, &
                              troplev_obs, airmass_tot, airmass_trop)
  use kinds
  use ufo_constants_mod, only: zero, gas_constant, avogadro
  use vert_interp_mod, only: vert_interp_weights, vert_interp_apply_tl
  implicit none
  integer, intent(in   ) :: nlayers_obs, nlayers_model
  real(kind_real), intent(in   ), dimension(nlayers_obs) :: avgkernel_obs
  real(kind_real), intent(in   ), dimension(nlayers_obs) :: prsl_obs
  real(kind_real), intent(in   ), dimension(nlayers_model) :: prsl_model
  real(kind_real), intent(in   ), dimension(nlayers_model) :: profile_model
  real(kind_real), intent(in   ), dimension(nlayers_model) :: temp_model
  real(kind_real), intent(in   ), dimension(nlayers_model+1) :: zi_model
  real(kind_real), intent(  out) :: hofx
  real(kind_real), intent(in   ), optional :: airmass_tot, airmass_trop
  integer, intent(in   ), optional :: troplev_obs
  real(kind_real) :: airmass_ratio, lnp_ob, wf, dz
  real(kind_real), dimension(nlayers_obs) :: avgkernel_use
  real(kind_real), dimension(nlayers_model) :: lnp_model, profile_model_use
  real(kind_real), dimension(nlayers_obs) :: profile_obslayers
  integer :: k, wi
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

  ! need to convert from kg/kg to molec/m3 for each layer's T and P
  profile_model_use = profile_model * ((avogadro*prsl_model)/(temp_model*gas_constant))

  ! need to convert from molec/m3 to molec/cm2 to sum later
  do k=1,nlayers_model
    dz = (zi_model(k+1) - zi_model(k))
    profile_model_use(k) = profile_model_use(k) * dz ! convert to molec/m2
    profile_model_use(k) = profile_model_use(k) / 1.0e4_kind_real ! to molec/cm2
  end do

  ! need to interpolate model profile layers to observation layers
  lnp_model = log(prsl_model)
  do k=1,nlayers_obs
    lnp_ob = log(prsl_obs(k))
    call vert_interp_weights(nlayers_model, lnp_ob, lnp_model, wi, wf)
    call vert_interp_apply_tl(nlayers_model, profile_model_use, &
                            profile_obslayers(k), wi, wf)
  end do

  ! compute hofx as column integral
  hofx = zero
  do k=1,nlayers_obs
    hofx = hofx + (avgkernel_use(k) * profile_obslayers(k))
  end do

  end subroutine simulate_column_ob_tl

subroutine simulate_column_ob_ad(nlayers_obs, nlayers_model, avgkernel_obs, &
                              prsl_obs, prsl_model, temp_model, zi_model, &
                              profile_model_ad, hofx_ad, &
                              troplev_obs, airmass_tot, airmass_trop)
  use kinds
  use ufo_constants_mod, only: zero, gas_constant, avogadro
  use vert_interp_mod, only: vert_interp_weights, vert_interp_apply_ad, vert_interp_apply_tl
  implicit none
  integer, intent(in   ) :: nlayers_obs, nlayers_model
  real(kind_real), intent(in   ), dimension(nlayers_obs) :: avgkernel_obs
  real(kind_real), intent(in   ), dimension(nlayers_obs) :: prsl_obs
  real(kind_real), intent(in   ), dimension(nlayers_model) :: prsl_model
  real(kind_real), intent(inout), dimension(nlayers_model) :: profile_model_ad
  real(kind_real), intent(in   ), dimension(nlayers_model) :: temp_model
  real(kind_real), intent(in   ), dimension(nlayers_model+1) :: zi_model
  real(kind_real), intent(in   ) :: hofx_ad
  real(kind_real), intent(in   ), optional :: airmass_tot, airmass_trop
  integer, intent(in   ), optional :: troplev_obs
  real(kind_real) :: airmass_ratio, lnp_ob, wf, dz
  real(kind_real), dimension(nlayers_obs) :: avgkernel_use
  real(kind_real), dimension(nlayers_model) :: lnp_model, profile_model_use_ad
  real(kind_real), dimension(nlayers_obs) :: profile_obslayers_ad
  integer :: k, wi
  logical :: troposphere
  integer :: nlevs

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

  profile_obslayers_ad = zero
  do k=nlayers_obs,1,-1
    profile_obslayers_ad(k) = profile_obslayers_ad(k) + avgkernel_use(k) * hofx_ad
  end do

  profile_model_use_ad = zero
  ! need to interpolate model profile layers to observation layers
  lnp_model = log(prsl_model)
  do k=nlayers_obs,1,-1
    lnp_ob = log(prsl_obs(k))
    call vert_interp_weights(nlayers_model, lnp_ob, lnp_model, wi, wf)
    call vert_interp_apply_ad(nlayers_model, profile_model_use_ad, &
                            profile_obslayers_ad(k), wi, wf)
  end do

  ! adjoint of unit conversions
  do k=nlayers_model,1,-1
    dz = (zi_model(k+1) - zi_model(k))
    profile_model_use_ad(k) = profile_model_use_ad(k) / 1.0e4_kind_real ! to molec/m2
    profile_model_use_ad(k) = profile_model_use_ad(k) * dz ! convert to molec/m3
  end do

  profile_model_ad = profile_model_use_ad * ((avogadro*prsl_model)/(temp_model*gas_constant))

  end subroutine simulate_column_ob_ad

end module satcolumn_mod
