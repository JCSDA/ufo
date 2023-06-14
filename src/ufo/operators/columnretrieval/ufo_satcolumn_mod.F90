module satcolumn_mod
use kinds
use ufo_constants_mod, only: zero, one, gas_constant, avogadro, M_dryair, grav
use vert_interp_mod, only: vert_interp_weights
contains

subroutine stretch_vertices(nzobs, nzmod, pobsin, pobsout, pmodin, pmodout, &
                            stretch)
  implicit none
  integer, intent(in   ) :: nzobs, nzmod
  real(kind_real), intent(in   ), dimension(nzobs+1) :: pobsin
  real(kind_real), intent(in   ), dimension(nzmod+1) :: pmodin
  real(kind_real), intent(inout), dimension(nzobs+1) :: pobsout
  real(kind_real), intent(inout), dimension(nzmod+1) :: pmodout
  character(len=:), intent(in   ), allocatable :: stretch

  ! stretch "bottom" makes sure we use the highest surface pressure
  ! (do not miss surface values or interpolate below surface)
  ! stretch "top" makes sure we use the lowest pressure of top vertice
  pobsout = pobsin
  pmodout = pmodin

  if (trim(stretch) == "bottom" .or. &
      trim(stretch) == "topbottom" ) then
     pmodout(nzmod+1) = max(pmodin(nzmod+1), pobsin(nzobs+1))
     pobsout(nzobs+1) = max(pmodin(nzmod+1), pobsin(nzobs+1))
  end if
  if (trim(stretch) == "top" .or. &
      trim(stretch) == "topbottom" ) then
     pmodout(1) = min(pmodin(1), pobsin(1))
     pobsout(1) = min(pmodin(1), pobsin(1))
  end if
end subroutine stretch_vertices

subroutine simulate_column_ob(nlayers_obs, nlayers_model, avgkernel_obs, &
                              prsi_obs, prsi_model, profile_model, hofx, &
                              stretch)
  implicit none
  integer, intent(in   ) :: nlayers_obs, nlayers_model
  real(kind_real), intent(in   ), dimension(nlayers_obs) :: avgkernel_obs
  real(kind_real), intent(in   ), dimension(nlayers_obs+1) :: prsi_obs
  real(kind_real), intent(in   ), dimension(nlayers_model+1) :: prsi_model
  real(kind_real), intent(in   ), dimension(nlayers_model) :: profile_model
  real(kind_real), intent(  out) :: hofx
  real(kind_real) :: airmass_ratio, wf_a, wf_b
  real(kind_real), dimension(nlayers_obs) :: profile_obslayers
  real(kind_real), dimension(nlayers_obs+1) :: pobs
  real(kind_real), dimension(nlayers_model+1) :: pmod
  integer, parameter :: max_string=800
  character(len=max_string) :: err_msg
  character(len=:), intent(in   ), allocatable :: stretch
  integer :: k, j, wi_a, wi_b

  call stretch_vertices(nlayers_obs, nlayers_model, prsi_obs, pobs, prsi_model, &
                        pmod, stretch)

  hofx = zero
  profile_obslayers = zero
  do k=1,nlayers_obs
     ! get obs layer bound model indexes and weights for staggered
     ! obs and geoval levels
     call vert_interp_weights(nlayers_model+1, pobs(k), pmod, wi_a, wf_a)
     call vert_interp_weights(nlayers_model+1, pobs(k+1), pmod, wi_b, wf_b)

     !check if pmod is monotonic and decreasing
     if ((pmod(wi_a+1) < pmod(wi_a)) .or. (pmod(wi_b+1) < pmod(wi_b))) then
       write(err_msg, *) "Error: inverted pressure coordinate in geovals, &
               &convention: top->bottom, decreasing pressures"
       call abor1_ftn(err_msg)
     end if

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
     hofx = hofx + (avgkernel_obs(k) * profile_obslayers(k))
  end do
end subroutine simulate_column_ob

subroutine simulate_column_ob_tl(nlayers_obs, nlayers_model, avgkernel_obs, &
                              prsi_obs, prsi_model, profile_model, hofx, &
                              stretch)
  implicit none
  integer, intent(in   ) :: nlayers_obs, nlayers_model
  real(kind_real), intent(in   ), dimension(nlayers_obs) :: avgkernel_obs
  real(kind_real), intent(in   ), dimension(nlayers_obs+1) :: prsi_obs
  real(kind_real), intent(in   ), dimension(nlayers_model+1) :: prsi_model
  real(kind_real), intent(in   ), dimension(nlayers_model) :: profile_model
  real(kind_real), intent(  out) :: hofx
  real(kind_real) :: airmass_ratio, wf_a, wf_b
  real(kind_real), dimension(nlayers_obs) :: profile_obslayers
  integer, parameter :: max_string=800
  real(kind_real), dimension(nlayers_obs+1) :: pobs
  real(kind_real), dimension(nlayers_model+1) :: pmod
  character(len=max_string) :: err_msg
  character(len=:), intent(in   ), allocatable :: stretch
  integer :: k, j, wi_a, wi_b

  call stretch_vertices(nlayers_obs, nlayers_model, prsi_obs, pobs, prsi_model, &
                        pmod, stretch)

  hofx = zero
  profile_obslayers = zero
  do k=1,nlayers_obs
     ! get obs layer bound model indexes and weights for staggered
     ! obs and geoval levels
     call vert_interp_weights(nlayers_model+1, pobs(k), pmod, wi_a, wf_a)
     call vert_interp_weights(nlayers_model+1, pobs(k+1), pmod, wi_b, wf_b)

     !check if pmod is monotonic and decreasing
     if ((pmod(wi_a+1) < pmod(wi_a)) .or. (pmod(wi_b+1) < pmod(wi_b))) then
       write(err_msg, *) "Error: inverted pressure coordinate in geovals, &
               &convention: top->bottom, decreasing pressures"
       call abor1_ftn(err_msg)
     end if

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
     hofx = hofx + (avgkernel_obs(k) * profile_obslayers(k))
  end do
end subroutine simulate_column_ob_tl

subroutine simulate_column_ob_ad(nlayers_obs, nlayers_model, avgkernel_obs, &
                              prsi_obs, prsi_model, profile_model_ad, hofx_ad, &
                              stretch)
  implicit none
  integer, intent(in   ) :: nlayers_obs, nlayers_model
  real(kind_real), intent(in   ), dimension(nlayers_obs) :: avgkernel_obs
  real(kind_real), intent(in   ), dimension(nlayers_obs+1) :: prsi_obs
  real(kind_real), intent(in   ), dimension(nlayers_model+1) :: prsi_model
  real(kind_real), intent(inout), dimension(nlayers_model) :: profile_model_ad
  real(kind_real), intent(in   ) :: hofx_ad
  real(kind_real) :: airmass_ratio, wf_a, wf_b
  real(kind_real), dimension(nlayers_obs) :: profile_obslayers_ad
  real(kind_real), dimension(nlayers_obs+1) :: pobs
  real(kind_real), dimension(nlayers_model+1) :: pmod
  integer, parameter :: max_string=800
  character(len=max_string) :: err_msg
  character(len=:), intent(in   ), allocatable :: stretch
  integer :: k, j, wi_a, wi_b

  call stretch_vertices(nlayers_obs, nlayers_model, prsi_obs, pobs, prsi_model, &
                        pmod, stretch)

  profile_obslayers_ad = zero
  do k=1,nlayers_obs
     !A.x ad
     profile_obslayers_ad(k) = profile_obslayers_ad(k) + avgkernel_obs(k) * hofx_ad

     ! get obs layer bound model indexes and weights for staggered
     ! obs and geoval levels
     call vert_interp_weights(nlayers_model+1, pobs(k), pmod, wi_a, wf_a)
     call vert_interp_weights(nlayers_model+1, pobs(k+1), pmod, wi_b, wf_b)

     !check if pmod is monotonic and decreasing
     if ((pmod(wi_a+1) < pmod(wi_a)) .or. (pmod(wi_b+1) < pmod(wi_b))) then
       write(err_msg, *) "Error: inverted pressure coordinate in geovals, &
               &convention: top->bottom, decreasing pressures"
       call abor1_ftn(err_msg)
     end if

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
