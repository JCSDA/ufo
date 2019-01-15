! (C) Copyright 2017-2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!>Stubbed  Fortran module to handle gnssro bending angle observations following 
!> the ROPP (2018 Aug) implementation

module ufo_gnssro_bndropp1d_mod
  use iso_c_binding
  use ufo_vars_mod
  use ufo_geovals_mod
  use ufo_geovals_mod_c, only: ufo_geovals_registry
  use ufo_basis_mod,     only: ufo_basis
  use vert_interp_mod
  use lag_interp_mod, only: lag_interp_const, lag_interp_smthWeights
  use obsspace_mod  

  use kinds
  implicit none
  public             :: ufo_gnssro_bndropp1d
  private

  !> Fortran derived type for gnssro trajectory
  type, extends(ufo_basis) :: ufo_gnssro_BndROPP1D
  contains
    procedure :: simobs    => ufo_gnssro_bndropp1d_simobs
  end type ufo_gnssro_BndROPP1D

  contains
! ------------------------------------------------------------------------------
  subroutine ufo_gnssro_bndropp1d_simobs(self, geovals, hofx, obss)
   !! use ropp_fm_types, only: State1dFM
   !! use ropp_fm_types, only: Obs1dBangle
   !! use datetimetypes, only: dp
      implicit none
      class(ufo_gnssro_BndROPP1D), intent(in):: self
      type(ufo_geovals), intent(in)          :: geovals
      real(kind_real),   intent(inout)       :: hofx(:)
      type(c_ptr), value, intent(in)         :: obss

   !!   type(State1dFM)                 :: x
   !!   type(Obs1dBangle)               :: y

      character(len=*), parameter     :: myname_="ufo_gnssro_bndropp1d_simobs_stub"
   !!   real(kind=dp)                   :: ob_time

      integer, parameter              :: max_string = 800

      character(max_string)           :: err_msg
      character(len=250)              :: record
      integer                         :: iobs
      integer                         :: nlev, nobs
      integer                         :: nvprof
      integer, allocatable, dimension(:)      :: ichk
      type(ufo_geoval), pointer       :: t, q, prs, z, z_sfc
      real(kind_real), allocatable    :: obsLat(:), obsLon(:), obsImpP(:), obsLocR(:), obsGeoid(:)
      real(kind_real), allocatable    :: obsYYYY(:), obsMM(:), obsDD(:), obsHH(:), obsMN(:), obsSS(:)
      integer :: obss_nobs

      write(*,*) "TRACE: ufo_gnssro_bndropp1d_simobs: begin"
      ! check if nobs is consistent in geovals & hofx
      if (geovals%nobs /= size(hofx)) then
        write(err_msg,*) myname_, ' error: nobs inconsistent!'
        call abor1_ftn(err_msg)
      endif
      ! check if prs (pressure at model levels)  variable is in geovals and get it
      call ufo_geovals_get_var(geovals, var_prs, prs)
      ! check if specific humidity variable is in geovals and get it
      call ufo_geovals_get_var(geovals, var_q, q)
      ! check if geopotential height variable is in geovals and get it
      call ufo_geovals_get_var(geovals, var_z, z)
      ! check if surface geopotential height variable is in geovals and get it
      call ufo_geovals_get_var(geovals, var_sfc_z, z_sfc)
      ! check if t variable is in geovals and get it
      call ufo_geovals_get_var(geovals, var_t, t)

      nlev  = q%nval ! number of model levels
      nobs  = geovals%nobs ! number of observations

      ! read observation vectors
      obss_nobs = obsspace_get_nobs(obss)
      allocate(obsLon(obss_nobs))
      allocate(obsLat(obss_nobs))
      allocate(obsImpP(obss_nobs))
      allocate(obsLocR(obss_nobs))
      allocate(obsGeoid(obss_nobs))

      call obsspace_get_db(obss, "Metadata", "longitude", obsLon)
      call obsspace_get_db(obss, "Metadata", "latitude", obsLat) 
      call obsspace_get_db(obss, "Metadata", "impact_parameter", obsImpP)
      call obsspace_get_db(obss, "Metadata", "earth_radius_of_curvature", obsLocR)
      call obsspace_get_db(obss, "Metadata", "geoid_height_above_reference_ellipsoid", obsGeoid)

      nvprof = 1 ! number of vertical profiles (occultation points)
      allocate(ichk(nvprof))
      ichk(:) = 0   ! this will hold QC values for observation from QC flags

      write(record,*) "DEBUG: ufo_gnssro_bndropp1d_simobs: begin observation loop  ", nobs
      obs_loop: do iobs = 1, nobs 

      !!  hofx(iobs) = y%bangle(nvprof)  ! nvprof is just one point
         hofx(iobs)  = 0.0  
      end do obs_loop
      
      deallocate(obsLat) !Note: to be removed
      deallocate(obsLon)
      deallocate(obsImpP)
      deallocate(obsLocR)
      deallocate(obsGeoid)
      write(*,*) "TRACE: ufo_gnssro_bndropp1d_simobs: completed"

   end subroutine ufo_gnssro_bndropp1d_simobs
! ------------------------------------------------------------------------------

end module ufo_gnssro_bndropp1d_mod
