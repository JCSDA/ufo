! (C) Copyright 2024 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module that wraps the GSI field of view code for easier use by JEDI/UFO code.
!> The main goal of this wrapper is to hide certain implementation details of GSI from
!> the calling JEDI/UFO code.

module ufo_fov_mod

  use kinds, only: kind_real

  implicit none
  private

  !> Fortran type that wraps implementation details of the GSI field of view
  type, public :: ufo_fov
    private

    integer :: instr  ! index that labels instruments
    integer :: ichan  ! a surface-sensitive channel, used for computing field of view
    logical :: is_crosstrack  ! is instrument scanning crosstrack (true) or conical (false)
    real(kind_real) :: expansion  ! wavelength-dependent field of view scaling

  contains
    procedure :: setup  => ufo_fov_setup
    procedure :: delete => ufo_fov_delete
    procedure :: fov_ellipse => ufo_fov_ellipse
    procedure :: antenna_power_within_fov => ufo_fov_inside_fov
  end type ufo_fov

contains

  ! ------------------------------------------------------------------------------
  !> Helper to convert from a JEDI sensor name to the internal GSI quantities held in the
  !> `ufo_fov` type.
  !>
  !> To add support for new sensors, look for GSI's values in the read_{airs,atms,etc} files.
  !> GSI sets values for these parameter while reading in the data. The mapping from sensor name
  !> to the `instr` integer can also be found in comments in the files calc_fov_{crosstrk,conical}.
  subroutine set_gsi_interface_quantities(self, sensor, platform)
    class(ufo_fov), intent(inout) :: self
    character(len=*), intent(in) :: sensor
    character(len=*), intent(in) :: platform  ! note GSI calls this "satellite"

    select case(trim(sensor))
      case('amsua')
        self%instr = 11
        self%ichan = 15
        self%is_crosstrack = .true.
        self%expansion = 2.9_kind_real
      case('atms')
        self%instr = 20
        self%ichan = 16
        self%is_crosstrack = .true.
        self%expansion = 2.9_kind_real
      case('cris')
        self%instr = 17
        self%ichan = -999
        self%is_crosstrack = .true.
        self%expansion = 1.0_kind_real
      case('cris-fsr')  ! same as 'cris' above
        self%instr = 17
        self%ichan = -999
        self%is_crosstrack = .true.
        self%expansion = 1.0_kind_real
      case('iasi')
        self%instr = 18
        self%ichan = -999
        self%is_crosstrack = .true.
        self%expansion = 1.0_kind_real
      case('ssmis')
        if (trim(platform) == 'f17') then
          self%instr = 27
        else
          call abor1_ftn("Error in ufo_fov_mod: unknown satellite/platform " &
                         // trim(platform) // " for sensor ssmis")
        end if
        self%ichan = 1
        self%is_crosstrack = .false.
        self%expansion = 2.9_kind_real
      case default
        call abor1_ftn("Error in ufo_fov_mod: unknown input sensor name " // trim(sensor))
    end select

  end subroutine set_gsi_interface_quantities

  ! ------------------------------------------------------------------------------
  subroutine ufo_fov_setup(self, sensor, platform, valid, npoly_out)
    use calc_fov_crosstrk, only: npoly, instrument_init_crosstrk => instrument_init
    use calc_fov_conical, only: instrument_init_conical => instrument_init

    class(ufo_fov), intent(inout) :: self
    character(len=*), intent(in) :: sensor
    character(len=*), intent(in) :: platform  ! note GSI calls this satellite
    logical, intent(out) :: valid
    integer, intent(out) :: npoly_out

    npoly_out = npoly

    call set_gsi_interface_quantities(self, sensor, platform)

    if (self%is_crosstrack) then
      call instrument_init_crosstrk(self%instr, platform, self%expansion, valid)
    else
      call instrument_init_conical(self%instr, platform, self%expansion, valid)
    end if

  end subroutine ufo_fov_setup

  ! ------------------------------------------------------------------------------
  subroutine ufo_fov_delete(self)
    use calc_fov_crosstrk, only: fov_cleanup

    class(ufo_fov), intent(inout) :: self

    if (self%is_crosstrack) then
      call fov_cleanup
    end if

  end subroutine ufo_fov_delete

  ! ------------------------------------------------------------------------------
  !> Compute the field of view ellipse, approximated as a polygon. The underlying GSI code uses a
  !> 29-sided polygon (30 vertices with first = last).
  !>
  !> Note that the GSI algorithm is incorrect for ellipses that overlap either the North or South
  !> poles (probably due to a bad branch choice in some trig function?). This is NOT checked, so
  !> the user should ensure the ellipse is not too close to either pole.
  subroutine ufo_fov_ellipse(self, sensor, scan_position, sat_azimuth_angle, &
                             fov_center_lon, fov_center_lat, fov_ellipse_lons, fov_ellipse_lats)
    use calc_fov_crosstrk, only: npoly, fov_check, fov_ellipse_crosstrk
    use calc_fov_conical, only: fov_ellipse_conical

    class(ufo_fov), intent(in) :: self
    character(len=*), intent(in) :: sensor
    integer, intent(in) :: scan_position
    real(kind_real), intent(in) :: sat_azimuth_angle
    real(kind_real), intent(in) :: fov_center_lon
    real(kind_real), intent(in) :: fov_center_lat
    real(kind_real), intent(out) :: fov_ellipse_lons(npoly)
    real(kind_real), intent(out) :: fov_ellipse_lats(npoly)

    logical :: valid

    if (self%is_crosstrack) then
      call fov_check(scan_position, self%instr, self%ichan, valid)
      if (.not. valid) then
        call abor1_ftn("Error in ufo_fov_mod: invalid field of view")
      end if

      call fov_ellipse_crosstrk(scan_position, self%instr, sat_azimuth_angle, &
                                fov_center_lat, fov_center_lon, fov_ellipse_lats, fov_ellipse_lons)
    else
      call fov_ellipse_conical(self%ichan, sat_azimuth_angle, &
                               fov_center_lat, fov_center_lon, fov_ellipse_lats, fov_ellipse_lons)
    end if

  end subroutine ufo_fov_ellipse

  ! ------------------------------------------------------------------------------
  !> Compute the antenna power for a test point within a field of view. The result is the relative
  !> antenna power in [0-1], or 0 if the test point is outside the field of view ellipse.
  subroutine ufo_fov_inside_fov(self, sensor, scan_position, sat_azimuth_angle, &
                                fov_center_lon, fov_center_lat, test_lon, test_lat, antenna_power)
    use calc_fov_crosstrk, only: inside_fov_crosstrk
    use calc_fov_conical, only: inside_fov_conical

    class(ufo_fov), intent(in) :: self
    character(len=*), intent(in) :: sensor
    integer, intent(in) :: scan_position
    real(kind_real), intent(in) :: sat_azimuth_angle
    real(kind_real), intent(in) :: fov_center_lon
    real(kind_real), intent(in) :: fov_center_lat
    real(kind_real), intent(in) :: test_lon
    real(kind_real), intent(in) :: test_lat
    real(kind_real), intent(out) :: antenna_power

    if (self%is_crosstrack) then
      call inside_fov_crosstrk(self%instr, scan_position, sat_azimuth_angle, &
                               fov_center_lat, fov_center_lon, test_lat, test_lon, &
                               self%expansion, self%ichan, antenna_power)
    else
      call inside_fov_conical(self%instr, self%ichan, sat_azimuth_angle, &
                              fov_center_lat, fov_center_lon, test_lat, test_lon, &
                              self%expansion, antenna_power)
    end if

  end subroutine ufo_fov_inside_fov

  ! ------------------------------------------------------------------------------

end module ufo_fov_mod
