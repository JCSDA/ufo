! (C) Copyright 2026 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Split ObsDiagnostics (hofxdiags) variable names into their parts.
!>
!> Shared by the CRTM and RTTOV operators. Built with the utils, so it does not depend on
!> which radiative transfer libraries were found.

module ufo_hofxdiags_parse_mod

use ufo_vars_mod,    only: MAXVARLEN
use ufo_geovals_mod, only: ufo_geovals

implicit none
private

public :: ufo_parse_hofxdiags

!> Separates the dependent and independent variable of a Jacobian name
character(len=*), parameter :: jacobianstr = "_jacobian_"

!> Channel of a name that does not end in a channel number
integer, parameter, public :: no_channel = -9999

contains

! ------------------------------------------------------------------------------
!> Split each hofxdiags variable name into dependent variable ystr, independent
!> variable xstr and channel ch:
!>
!>   <ystr>_jacobian_<xstr>_<chstr>  Jacobian dy/dx of channel ch
!>   <ystr>_jacobian_<xstr>          Jacobian dy/dx, no channel
!>   <ystr>_<chstr>                  forward diagnostic of channel ch
!>   <ystr>                          forward diagnostic, no channel
!>
!> ch is the part after the last "_" when that part is an unsigned integer, otherwise
!> no_channel. xstr is "" for forward diagnostics. Empty names give ystr = xstr = ""
!> and ch = no_channel. jacobian_needed is true if any name contains "_jacobian_".
!>
!> The parser only splits the hofxdiags variable names; it does not check that
!> ystr/xstr are diagnostics or Jacobians the calling operator can compute.
!> Bias correction terms (<predictor>_<ch>, e.g. constant_1), which ObsBias
!> reserves in hofxdiags, parse as forward diagnostics.
!>
!> Aborts if ystr, xstr or ch is shorter than hofxdiags%nvar, if a name starts with
!> "_jacobian_" or has nothing after it, or if a part does not fit in ystr/xstr.
!> caller (e.g. "ufo_radiancecrtm_simobs", "ufo_rttov_parse_hofxdiags") is prepended
!> to the abort messages.
subroutine ufo_parse_hofxdiags(hofxdiags, ystr, xstr, ch, jacobian_needed, caller)
type(ufo_geovals),  intent(in)  :: hofxdiags
character(len=MAXVARLEN), intent(out) :: ystr(:), xstr(:)
integer,            intent(out) :: ch(:)
logical,            intent(out) :: jacobian_needed
character(len=*),   intent(in)  :: caller

character(len=MAXVARLEN) :: name
integer :: jvar, nlen, last_underscore, ijac, xbeg, xend, yend, ierr

if (size(ystr) < hofxdiags%nvar .or. size(xstr) < hofxdiags%nvar .or. &
    size(ch) < hofxdiags%nvar) then
  call abor1_ftn(trim(caller) // ": hofxdiags name arrays are smaller than hofxdiags%nvar")
end if

ystr = ""
xstr = ""
ch = no_channel
jacobian_needed = .false.

do jvar = 1, hofxdiags%nvar
  name = hofxdiags%variables(jvar)
  nlen = len_trim(name)
  if (nlen < 1) cycle

  ! channel: if the text after the last "_" is all digits, read it into ch as an integer.
  ! A failed read (e.g. overflow) leaves ch = no_channel.
  last_underscore = index(name(1:nlen), "_", back=.true.)
  if (last_underscore > 0 .and. last_underscore < nlen) then
    if (verify(name(last_underscore+1:nlen), "0123456789") == 0) then
      read(name(last_underscore+1:nlen), *, iostat=ierr) ch(jvar)
      if (ierr /= 0) ch(jvar) = no_channel
    end if
  end if

  ijac = index(name(1:nlen), jacobianstr)
  if (ijac == 1) then
    call abor1_ftn(trim(caller) // ": _jacobian_ must be preceded by dependent " // &
                   "variable in config: " // trim(name))
  else if (ijac > 1) then
    ! Jacobian: <ystr>_jacobian_<xstr>[_<chstr>]
    jacobian_needed = .true.
    yend = ijac - 1
    xbeg = ijac + len(jacobianstr)
    xend = nlen
    if (ch(jvar) /= no_channel) xend = last_underscore - 1
    if (xend < xbeg) then
      call abor1_ftn(trim(caller) // ": _jacobian_ must be followed by independent " // &
                     "variable in config: " // trim(name))
    end if
    call set_part(ystr(jvar), name(1:yend))
    call set_part(xstr(jvar), name(xbeg:xend))
  else
    ! forward diagnostic: <ystr>[_<ch>]
    yend = nlen
    if (ch(jvar) /= no_channel) yend = last_underscore - 1
    call set_part(ystr(jvar), name(1:yend))
  end if
end do

contains

  !> part = str, aborting if str does not fit
  subroutine set_part(part, str)
    character(len=MAXVARLEN), intent(out) :: part
    character(len=*),         intent(in)  :: str
    if (len_trim(str) > len(part)) then
      call abor1_ftn(trim(caller) // ": hofxdiags name part too long: " // trim(str))
    end if
    part = str
  end subroutine set_part

end subroutine ufo_parse_hofxdiags

! ------------------------------------------------------------------------------

end module ufo_hofxdiags_parse_mod
