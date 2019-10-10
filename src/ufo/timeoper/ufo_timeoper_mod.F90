! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module for timeoper observation operator

module ufo_timeoper_mod

use iso_c_binding
use fckit_configuration_module, only: fckit_configuration
use fckit_log_module, only : fckit_log
use kinds
use datetime_mod
use duration_mod

use ufo_geovals_mod, only: ufo_geovals, ufo_geoval, ufo_geovals_get_var
use ufo_basis_mod, only: ufo_basis
use ufo_vars_mod
use obsspace_mod

implicit none
private

!> Fortran derived type for the observation type
type, public :: ufo_timeoper
  integer                       :: c_size
  character(:), allocatable     :: optypename
  integer                       :: time_stencil
  real(kind_real), allocatable  :: time_weight(:)
contains
  procedure :: setup  => ufo_timeoper_setup
  procedure :: set_timeweight => ufo_timeoper_set_timeweight
  procedure :: simobs => ufo_timeoper_simobs
  final :: destructor
end type ufo_timeoper

contains

! ------------------------------------------------------------------------------
subroutine ufo_timeoper_setup(self, f_conf, c_size, c_tstencil)

implicit none
class(ufo_timeoper),      intent(inout) :: self
type(fckit_configuration), intent(in)   :: f_conf
integer,                     intent(in) :: c_size ! 1d obsspace vector length
integer,                     intent(in) :: c_tstencil

self%c_size = c_size
self%time_stencil = c_tstencil

call f_conf%get_or_die("name", self%optypename)

allocate(self%time_weight(c_size*self%time_stencil))
self%time_weight = 0.0

end subroutine ufo_timeoper_setup

! ------------------------------------------------------------------------------
subroutine ufo_timeoper_set_timeweight(self, f_conf, obss, t0, t3, state_time, &
                                       t_step)

implicit none
class(ufo_timeoper),      intent(inout) :: self
type(fckit_configuration), intent(in)   :: f_conf
type(c_ptr),  value, intent(in) :: obss
type(datetime),      intent(in) :: t0, t3, state_time
type(duration),      intent(in) :: t_step

! local variables
character(len=*),parameter    :: myname = "ufo_timeoper_locs_init"

integer :: i, tw_nlocs, nlocs
integer,         dimension(:), allocatable :: tw_indx, date_indx
real(kind_real), dimension(:), allocatable :: lon_in, lat_in
type(datetime),  dimension(:), allocatable :: date_time_in
type(duration) :: dt

! Local copies pre binning
nlocs = obsspace_get_nlocs(obss)

allocate(date_time_in(nlocs), lon_in(nlocs), lat_in(nlocs))

call obsspace_get_db(obss, "MetaData", "datetime", date_time_in)
call obsspace_get_db(obss, "MetaData", "longitude", lon_in)
call obsspace_get_db(obss, "MetaData", "latitude", lat_in)

allocate(tw_indx(nlocs))
allocate(date_indx(nlocs))
tw_nlocs = 0
do i = 1, nlocs
  if (date_time_in(i) >= t0 .and. date_time_in(i) < t3 &
     .and. date_time_in(i) >= state_time) then
    tw_nlocs = tw_nlocs + 1
    tw_indx(tw_nlocs) = 2 * i - 1
    date_indx(tw_nlocs) = i
  elseif  (date_time_in(i) >= t0 .and. date_time_in(i) < t3 &
    .and. date_time_in(i) < state_time) then
    tw_nlocs = tw_nlocs + 1
    tw_indx(tw_nlocs) = 2 * i
    date_indx(tw_nlocs) = i
  endif
enddo

do i = 1, tw_nlocs
  call datetime_diff(date_time_in(date_indx(i)), state_time, dt)
  self%time_weight(tw_indx(i)) = 1.0 - &
    abs(dble(duration_seconds(dt)) / &
    dble(duration_seconds(t_step)))
end do

end subroutine ufo_timeoper_set_timeweight


! ------------------------------------------------------------------------------
subroutine destructor(self)
implicit none
type(ufo_timeoper), intent(inout) :: self

if (allocated(self%time_weight)) deallocate(self%time_weight)

end subroutine destructor

! ------------------------------------------------------------------------------
subroutine ufo_timeoper_simobs(self, gv, obss)
implicit none
class(ufo_timeoper), intent(in) :: self
type(ufo_geovals),  intent(inout) :: gv
type(c_ptr), value, intent(in) :: obss

! Local variables
integer :: jl, jv, jo, nlocs, i, nval_tot
logical :: l_timeinterp

nlocs = obsspace_get_nlocs(obss)

! we need to be able not run time interpolation once it has already been done
! we can't do this by having a switch in the setup and that will not work
! in the tlad case. So the best solution is to count the number of missing-value
! indicators. If they are greater than half of the whole geovals then we switch
! the code off.
l_timeinterp = .true.
i = 0
nval_tot = 0
do jv = 1, gv%nvar
  do jo = 1, nlocs
    do jl = 1, gv%geovals(jv)%nval
      if (gv%geovals(jv)%vals(jl,nlocs+jo) <= gv%missing_value) i = i + 1
    end do
  end do
  nval_tot = nval_tot + nlocs * gv%geovals(jv)%nval
enddo
if (i >= nval_tot) l_timeinterp = .false.

if (l_timeinterp) then

  ! convolve with time-weights
  do jv = 1, gv%nvar
    do jo = 1, gv%nlocs
      do jl = 1, gv%geovals(jv)%nval
        gv%geovals(jv)%vals(jl,jo) = gv%geovals(jv)%vals(jl,jo) * &
                                     self%time_weight(jo)
      enddo
    enddo
  enddo

  ! time-interpolation
  do jo = 1, nlocs

    do jv = 1, gv%nvar
      gv%geovals(jv)%vals(:,jo) = &
        SUM(gv%geovals(jv)%vals(:, &
            self%time_stencil*(jo-1)+1:self%time_stencil*jo), DIM=2)
    enddo

  enddo

  !missing flag on rest of geovals
  do jv = 1, gv%nvar
    gv%geovals(jv)%vals(:,nlocs+1:self%time_stencil*nlocs) = gv%missing_value
  enddo

end if

end subroutine ufo_timeoper_simobs


! ------------------------------------------------------------------------------

end module ufo_timeoper_mod
