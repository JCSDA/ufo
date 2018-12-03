! (C) Copyright 2017-2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

module ufo_basis_mod

  use iso_c_binding
  use ufo_geovals_mod
  use ufo_geovals_mod_c,   only: ufo_geovals_registry
  use ufo_locs_mod
  use ufo_locs_mod_c, only : ufo_locs_registry
  use datetime_mod

  type, abstract :: ufo_basis
    private
  contains
    procedure, non_overridable :: opr_simobs  => opr_simobs_
    procedure(simobs_), deferred  :: simobs
    procedure, non_overridable :: opr_locateobs  => opr_locateobs_
    procedure :: locateobs => ufo_basis_locateobs !Overridable in extending types
!    procedure(locateobs_), deferred  :: locateobs 
  end type ufo_basis

  abstract interface

  ! ------------------------------------------------------------------------------

    subroutine simobs_(self, geovals, hofx, obss)
      use iso_c_binding
      import ufo_basis, ufo_geovals
      implicit none
      class(ufo_basis),         intent(in)    :: self
      type(ufo_geovals),        intent(in)    :: geovals
      real(c_double),           intent(inout) :: hofx(:)
      type(c_ptr), value,       intent(in)    :: obss
    end subroutine simobs_

  ! ------------------------------------------------------------------------------

    subroutine locateobs_(self, obss, t1, t2, locs)
!    subroutine locateobs_(obss, t1, t2, locs)
      use iso_c_binding
      use datetime_mod
      import ufo_basis, ufo_locs!, datetime_mod
      implicit none
      class(ufo_basis),         intent(in)    :: self
      type(c_ptr), value,       intent(in)    :: obss
      type(datetime),           intent(in)    :: t1, t2
      type(ufo_locs),           intent(inout) :: locs
    end subroutine locateobs_

  ! ------------------------------------------------------------------------------

  end interface
  
contains

! ------------------------------------------------------------------------------
    
  subroutine opr_simobs_(self, c_key_geovals, c_obsspace, c_hofx)
    implicit none

    class(ufo_basis), intent(in)    :: self
    integer(c_int), intent(in) :: c_key_geovals
    type(c_ptr), value, intent(in) :: c_obsspace
    real(c_double), intent(inout) :: c_hofx(:)

    type(ufo_geovals), pointer :: geovals

    call ufo_geovals_registry%get(c_key_geovals,geovals)

    call self%simobs(geovals, c_hofx, c_obsspace)

  end subroutine opr_simobs_

! ------------------------------------------------------------------------------
    
  subroutine opr_locateobs_(self, c_obsspace, c_t1, c_t2, c_locs)
    implicit none

    class(ufo_basis), intent(in)   :: self
    type(c_ptr), value, intent(in) :: c_obsspace
    type(c_ptr), intent(in)        :: c_t1, c_t2
    integer(c_int), intent(inout)  :: c_locs

    type(datetime) :: t1, t2
    type(ufo_locs), pointer :: locs

    call c_f_datetime(c_t1, t1)
    call c_f_datetime(c_t2, t2)

    call ufo_locs_registry%init()
    call ufo_locs_registry%add(c_locs)
    call ufo_locs_registry%get(c_locs,locs)

    call self%locateobs(c_obsspace, t1, t2, locs)

  end subroutine opr_locateobs_

! ------------------------------------------------------------------------------

subroutine ufo_basis_locateobs(self, obss, t1, t2, locs)
!subroutine ufo_basis_locateobs(obss, t1, t2, locs)
  use kinds
!  use datetime_mod
  use twindow_utils_mod
  use fckit_log_module, only : fckit_log
  use obsspace_mod

  implicit none

  class(ufo_basis), intent(in) :: self
  type(c_ptr), value, intent(in)              :: obss
  type(datetime), intent(in)                  :: t1, t2
  type(ufo_locs), intent(inout)               :: locs

  integer :: nlocs
  type(datetime) :: refdate

  character(len=*),parameter:: &
     myname = "ufo_basis_locateobs"
  character(len=255) :: record
  integer :: i
  integer :: tw_nlocs
  integer, dimension(:), allocatable :: tw_indx
  real(kind_real), dimension(:), allocatable :: time, lon, lat

  ! Local copies pre binning
  nlocs = obsspace_get_nlocs(obss)
  refdate = obsspace_get_refdate(obss)  !!FUNCTION
!  call obsspace_get_refdate(obss,refdate)  !!SUBROUTINE

  allocate(time(nlocs), lon(nlocs), lat(nlocs))

!TODO(JG): Add "MetaData" or similar attribute to all ObsSpace's 
!  call obsspace_get_db(obss, "MetaData", "time", time)
  call obsspace_get_db(obss, "", "time", time)

  ! Generate the timing window indices
  allocate(tw_indx(nlocs))
  call gen_twindow_index(refdate, t1, t2, nlocs, time, tw_indx, tw_nlocs)

  !!Each operator may have its own way to derive lon, lat from MetaData
  !!BEGIN THIS PART CAN BE UNIQUE FOR SOME OBS OPERATORS
!TODO(JG): Add "MetaData" or similar attribute to all ObsSpace's 
!  call obsspace_get_db(obss, "MetaData", "longitude", lon)
!  call obsspace_get_db(obss, "MetaData", "latitude", lat)
  call obsspace_get_db(obss, "", "longitude", lon)
  call obsspace_get_db(obss, "", "latitude", lat)

  !Setup ufo locations
  call ufo_locs_setup(locs, tw_nlocs)
  do i = 1, tw_nlocs
    locs%lon(i)  = lon(tw_indx(i))
    locs%lat(i)  = lat(tw_indx(i))
    locs%time(i) = time(tw_indx(i))
  enddo
  locs%indx = tw_indx(1:tw_nlocs)
  !!END THIS PART CAN BE UNIQUE FOR SOME OBS OPERATORS

  deallocate(time, lon, lat, tw_indx)

  write(record,*) myname,': allocated/assigned obs locations'
  call fckit_log%info(record)

end subroutine ufo_basis_locateobs

! ------------------------------------------------------------------------------

end module ufo_basis_mod
