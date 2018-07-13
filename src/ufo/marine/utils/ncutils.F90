! (C) Copyright 2017-2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Fortran module to handle temperature profile observations

module ufo_marine_ncutils

  use ncd_kinds
  
  implicit none
  private

  ! General type definition for marine obs
  type, public :: simple_marine_obs
     integer      :: Station_ID
     real(r_kind) :: Observation_Type
     real(r_kind) :: Latitude
     real(r_kind) :: Longitude
     real(r_kind) :: Depth
     real(r_kind) :: Time
     real(r_kind) :: Observation
     real(r_kind) :: Observation_error     
     real(r_kind) :: Obs_Minus_Forecast
  end type simple_marine_obs

  type, public :: diag_marine_obs
     integer                              :: nobs           !< Number of obs
     character(len=120)                   :: filename       !< Netcdf output filename
     logical                              :: append=.false. !< If file exist, append to it
     type(simple_marine_obs), allocatable :: diag(:)        !< Data type to hold obs space diagnostics
   contains
     procedure :: init
     procedure :: write_diag
     procedure :: write_geoval     
     procedure :: finalize
  end type diag_marine_obs
  
contains

  ! ------------------------------------------------------------------------------

  subroutine init(self, nobs, filename)

    implicit none

    class(diag_marine_obs) ,intent(out) :: self     !< Obs space diagnostics
    integer                 ,intent(in) :: nobs     !< Number of obs    
    character(len=120)      ,intent(in) :: filename !< Filename for netcdf output
    
    self%nobs=nobs
    self%filename=filename
    allocate(self%diag(nobs))
    
  end subroutine init

  ! ------------------------------------------------------------------------------

  subroutine finalize(self)

    implicit none

    class(diag_marine_obs), intent(inout) :: self  !< Obs space diagnostics

    self%nobs=0
    self%filename=''
    deallocate(self%diag)
    
  end subroutine finalize

  ! ------------------------------------------------------------------------------

  subroutine write_diag(self)

    use netcdf
    use ufo_vars_mod

    implicit none

    class(diag_marine_obs), intent(inout) :: self  !< Obs space diagnostics

    !netcdf stuff
    integer(kind=4) :: iNcid,i,iStationNo
    integer(kind=4) :: iDimStation_ID,iDimTime_ID,iDimLenStringName_ID, iDimLev_ID
    integer(kind=4) :: iVarLON_ID, iVarLAT_ID, iVarLev_ID, iVarObs_ID
    integer(kind=4) :: iVarOMF_ID, iVarOMA_ID, iVarSIGO_ID, iVarOBSID_ID

    integer, allocatable :: obsid(:)
    integer :: nlev,nobs,iobs

    ! Create file.
    call check(nf90_create(self%filename, NF90_CLOBBER, iNcid))
    call check(nf90_def_dim(iNcid, "nobs", NF90_UNLIMITED, iDimStation_ID))
    nobs = self%nobs    

    ! Define of variables.
    call check( nf90_def_var(iNcid, "OBSID", NF90_INT, (/ iDimStation_ID /), iVarOBSID_ID) )
    call check( nf90_def_var(iNcid, "lon", NF90_REAL, (/ iDimStation_ID /), iVarLON_ID) )
    call check( nf90_def_var(iNcid, "lat", NF90_REAL, (/ iDimStation_ID /), iVarLAT_ID) )
    call check( nf90_def_var(iNcid, "lev", NF90_REAL, (/ iDimStation_ID /), iVarLev_ID) )
    call check( nf90_def_var(iNcid, "obs", NF90_REAL, (/ iDimStation_ID /), iVarObs_ID) )
    call check( nf90_def_var(iNcid, "sigo", NF90_REAL, (/ iDimStation_ID /), iVarSIGO_ID) )
    call check( nf90_def_var(iNcid, "omf", NF90_REAL, (/ iDimStation_ID /), iVarOMF_ID) )

    ! End define mode.
    call check(nf90_enddef(iNcid))

    ! Writing
    call check(nf90_put_var(iNcid, iVarLON_ID , self%diag(:)%Station_ID))
    call check(nf90_put_var(iNcid, iVarLON_ID , self%diag(:)%Longitude))
    call check(nf90_put_var(iNcid, iVarLAT_ID , self%diag(:)%Latitude))
    call check(nf90_put_var(iNcid, iVarLev_ID , self%diag(:)%Depth))
    call check(nf90_put_var(iNcid, iVarSIGO_ID , self%diag(:)%Observation_error))        
    call check(nf90_put_var(iNcid, iVarObs_ID , self%diag(:)%Observation))
    call check(nf90_put_var(iNcid, iVarOMF_ID , self%diag(:)%Obs_Minus_Forecast))

    ! Close file.
    call check(nf90_close(iNcid))
    self%append=.true.
  end subroutine write_diag

  ! ------------------------------------------------------------------------------

  subroutine write_geoval(self,varname,geoval,arg_dim_name)

    use netcdf
    use ufo_vars_mod
    use ufo_geovals_mod
    
    implicit none

    class(diag_marine_obs)             , intent(in) :: self     !< Obs space diagnostics
    character(len=MAXVARLEN)           , intent(in) :: varname  !< One of var_ from ufo_vars_mod
    type(ufo_geoval)          , pointer, intent(in) :: geoval   !< 2D array for 1 geoval
    character(len=MAXVARLEN), optional , intent(in) :: arg_dim_name  !< Name for the second dimension
    !netcdf stuff    
    integer(kind=4) :: iNcid
    integer(kind=4) :: iDimStation_ID, iDimLev_ID
    integer(kind=4) :: iVargeoval_ID
    integer :: nlev, ndims
    character(len=MAXVARLEN) :: dim_name
    
    dim_name="nlev"
    if (present(arg_dim_name)) dim_name=arg_dim_name
    
    if (self%append) then  ! If file exists, append to it
       call check( nf90_open(self%filename, NF90_WRITE, iNcid) )
       call check( nf90_inquire(iNcid, nDimensions = ndims) )
       call check( nf90_inq_dimid(iNcid, "nobs", iDimStation_ID) )
       if (ndims.eq.2) then
          call check( nf90_inq_dimid(iNcid, dim_name, iDimLev_ID) )
       end if
       call check( nf90_redef(iNcid) )
       if (ndims.eq.1) then
          nlev = geoval%nval          
          call check(nf90_def_dim(iNcid, dim_name, nlev, iDimLev_ID))
       end if
    else    
       call check(nf90_create(self%filename, NF90_CLOBBER, iNcid))
       call check(nf90_def_dim(iNcid, "nobs", NF90_UNLIMITED, iDimStation_ID))
       nlev = geoval%nval    
       call check(nf90_def_dim(iNcid, dim_name, nlev, iDimLev_ID))
    end if

    ! Define of variables.
    call check( nf90_def_var(iNcid, varname, NF90_REAL,  (/ iDimLev_ID, iDimStation_ID /), iVargeoval_ID) )

    ! End define mode.
    call check(nf90_enddef(iNcid))

    ! Writing
    call check(nf90_put_var(iNcid, iVargeoval_ID , geoval%vals))     

    ! Close file.
    call check(nf90_close(iNcid))
    
  end subroutine write_geoval

  ! ------------------------------------------------------------------------------  
  
  subroutine check(status)

    use netcdf
    IMPLICIT NONE
    integer(4), intent ( in) :: status
    if(status /= nf90_noerr) then
       print *, trim(nf90_strerror(status))
       stop "Stopped"
    end if
  end subroutine check

end module ufo_marine_ncutils
