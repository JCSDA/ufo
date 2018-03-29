!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    read_profiles                       read profiles from gsi
!   prgmmr: tahara           org: np20                date: 2003-01-01
!
!   2017-11-10 pagowski - converted radiance to aod 
!

MODULE read_profiles_mod

  USE ncd_kinds, ONLY:  i_kind,r_single,r_kind
  IMPLICIT NONE

! Declare public and private
  PRIVATE
  PUBLIC :: max_name_length,max_vars
  PUBLIC :: read_profiles_header, read_profiles

!same as in ufo
  INTEGER(i_kind),PARAMETER :: max_name_length=56,max_vars=50

CONTAINS

  SUBROUTINE read_profiles_header(ftin,idate,nsig,nvars,naeros,varnames,iflag,lverbose)
!                .      .    .                                       .
! Declare passed arguments
    INTEGER,INTENT(in)             :: ftin
    INTEGER,INTENT(out) :: idate,nsig,nvars,naeros
    CHARACTER(len=max_name_length), DIMENSION(max_vars) :: &
         &varnames
    INTEGER,INTENT(out)            :: iflag
    LOGICAL,OPTIONAL,INTENT(in)            :: lverbose    

    LOGICAL loutall

    loutall=.TRUE.
    IF(PRESENT(lverbose)) loutall=lverbose

! Read header (fixed_part).

    READ(ftin,IOSTAT=iflag)  nsig,nvars,naeros,idate
    READ(ftin,IOSTAT=iflag)  varnames(1:nvars)

  END SUBROUTINE read_profiles_header

  SUBROUTINE read_profiles(ftin,nsig,nvarsphys,naeros,&
       &tvp,qvp,prsltmp,prsitmp,aeros,iflag,lverbose)
!                .      .    .                                       .
! Declare passed arguments
    INTEGER,INTENT(in)             :: ftin
    INTEGER,INTENT(in) :: nsig,nvarsphys,naeros
    REAL(r_single), DIMENSION(nsig) :: tvp,qvp,prsltmp
    REAL(r_single), DIMENSION(nsig+1) :: prsitmp
    REAL(r_single), DIMENSION(nsig,naeros) :: aeros
    INTEGER,INTENT(out)            :: iflag
    LOGICAL,OPTIONAL,INTENT(in)            :: lverbose    

    LOGICAL loutall

    loutall=.TRUE.
    IF(PRESENT(lverbose)) loutall=lverbose

    READ(ftin,IOSTAT=iflag) tvp,qvp,prsltmp,prsitmp
    READ(ftin,IOSTAT=iflag) aeros

  END SUBROUTINE read_profiles

END MODULE read_profiles_mod

