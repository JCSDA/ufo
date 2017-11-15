module obs_read_mod

use iso_c_binding
use radDiag_mod, only: RadDiag
use radDiag_mod, only: RadDiag_read
use raobDiag_mod, only: RaobDiag
use raobDiag_mod, only: RaobDiag_read

implicit none
private
public :: obs_read_setup
public :: obs_read_delete

type all_obs_type
   type(RadDiag), pointer :: Rad
   type(RaobDiag),pointer :: Raob
end type all_obs_type

interface obs_read_setup; procedure setup_; end interface
interface obs_read_delete; procedure delete_; end interface

contains
   subroutine setup_(filein,obstype,nobs)
   character(len=*),intent(in) :: filein
   character(len=*),intent(in) :: obstype
   integer(c_int),intent(inout) :: nobs

   type(all_obs_type) :: obs_type

   if (trim(obstype) == "Radiance"   ) then
       allocate(obs_type%Rad)
       call radDiag_read (obs_type%Rad,filein,nobs)
   endif
   if (trim(obstype) == "Radiosonde" ) then
       allocate(obs_type%Raob)
      call raobDiag_read(obs_type%Raob,filein,nobs)
   endif
   end subroutine setup_

   subroutine delete_()
   type(all_obs_type) :: obs_type
   if(associated(obs_type%Rad))  nullify(obs_type%Rad)
   if(associated(obs_type%Raob)) nullify(obs_type%Raob)
   end subroutine delete_
end module obs_read_mod
