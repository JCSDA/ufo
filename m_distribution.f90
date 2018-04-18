module type_distribution
   use fckit_module
   use fckit_log_module, only: fckit_log

   implicit none

   private
   public:: random_distribution

   type :: random_distribution
      type(fckit_mpi_comm),private :: comm    !> parallel communicator
      integer, private             :: nobs    !> the counter of the number of jobs 
      integer, private             :: myproc  !> size of parallel communicator
      integer, private             :: nproc   !> my rank
      integer, allocatable, public :: indx(:) !> distribution layout on each PE
   contains
      procedure :: received
      procedure :: reset
      procedure :: nobs_pe
      procedure :: mype
      final     :: destructor
   end type

   interface random_distribution
      module procedure constructor
   end interface

   contains

      type(random_distribution) function constructor(gnobs)
         integer, intent(in) :: gnobs

         ! Local variables
         integer :: info, i, ic
         logical :: init

         call fckit_log%debug('random_distribution object created')
         constructor%nobs=0
         
         ! Get the nproc and myproc
         constructor%nproc=constructor%comm%size()
         constructor%myproc=constructor%comm%rank()

         ! Count the obs. number on each PE
         do i=1,gnobs
            if (mod(i-1, constructor%nproc) .eq. constructor%myproc) then
               constructor%nobs=constructor%nobs+1
            end if
         end do

         allocate(constructor%indx(constructor%nobs))
        
         ! Create a default layout
         ic=0
         do i=1,gnobs
            if (mod(i-1, constructor%nproc) .eq. constructor%myproc) then
               ic=ic+1
               constructor%indx(ic)=i
            end if
         end do
      end function

      subroutine destructor(this)
         type(random_distribution), intent(inout):: this
         deallocate(this%indx) 
         call fckit_log%debug('random_distribution object finalized')
      end subroutine

      logical function received(this, seqno)
         class(random_distribution), intent(inout):: this
         integer, intent(in) :: seqno
       
         received=.False.
         !> seqno should start from 1
         if (mod(seqno-1, this%nproc) .eq. this%myproc) then
             this%nobs=this%nobs+1
             received=.True.
         end if
         call fckit_log%debug('distribute Obs Data Functionality')
      end function

      integer function nobs_pe(this)
         class(random_distribution), intent(in):: this
       
         call fckit_log%debug('nobs_pe Data Functionality')
         nobs_pe=this%nobs
      end function

      integer function mype(this)
         class(random_distribution), intent(in):: this

         mype=this%myproc
      end function

      subroutine reset(this)
         class(random_distribution), intent(inout):: this
       
         call fckit_log%debug('reset Data Functionality')
         this%nobs=0
      end subroutine

end module type_distribution
