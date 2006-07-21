! Data Assimilation Research Testbed -- DART
! Copyright 2004-2006, Data Assimilation Research Section
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

module mpi_utilities_mod

! <next five lines automatically updated by CVS, do not edit>
! $Source$ 
! $Revision$ 
! $Date$ 
! $Author$ 
! $Name$ 

!-----------------------------------------------------------------------------
!
!   A collection of interfaces to the MPI (Message Passing Interface)
!   multi-processor communication library routines.
!
!    # initialize_mpi_utilities()  Subroutine that initializes MPI and sets
!                                  local values needed later.  Must be called
!                                  before any other routine here.
!
!    # finalize_mpi_utilities()  Subroutine that shuts down MPI cleanly.
!                                Must be called before program exits, and no
!                                other routines here can be used afterwards.
!
!    # task_count()       Function that returns the total number of MPI tasks.
!
!    # my_task_id()       Function that returns my task number.  Note that
!                         in the MPI world task numbers run from 0 to N-1.
!
!    # send_to()          Subroutine which sends a 1D data array
!                         synchronously to another task (point-to-point).
!
!    # receive_from()     Subroutine which receives a 1D data array
!                         synchronously from another task (point-to-point).
!
!    # task_sync()        Subroutine that only returns after every task has
!                         reached this same location in the code.
!        
!    # array_broadcast()  Subroutine that sends a copy of the entire data 
!                         array to all other tasks. 
!                
!    * transpose_array()  Subroutine that transposes a 2D array
!                         from column-major to row-major or back.
!
!   ** array_distribute() Subroutine that distributes a data array across the
!                         other tasks, so each task gets a non-overlapping 
!                         subset of the data.
!                
!   MPI cover routines more specific for DART and hopefully more useful.
!
!    # iam_task0()        Function which returns .TRUE. if task id is 0,
!                         .FALSE. for anything else.
!
!    # broadcast_send()   Subroutine which takes two r8 arrays and broadcasts
!                         them to all other tasks.  If sending ID is not the
!                         same as the local task ID, an error is returned.
!                         Does not return until all other tasks have called
!                         recv to pick up the data.
!
!    # broadcast_recv()   Subroutine which receives two r8 arrays from the 
!                         sending task ID.  If the sending ID is the same as
!                         the local task ID, an error is returned.  All other
!                         tasks must call recv before this routine returns.
!
!    * sum_across_tasks() Subroutine which takes a single integer argument
!                         from each task, and returns the sum of all integers
!                         across all tasks back to all tasks.  All tasks must
!                         call this routine before it can compute and return
!                         the value.
!
!   Lower level utility routines which interact with the utilities_mod.f90 
!   code to open a named pipe per MPI task, read and write from them, and
!   close and/or remove them.
!
!  *** make_pipe()        Function that creates a named pipe (fifo), opens it,
!                         and returns the unit number.  Ok to call if the pipe
!                         already exists or is already open; it will skip
!                         those steps and just return the unit number.  The 
!                         name argument is used as a base and a filename
!                         in the form 'base.NNNN' is generated, where the N's
!                         are the MPI rank number, 0 padded.
!
!  *** destroy_pipe()     The unit number is closed and the pipe file is 
!                         removed.
!
!    * read_pipe()        The character string is read from the pipe.
!                         (Can be overloaded to read ints if time or status
!                         info is useful to exchange between processes.) 
!                         This routine blocks until data is available.
!
!    * write_pipe()       The character string is written to the pipe.
!                         (Can be overloaded to write ints if time or status
!                         info is useful to exchange between processes.) 
!                         This routine writes and returns immediately.
!
!   Wrappers for system functions.  Covers differences if you run with
!   or without MPI.
!
!  *** shell_execute()    Use the system() command to execute a command string.
!                         Will wait for the command to complete and returns an
!                         error code unless you end the command with & to put
!                         it into background.   Function which returns the rc
!                         of the command, 0 being all is ok.
!
!  *** sleep_seconds()    Wrapper for the sleep command.  Argument is a real
!                         in seconds.  Different systems have different lower
!                         resolutions for the minimum time it will sleep.
!                         Subroutine, no return value.
!
!
!   # code done and tested
! *** both code and interface are done (but untested so far)
!  ** interface with proposed arguments exists but code not complete
!   * interface name only; no arg list devised yet 
!
!-----------------------------------------------------------------------------
! 
! these do not exist - i believe a single transpose will work.  but if not,
! they can be separated into these two, which can either work on a real
! 2D array or a single linearized array which is logically 2D but in reality
! stored in a 1D fortran array:
!
!      transpose_row_major()  Subroutine which transposes a logical 2D array
!                             from column-major to row-major.  The source and
!                             destination arrays must be stored in 1D arrays 
!                             of length (nrows * ncols).
!
!      transpose_col_major()  Subroutine which transposes a logical 2D array
!                             from row-major to column-major.  The source and
!                             destination arrays must be stored in 1D arrays 
!                             of length (nrows * ncols).
!
!-----------------------------------------------------------------------------

use types_mod, only : r8
use utilities_mod, only : register_module, error_handler, & 
                          E_ERR, E_WARN, E_MSG, E_DBG, get_unit, close_file, &
                          do_output
use time_manager_mod, only : time_type, get_time, set_time

!
! Some MPI installations have an MPI module; if one is present, use that.
! If not, there will be an MPI include file which defines the parameters.
! Use one but not both.  For help on compiling a module which uses MPI
! see the $DART/doc/mpi directory.

!use mpi

implicit none
private

include "mpif.h"


!   ---- private data for mpi_utilities ----

integer :: myrank          ! my mpi number
integer :: total_tasks     ! total mpi tasks/procs
integer :: my_local_comm   ! duplicate communicator private to this file
integer :: comm_size       ! if ens count < tasks, only the first N participate


public :: task_count, my_task_id, transpose_array, &
          initialize_mpi_utilities, finalize_mpi_utilities, &
          make_pipe, destroy_pipe
public :: task_sync, array_broadcast, array_distribute, &
          send_to, receive_from, iam_task0, broadcast_send, broadcast_recv, &
          shell_execute, sleep_seconds, sum_across_tasks

! CVS Generated file description for error handling, do not edit
character(len=128) :: &
source   = "$Source$", &
revision = "$Revision$", &
revdate  = "$Date$"

logical, save :: module_initialized = .false.

character(len = 129) :: errstring

! Namelist input - placeholder for now.
!namelist /mpi_utilities_nml/ x

contains

!-----------------------------------------------------------------------------
! mpi cover routines
!-----------------------------------------------------------------------------

subroutine initialize_mpi_utilities()

! Initialize MPI and query it for global information.  Make a duplicate
! communicator so that any user code which wants to call MPI will not 
! interfere with any outstanding asynchronous requests, accidental tag
! matches, etc.  This routine must be called before any other routine in
! this file, and it should not be called more than once (but it does have
! defensive code in case that happens.)

integer :: errcode
logical :: already

if ( module_initialized ) then
   ! return without calling the code below multiple times
   write(errstring, *) 'initialize_mpi_utilities has already been called'
   call error_handler(E_WARN,'initialize_mpi_utilities', errstring, source, revision, revdate)
   return
endif

if ( .not. module_initialized ) then
   ! Initialize the module with utilities
   call register_module(source, revision, revdate)
   module_initialized = .true.
endif

! allow for the possibility that the user code has already initialized
! the MPI libs and do not try to call initialize again.

call MPI_Initialized(already, errcode)
if (errcode /= MPI_SUCCESS) then
   write(errstring, '(a,i8)') 'MPI_Initialized returned error code ', errcode
   call error_handler(E_ERR,'initialize_mpi_utilities', errstring, source, revision, revdate)
endif

! generally this will be the first use of mpi and this code will be executed.
if (.not.already) then
   errcode = -999
   call MPI_Init(errcode)
   if (errcode /= MPI_SUCCESS) then
      write(errstring, '(a,i8)') 'MPI_Init returned error code ', errcode
      call error_handler(E_ERR,'initialize_mpi_utilities', errstring, source, revision, revdate)
   endif
endif

! duplicate the world communicator to isolate us from any other user
! calls to MPI.  All subsequent mpi calls here will use the local communicator
! and not the global world comm.
call MPI_Comm_dup(MPI_COMM_WORLD, my_local_comm, errcode)
if (errcode /= MPI_SUCCESS) then
   write(errstring, '(a,i8)') 'MPI_Comm_dup returned error code ', errcode
   call error_handler(E_ERR,'initialize_mpi_utilities', errstring, source, revision, revdate)
endif

! find out who we are (0 to N-1).
call MPI_Comm_rank(my_local_comm, myrank, errcode)
if (errcode /= MPI_SUCCESS) then
   write(errstring, '(a,i8)') 'MPI_Comm_rank returned error code ', errcode
   call error_handler(E_ERR,'initialize_mpi_utilities', errstring, source, revision, revdate)
endif

! number of tasks (if 10, returns 10.  task id numbers go from 0-9)
call MPI_Comm_size(my_local_comm, total_tasks, errcode)
if (errcode /= MPI_SUCCESS) then
   write(errstring, '(a,i8)') 'MPI_Comm_size returned error code ', errcode
   call error_handler(E_ERR,'initialize_mpi_utilities', errstring, source, revision, revdate)
endif

! TODO: if there are fewer ensembles than tasks, all the collective routines
! need to take that into account and not participate if they are > comm_size.
comm_size = total_tasks

! Turn off non-critical log messages from all but task 0, for performance
! TODO: this should be controlled by a namelist option, to enable selected
!       mpi tasks to default to printing for debugging.
if (myrank /= 0) call do_output(.FALSE.)

! MPI successfully initialized.

end subroutine initialize_mpi_utilities

!-----------------------------------------------------------------------------

subroutine finalize_mpi_utilities(callfinalize)
 logical, intent(in), optional :: callfinalize

! Shut down MPI cleanly.  This must be done before the program exits; on
! some implementations of MPI the final I/O flushes are not done until this
! is called.  The optional argument can prevent us from calling MPI_Finalize,
! so that user code can continue to use MPI after this returns.  Calling other
! routines in this file after calling finalize will invalidate your warranty.

integer :: errcode
logical :: dofinalize

if ( .not. module_initialized ) then
   write(errstring, *) 'initialize_mpi_utilities() must be called first'
   call error_handler(E_ERR,'finalize_mpi_utilities', errstring, source, revision, revdate)
endif

! Release the private communicator we created at init time.
call MPI_Comm_free(my_local_comm, errcode)
if (errcode /= MPI_SUCCESS) then
   write(errstring, '(a,i8)') 'MPI_Comm_free returned error code ', errcode
   call error_handler(E_ERR,'finalize_mpi_utilities', errstring, source, revision, revdate)
endif

! If the optional argument is not given, or is given and is true, 
! shut down mpi.  Only if the argument is specified and is false do we
! skip the finalization.
if (.not.present(callfinalize)) then
   dofinalize = .TRUE.
else
   dofinalize = callfinalize
endif

! Normally we shut down MPI here.  If the user tells us not to shut down MPI
! they must call this routine from their own code before exiting.
if (dofinalize) then
   call MPI_Finalize(errcode)
   if (errcode /= MPI_SUCCESS) then
      write(errstring, '(a,i8)') 'MPI_Finalize returned error code ', errcode
      call error_handler(E_ERR,'finalize_mpi_utilities', errstring, source, revision, revdate)
   endif
endif


end subroutine finalize_mpi_utilities


!-----------------------------------------------------------------------------

function task_count()

! Return the total number of MPI tasks.  e.g. if the number of tasks is 4,
! it returns 4.  (The actual task numbers are 0-3.)

integer :: task_count

if ( .not. module_initialized ) then
   write(errstring, *) 'initialize_mpi_utilities() must be called first'
   call error_handler(E_ERR,'task_count', errstring, source, revision, revdate)
endif

task_count = total_tasks

end function task_count


!-----------------------------------------------------------------------------

function my_task_id()

! Return my unique task id.  Values run from 0 to N-1 (where N is the
! total number of MPI tasks.

integer :: my_task_id

if ( .not. module_initialized ) then
   write(errstring, *) 'initialize_mpi_utilities() must be called first'
   call error_handler(E_ERR,'my_task_id', errstring, source, revision, revdate)
endif

my_task_id = myrank

end function my_task_id


!-----------------------------------------------------------------------------

subroutine task_sync()

! Synchronize all tasks.  This subroutine does not return until all tasks
! execute this line of code.

integer :: errcode

if ( .not. module_initialized ) then
   write(errstring, *) 'initialize_mpi_utilities() must be called first'
   call error_handler(E_ERR,'task_sync', errstring, source, revision, revdate)
endif

call MPI_Barrier(my_local_comm, errcode)
if (errcode /= MPI_SUCCESS) then
   write(errstring, '(a,i8)') 'MPI_Barrier returned error code ', errcode
   call error_handler(E_ERR,'task_sync', errstring, source, revision, revdate)
endif

end subroutine task_sync


!-----------------------------------------------------------------------------

subroutine send_to(dest_id, srcarray, time)
 integer, intent(in) :: dest_id
 real(r8), intent(in) :: srcarray(:)
 type(time_type), intent(in), optional :: time

! Send the srcarray to the destination id.
! If time is specified, it is also sent in a separate communications call.  
! This is a synchronous call; it will not return until the destination has 
! called receive to accept the data.  If the send_to/receive_from calls are 
! not paired correctly the code will hang.

integer :: i, tag, errcode
integer :: datasize
integer :: itime(2)

if ( .not. module_initialized ) then
   write(errstring, *) 'initialize_mpi_utilities() must be called first'
   call error_handler(E_ERR,'send_to', errstring, source, revision, revdate)
endif

! simple idiotproofing
if ((dest_id < 0) .or. (dest_id >= total_tasks)) then
   write(errstring, '(a,i8,a,i8)') "destination task id ", dest_id, &
                                   "must be >= 0 and < ", total_tasks
   call error_handler(E_ERR,'send_to', errstring, source, revision, revdate)
endif

!print *, "kind = ", kind(srcarray(1))
!print *, "digits = ", digits(srcarray(1))
! TODO: FIX THIS based on what kind and digits say
datasize = MPI_DOUBLE_PRECISION

! use my task id as the tag; unused at this point.
tag = myrank

! call MPI to send the data to the remote task
call MPI_Ssend(srcarray, size(srcarray), datasize, dest_id, tag, &
              my_local_comm, errcode)
if (errcode /= MPI_SUCCESS) then
   write(errstring, '(a,i8)') 'MPI_Ssend returned error code ', errcode
   call error_handler(E_ERR,'send_to', errstring, source, revision, revdate)
endif

! if time specified, call MPI again to send the 2 time ints.
if (present(time)) then
   call get_time(time, itime(1), itime(2))
   call MPI_Ssend(itime, 2, MPI_INTEGER, dest_id, tag*2, my_local_comm, errcode)
   if (errcode /= MPI_SUCCESS) then
      write(errstring, '(a,i8)') 'MPI_Ssend returned error code ', errcode
      call error_handler(E_ERR,'send_to', errstring, source, revision, revdate)
   endif
endif


end subroutine send_to


!-----------------------------------------------------------------------------

subroutine receive_from(src_id, destarray, time)
 integer, intent(in) :: src_id
 real(r8), intent(out) :: destarray(:)
 type(time_type), intent(out), optional :: time

! Receive data into the destination array from the src task.
! If time is specified, it is received in a separate communications call.  
! This is a synchronous call; it will not return until the source has 
! sent the data.  If the send_to/receive_from calls are not paired correctly 
! the code will hang.

integer :: i, tag, errcode
integer :: datasize
integer :: itime(2)
real(r8), pointer :: temparray(:)
integer :: status(MPI_STATUS_SIZE)

if ( .not. module_initialized ) then
   write(errstring, *) 'initialize_mpi_utilities() must be called first'
   call error_handler(E_ERR,'receive_from', errstring, source, revision, revdate)
endif

! simple idiotproofing
if ((src_id < 0) .or. (src_id >= total_tasks)) then
   write(errstring, '(a,i8,a,i8)') "source task id ", src_id, &
                                   "must be >= 0 and < ", total_tasks
   call error_handler(E_ERR,'receive_from', errstring, source, revision, revdate)
endif

!print *, "kind = ", kind(destarray(1))
!print *, "digits = ", digits(destarray(1))
! TODO: FIX THIS based on what kind and digits say
datasize = MPI_DOUBLE_PRECISION

! send_to uses its own id as the tag.
tag = src_id

! call MPI to receive the data from the remote task
call MPI_Recv(destarray, size(destarray), datasize, src_id, MPI_ANY_TAG, &
              my_local_comm, status, errcode)
if (errcode /= MPI_SUCCESS) then
   write(errstring, '(a,i8)') 'MPI_Recv returned error code ', errcode
   call error_handler(E_ERR,'receive_from', errstring, source, revision, revdate)
endif

! if time specified, call MPI again to send the 2 time ints.
if (present(time)) then
   call MPI_Recv(itime, 2, MPI_INTEGER, src_id, tag*2, &
                 my_local_comm, status, errcode)
   if (errcode /= MPI_SUCCESS) then
      write(errstring, '(a,i8)') 'MPI_Recv returned error code ', errcode
      call error_handler(E_ERR,'receive_from', errstring, source, revision, revdate)
   endif
   if (itime(2) < 0) then
      write(errstring, '(a,i8)') 'seconds in settime were < 0; using 0'
      call error_handler(E_MSG,'receive_from', errstring, source, revision, revdate)
      time = set_time(itime(1), 0)
   else
      time = set_time(itime(1), itime(2))
   endif
endif


end subroutine receive_from


!-----------------------------------------------------------------------------
! TODO: do i need to overload this for both integer and real?
!       do i need to handle 1D, 2D, 3D inputs?


subroutine transpose_array

! not implemented here yet.  will have arguments -- several of them.

if ( .not. module_initialized ) then
   write(errstring, *) 'initialize_mpi_utilities() must be called first'
   call error_handler(E_ERR,'transpose_array', errstring, source, revision, revdate)
endif

write(errstring, *) 'not implemented yet'
call error_handler(E_ERR,'transpose_array', errstring, source, revision, revdate)

end subroutine transpose_array


!-----------------------------------------------------------------------------
! TODO: do i need to overload this for both integer and real?
!       do i need to handle 2D inputs?

subroutine array_broadcast(array, root)
 real(r8), intent(inout) :: array(:)
 integer, intent(in) :: root

! The data array values on the root task will be broadcast to every other
! task.  When this routine returns, all tasks will have the contents of the
! root array in their own arrays.  Thus 'array' is intent(in) on root, and
! intent(out) on all other tasks.

integer :: itemcount, datasize, errcode

if ( .not. module_initialized ) then
   write(errstring, *) 'initialize_mpi_utilities() must be called first'
   call error_handler(E_ERR,'array_broadcast', errstring, source, revision, revdate)
endif

! simple idiotproofing
if ((root < 0) .or. (root >= total_tasks)) then
   write(errstring, '(a,i8,a,i8)') "root task id ", root, &
                                   "must be >= 0 and < ", total_tasks
   call error_handler(E_ERR,'array_broadcast', errstring, source, revision, revdate)
endif

itemcount = size(array)
!print *, "kind = ", kind(array(1))
!print *, "digits = ", digits(array(1))
! TODO: FIX THIS based on what kind and digits say
datasize = MPI_DOUBLE_PRECISION

call MPI_Bcast(array, itemcount, datasize, root, my_local_comm, errcode)
if (errcode /= MPI_SUCCESS) then
   write(errstring, '(a,i8)') 'MPI_Bcast returned error code ', errcode
   call error_handler(E_ERR,'array_broadcast', errstring, source, revision, revdate)
endif

end subroutine array_broadcast


!-----------------------------------------------------------------------------
! TODO: do i need to overload this for both integer and real?
!       do i need to handle 2D inputs?

subroutine array_distribute(srcarray, root, dstarray, dstcount, how, which)
 real(r8), intent(in) :: srcarray(:)
 integer, intent(in) :: root
 real(r8), intent(out) :: dstarray(:)
 integer, intent(out) :: dstcount
 integer, intent(in) :: how
 integer, intent(out) :: which(:)

! 'srcarray' on the root task will be distributed across all the tasks
! into 'dstarray'.  dstarray must be large enough to hold each task's share
! of the data.  The actual number of values returned on each task will be
! passed back in the 'count' argument.  'how' is a flag to select how to
! distribute the data (round-robin, contiguous chunks, etc).  'which' is an
! integer index array which lists which of the original values were selected
! and put into 'dstarray'.

real(r8), allocatable :: localchunk(:)
integer :: srccount, datasize, leftover
integer :: i, tag, errcode
logical :: iamroot
integer :: status(MPI_STATUS_SIZE)

if ( .not. module_initialized ) then
   write(errstring, *) 'initialize_mpi_utilities() must be called first'
   call error_handler(E_ERR,'array_distribute', errstring, source, revision, revdate)
endif

! simple idiotproofing
if ((root < 0) .or. (root >= total_tasks)) then
   write(errstring, '(a,i8,a,i8)') "root task id ", root, &
                                   "must be >= 0 and < ", total_tasks
   call error_handler(E_ERR,'array_broadcast', errstring, source, revision, revdate)
endif

iamroot = (root == myrank)
tag = 1

srccount = size(srcarray)
print *, "kind = ", kind(srcarray(1))
print *, "digits = ", digits(srcarray(1))
! TODO: FIX THIS based on what kind and digits say
datasize = MPI_DOUBLE_PRECISION

! TODO: right now this code does contig chunks only
! TODO: it should select on the 'how' argument
dstcount = srccount / total_tasks
leftover = srccount - (dstcount * total_tasks)
if (myrank == total_tasks-1) dstcount = dstcount + leftover


! idiotproofing, continued...
if (size(dstarray) < dstcount) then
   write(errstring, '(a,i8,a,i8)') "size of dstarray is", size(dstarray), & 
                      " but must be >= ", dstcount
   call error_handler(E_ERR,'array_broadcast', errstring, source, revision, revdate)
endif
if (size(which) < dstcount) then
   write(errstring, '(a,i8,a,i8)') "size of which is", size(which), & 
                      " but must be >= ", dstcount
   call error_handler(E_ERR,'array_broadcast', errstring, source, revision, revdate)
endif

! TODO: this code is separate from the 'dstcount' computation because we
! need to test to be sure the user has passed us in arrays large enough to
! hold the data, but then this section needs to have a select (how) and set
! the corresponding index numbers accordingly.
which(1:dstcount) = (/ (i, i= myrank *dstcount, (myrank+1)*dstcount - 1) /)
if (size(which) > dstcount) which(dstcount+1:) = -1
   

if (.not.iamroot) then

   ! my task is receiving data.
   call MPI_Recv(dstarray, dstcount, datasize, root, MPI_ANY_TAG, &
	          my_local_comm, status, errcode)
   if (errcode /= MPI_SUCCESS) then
      write(errstring, '(a,i8)') 'MPI_Recv returned error code ', errcode
      call error_handler(E_ERR,'array_broadcast', errstring, source, revision, revdate)
   endif

else
   ! my task must send to everyone else and copy to myself.
   allocate(localchunk(dstcount), stat=errcode)  
   if (errcode /= 0) then
      write(errstring, *) 'allocation error of allocatable array'
      call error_handler(E_ERR,'array_broadcast', errstring, source, revision, revdate)
   endif

   do i=0, total_tasks-1
      ! copy correct data from srcarray to localchunk for each destination
      if (i == myrank) then
         ! this is my task, so do a copy from localchunk to dstarray
         dstarray(1:dstcount) = localchunk(1:dstcount)
      else
         ! call MPI to send the data to the remote task
         call MPI_Ssend(localchunk, dstcount, datasize, i, tag, &
   	               my_local_comm, errcode)
         if (errcode /= MPI_SUCCESS) then
            write(errstring, '(a,i8)') 'MPI_Ssend returned error code ', errcode
            call error_handler(E_ERR,'array_broadcast', errstring, source, revision, revdate)
         endif
      endif
      tag = tag + 1
   enddo

   deallocate(localchunk, stat=errcode)  
   if (errcode /= 0) then
      write(errstring, *) 'deallocation error of allocatable array'
      call error_handler(E_ERR,'array_broadcast', errstring, source, revision, revdate)
   endif
endif
   
! set any additional space which wasn't filled with zeros.
if (size(dstarray) > dstcount) dstarray(dstcount+1:) = 0.0

end subroutine array_distribute

!-----------------------------------------------------------------------------
! DART-specific cover utilities
!-----------------------------------------------------------------------------

function iam_task0()

! Return .TRUE. if my local task id is 0, .FALSE. otherwise.
! (Task numbers in MPI start at 0, contrary to the rules of polite fortran.)

logical :: iam_task0

if ( .not. module_initialized ) then
   write(errstring, *) 'initialize_mpi_utilities() must be called first'
   call error_handler(E_ERR,'iam_task0', errstring, source, revision, revdate)
endif

iam_task0 = (myrank == 0)

end function iam_task0

!-----------------------------------------------------------------------------
subroutine broadcast_send(from, array1, array2)
 integer, intent(in) :: from
 ! really only intent(in) here, but must match array_broadcast() call.
 real(r8), intent(inout) :: array1(:), array2(:)

! cover routine for array broadcast.  one additional sanity check -- make 
! sure the 'from' matches my local task id.  also, these arrays are
! intent(in) here, but they call a routine which is intent(inout) so they
! must be the same here.

if ( .not. module_initialized ) then
   write(errstring, *) 'initialize_mpi_utilities() must be called first'
   call error_handler(E_ERR,'broadcast_send', errstring, source, revision, revdate)
endif

! simple idiotproofing
if (from /= myrank) then
   write(errstring, '(a,i8,a,i8)') "'from' task id ", from, &
                                   "must be same as current task id ", myrank
   call error_handler(E_ERR,'broadcast_send', errstring, source, revision, revdate)
endif

! this must be paired with broadcast_recv() on all other tasks. 
! it will not return until all tasks in the communications group have
! made the call.
call array_broadcast(array1, from)
call array_broadcast(array2, from)

end subroutine broadcast_send

!-----------------------------------------------------------------------------
subroutine broadcast_recv(from, array1, array2)
 integer, intent(in) :: from
 ! really only intent(out) here, but must match array_broadcast() call.
 real(r8), intent(inout) :: array1(:), array2(:)

! cover routine for array broadcast.  one additional sanity check -- make 
! sure the 'from' is not the same as my local task id.  these arrays are
! intent(out) here, but they call a routine which is intent(inout) so they
! must be the same here.

if ( .not. module_initialized ) then
   write(errstring, *) 'initialize_mpi_utilities() must be called first'
   call error_handler(E_ERR,'broadcast_recv', errstring, source, revision, revdate)
endif

! simple idiotproofing
if (from == myrank) then
   write(errstring, '(a,i8,a,i8)') "'from' task id ", from, &
                                   "cannot be same as current task id ", myrank
   call error_handler(E_ERR,'broadcast_recv', errstring, source, revision, revdate)
endif

! this must be paired with a single broadcast_send() on the 'from' task.
! it will not return until all tasks in the communications group have
! made the call.
call array_broadcast(array1, from)
call array_broadcast(array2, from)

end subroutine broadcast_recv

!-----------------------------------------------------------------------------
subroutine sum_across_tasks(addend, sum)
 integer, intent(in) :: addend
 integer, intent(out) :: sum

 integer :: errcode
 integer :: localaddend(1), localsum(1)

! cover routine for MPI all-reduce

if ( .not. module_initialized ) then
   write(errstring, *) 'initialize_mpi_utilities() must be called first'
   call error_handler(E_ERR,'sum_across_tasks', errstring, source, revision, revdate)
endif

localaddend(1) = addend

call MPI_Allreduce(localaddend, localsum, 1, MPI_INTEGER, MPI_SUM, &
                   my_local_comm, errcode)
if (errcode /= MPI_SUCCESS) then
   write(errstring, '(a,i8)') 'MPI_Allreduce returned error code ', errcode
   call error_handler(E_ERR,'sum_across_tasks', errstring, source, revision, revdate)
endif

sum = localsum(1)

end subroutine sum_across_tasks


!-----------------------------------------------------------------------------
! pipe utilities
!-----------------------------------------------------------------------------
!    * make_pipe()        Function that creates a named pipe (fifo), opens it,
!                         and returns the unit number.  Ok to call if the pipe
!                         already exists or is already open; it will skip
!                         those steps and just return the unit number.  The 
!                         name argument is used as a base and a filename
!                         in the form 'base.NNNN' is generated, where the N's
!                         are the MPI rank number, 0 padded.
!
function make_pipe(pipename, exists) result (iunit)
 character(len=*), intent(in) :: pipename
 logical, intent(in), optional :: exists
 integer :: iunit

 ! forward declaration of external function
 interface
  function system(string)
   character(len=*) :: string
   integer :: system
  end function system
 end interface

! Create, open, and return a fortran unit number for a named pipe.
! The local MPI rank number will be appended to the given name to create
! a file of the form 'base.NNNN', where N's are the MPI rank number, 0 padded.
! TODO: based on the total number of tasks get extra style points for
! creating the shortest name necessary; e.g. base.N, base.NN, base.NNN, etc.
!
! If the optional 'exists' flag is not present, then it is not an error
! whether the pipe already exists or not.  It is made if it does not exist, 
! it is opened if not already opened, and the fortran unit number is returned.
! If 'exists' is present then it forces the issue of whether the pipe file
! must exist already or not.  The error handler is called if things aren't 
! as expected.  apologies to tim hoar for the intentional imitation of
! the generic file utilities_mod.f90 code.

logical :: open, there
character(len=128) :: fname
character(len=11) :: format
integer :: rc

if ( .not. module_initialized ) then
   write(errstring, *) 'initialize_mpi_utilities() must be called first'
   call error_handler(E_ERR,'make_pipe', errstring, source, revision, revdate)
endif

write(fname, "(a,i4.4)") trim(pipename)//".", myrank
print *, "fname now = ", trim(fname)

! check to see if the pipe already exists; if so, we've got the unit number
! (directly into the output value) and we're done.  otherwise, make it and
! open it.
inquire (file=fname, exist=there, opened=open, number=iunit, form=format)

if (.not. open) then

   if (.not. there) then
      ! make pipe file; mkfifo should be standard on any unix/linux system.
      rc = system('mkfifo '//trim(fname)//' '//char(0))

      ! and check to be sure it was made
      inquire (file=fname, exist=there)

      if (.not. there) then
        write(errstring, *) 'mkfifo command failed to create '//trim(fname)
        call error_handler(E_ERR,'make_pipe', errstring, source, revision, revdate)
      endif
   endif

   ! open pipe using an available unit number
   iunit = get_unit()
   open(unit=iunit, file=fname)

endif

! iunit contains the function return value.

end function make_pipe


!-----------------------------------------------------------------------------
!    * destroy_pipe()     The unit number is closed and the pipe file is 
!                         removed.
!
subroutine destroy_pipe(iunit)
 integer, intent(in) :: iunit

 ! forward declaration of external function
 interface
  function system(string)
   character(len=*) :: string
   integer :: system
  end function system
 end interface

character(len=128) :: pipename
integer :: ios, rc

if ( .not. module_initialized ) then
   write(errstring, *) 'initialize_mpi_utilities() must be called first'
   call error_handler(E_ERR,'destroy_pipe', errstring, source, revision, revdate)
endif

write(errstring, *) 'not implemented yet'
call error_handler(E_ERR,'destroy_pipe', errstring, source, revision, revdate)


! general idea is:

! call inquire to get name
inquire(unit=iunit, name=pipename, iostat=ios)
if (ios /= 0) then
   write(errstring, '(a,i4)') 'failure trying to inquire about unit ', iunit
   call error_handler(E_ERR,'destroy_pipe', errstring, source, revision, revdate)
endif

call close_file(iunit)

! remove echo when we trust this command.
rc = system('echo rm -f '//trim(pipename)//' '//char(0))


end subroutine destroy_pipe

!-----------------------------------------------------------------------------
!    * read_pipe()        The character string is read from the pipe.
!                         (Can be overloaded to read ints if time or status
!                         info is useful to exchange between processes.) 
!                         This routine blocks until data is available.
!
subroutine read_pipe()
end subroutine

!-----------------------------------------------------------------------------
!    * write_pipe()       The character string is written to the pipe.
!                         (Can be overloaded to write ints if time or status
!                         info is useful to exchange between processes.) 
!                         This routine writes and returns immediately.
!
subroutine write_pipe()
end subroutine


!-----------------------------------------------------------------------------
! general system util wrappers.
!-----------------------------------------------------------------------------
function shell_execute(execute_string, serialize)
 character(len=*), intent(in) :: execute_string
 logical, intent(in), optional :: serialize
 integer :: shell_execute

 ! forward declaration of external function
 interface
  function system(string)
   character(len=*) :: string
   integer :: system
  end function system
 end interface

! Use the system() command to execute a command string.
! Will wait for the command to complete and returns an
! error code unless you end the command with & to put
! it into background.   Function which returns the rc
! of the command, 0 being all is ok.

! allow code to test the theory that maybe the system call is
! not reentrant on some platforms.  if serialize is set and
! is true, do each call serially.

character(len=255) :: doit
logical :: ripit
integer :: i, errcode, dummy
integer :: status(MPI_STATUS_SIZE)

   ! default to everyone running concurrently, but if set and not true,
   ! serialize the calls to system() so they do not step on each other.
   if (present(serialize)) then
      ripit = .not. serialize
   else
      ripit = .TRUE.
   endif

   print *, "in-string is: ", trim(execute_string)

   !write(doit, "(a, 1x, a, 1x, a, a1)") "sh -c '", trim(execute_string), &
   !            "' < /dev/null >> fred.out", char(0)

   write(doit, "(a, 1x, a1)") trim(execute_string), char(0)

   print *, "about to run: ", trim(doit)
   print *, "input string length = ", len(trim(doit))

   if (ripit) then
      shell_execute = system(doit)
      print *, "execution returns, rc = ", shell_execute
   else
      ! make each task wait for a message from the (N-1)th task.
      if (myrank == 0) then

         ! my turn to execute
         shell_execute = system(doit)
         print *, "execution returns, rc = ", shell_execute, "iam = ", myrank

         call MPI_Send(dummy, 1, MPI_INTEGER, i+1, i+1, my_local_comm, errcode)
         if (errcode /= MPI_SUCCESS) then
            write(errstring, '(a,i8)') 'MPI_Send returned error code ', &
                                        errcode
            call error_handler(E_ERR,'send_to', errstring, source, &
                               revision, revdate)
         endif

      else if (myrank /= (total_tasks-1)) then
         ! wait for a message from task (me-1) that i can continue.
         call MPI_Recv(dummy, 1, MPI_INTEGER, myrank-1, myrank, &
                       my_local_comm, status, errcode)
         if (errcode /= MPI_SUCCESS) then
            write(errstring, '(a,i8)') 'MPI_Recv returned error code ', errcode
            call error_handler(E_ERR,'receive_from', errstring, source, &
                               revision, revdate)
         endif

         ! my turn to execute
         shell_execute = system(doit)
         print *, "execution returns, rc = ", shell_execute, "iam = ", myrank

	 ! send a message to the next task that it is good to go
         ! (unless there is no "next task").
         call MPI_Send(dummy, 1, MPI_INTEGER, i+1, i+1, my_local_comm, errcode)
         if (errcode /= MPI_SUCCESS) then
            write(errstring, '(a,i8)') 'MPI_Send returned error code ', &
                                        errcode
            call error_handler(E_ERR,'send_to', errstring, source, &
                               revision, revdate)
         endif
      else
         ! last task, no one else to send to.
         call MPI_Recv(dummy, 1, MPI_INTEGER, myrank-1, myrank, &
                       my_local_comm, status, errcode)
         if (errcode /= MPI_SUCCESS) then
            write(errstring, '(a,i8)') 'MPI_Recv returned error code ', errcode
            call error_handler(E_ERR,'receive_from', errstring, source, &
                               revision, revdate)
         endif

         ! my turn to execute
         shell_execute = system(doit)
         print *, "execution returns, rc = ", shell_execute, "iam = ", myrank

      endif
       

   endif

   ! debug only.
   !print *, "*** BYPASSING call to SYSTEM() *** "
   !shell_execute = 0


end function shell_execute

!-----------------------------------------------------------------------------
subroutine sleep_seconds(naplength)
 real(r8), intent(in) :: naplength

! Wrapper for the sleep command.  Argument is a real
! in seconds.  Different systems have different lower
! resolutions for the minimum time it will sleep.
! Subroutine, no return value.

 integer :: sleeptime

 sleeptime = floor(naplength)
 if (sleeptime <= 0) sleeptime = 1

 call sleep(sleeptime)

end subroutine sleep_seconds

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

end module mpi_utilities_mod

