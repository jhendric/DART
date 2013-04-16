!----------------------------------------------------
! Timing of transpose
!
!> @author
!> Helen 
!
!> DESCRIPTION
!> Aim to look at the scaling of the transpose 
!> all_vars_to_all_copies and all_copies_all_vars
!> my_task_id() 0 timings output to standard out
!> All my_task_id()s timings output to standard error
!> can use cpu time or MPI_WTIME (namelist flag)
!> can switch on loop timingsi (MPI_WTIME) with printme=1
!-----------------------------------------------------      
program timing
      use mpi_utilities_mod,    only : initialize_mpi_utilities,finalize_mpi_utilities,                &
                                       my_task_id, task_sync
      use utilities_mod,        only: find_namelist_in_file, check_namelist_read,                      &
                                      set_nml_output
      use ensemble_manager_mod, only: init_ensemble_manager, ensemble_type, all_vars_to_all_copies,    &
                                       all_copies_to_all_vars, end_ensemble_manager, get_my_num_vars, get_my_num_copies

       use mpi

implicit none

type(ensemble_type) :: ens_handle, obs_ens_handle, forward_op_ens_handle !> ensemble type

!namelist with default values
integer :: debug_flag = 0 !> default is 0
integer :: timing_ens_size = 80 !> ensemble size 
integer :: ens_size_extras = 6 !> number of extra columns to add on to ens_size
integer :: n_repeats = 1 !> number of repetitions of loop for a given array size
integer :: min_array_length = 20 !> minimum number of vars
integer :: inc_array_length = 10 !> step size for number of vars loop
integer :: max_array_length = 50 !> max number of vars
logical :: use_cpu_time = .false. !> default is to use MPI_WTIME

!timing variables
integer :: array_length !> length of array to transpose
integer :: i !> for repeats loop
real:: cpu_start, cpu_finish !> for cpu_time
double precision:: start, finish !> for MPI_Wtime
integer :: printme = 0!> for loop timing inside transpose routines, 0 is off, 1 is on

! utilites
integer :: iunit !> for reading in namelist from flle
integer :: async = 0 !> what does this do for mpi_finalize
integer :: io !> for iostat
integer :: stderr = 0 !> for output to standard error
double precision, allocatable :: memory_splat(:,:) 
integer :: ierr !> for MPI_BCAST

! namelist : debug flag,  ensemble size (80), number of extras (4 or 6 normally), array length( min, inc, max), number of repitions
! use_cpu_time (true/false) 
namelist /timing_nml/ debug_flag, &
 timing_ens_size, ens_size_extras, &
 min_array_length, inc_array_length, max_array_length, n_repeats, &
 use_cpu_time

! initialize mpi
 call initialize_mpi_utilities('timing')

 call find_namelist_in_file("input.nml", "timing_nml", iunit)
 read(iunit, nml = timing_nml, iostat = io)

 call check_namelist_read(iunit, io, "timing_nml")

! print run info
 if (my_task_id() == 0) print*, 'ens_size', timing_ens_size, '+ extras', ens_size_extras, 'use_cpu_time', use_cpu_time

! print output format to screen
 if (my_task_id() == 0) then
   print*, 'Timing output: message  array_length  repetition#: time'
   write(stderr, *) 'Timing output: my_task_id() message array_length repetition# : time'
 end if
 
       numvar : do array_length = min_array_length, max_array_length, inc_array_length

            numRep : do i=1, n_repeats

                ! HK loop timing for 
                 if (i==1) printme =0
                 if (i==2) printme =0

             call task_sync()

            ! setup problem: initialize the ensemble manager storage - Need enough copies for ensemble plus extras
            call init_ensemble_manager(ens_handle, timing_ens_size + ens_size_extras, array_length, 1)    
                if (my_task_id() == 0 ) then
                    if (debug_flag ==1 ) print*, 'inialized ensemble', array_length, i
                end if     

            ! fill up ensemble
            ens_handle%my_copies(:) = 42
           !ens_handles%time(:) ! what to fill this with
            ens_handle%my_vars(:) = 43
            ens_handle%vars(:,:) = 44.d0
            ens_handle%copies(:,:) = 45.d0
  
            if (i==1) then
               if( my_task_id()==1) print *, 'task', my_task_id(), 'num copies = ',  get_my_num_copies(ens_handle), 'my num vars = ', get_my_num_vars(ens_handle)
           endif

            call task_sync()

          if (use_cpu_time .eqv. .true.) then
   
         ! call transpose all_vars_to_all_copies
                    call cpu_time(start)
                    call all_vars_to_all_copies(ens_handle,printme)
                    call cpu_time(finish)
                    if (my_task_id() == 0 ) print *, 'c all_vars_to_all_copies ', array_length, i,' : ', finish - start
                   write(stderr, *)my_task_id(),' vars_to_copies',  array_length, i,' : ', finish - start

                call task_sync()
  
        ! call transpose all_copies_to_all_vars
 
                    call cpu_time(start)
                    call all_copies_to_all_vars(ens_handle, printme)
                    call cpu_time(finish)
                    if (my_task_id() == 0 ) print *, 'c all_copies_to_all_vars ', array_length, i,' : ', finish - start
                    write(stderr, *)my_task_id(),' copies_to_vars',  array_length, i,' : ', finish - start
  
              else   

          ! call transpose all_vars_to_all_copies
                    start = MPI_WTIME() 
                   call all_vars_to_all_copies(ens_handle,printme)
                    finish = MPI_WTIME() 
                    if (my_task_id() == 0 ) print *, 'c all_vars_to_all_copies ', array_length, i,' : ', finish - start
                    write(stderr, *)my_task_id(),' vars_to_copies',  array_length, i,' : ', finish - start

                call task_sync()

         ! call transpose all_copies_to_all_vars
                   start = MPI_WTIME()  
                   call all_copies_to_all_vars(ens_handle, printme)
                   finish = MPI_WTIME() 
                   if (my_task_id() == 0 ) print *, 'c all_copies_to_all_vars ', array_length, i,' : ', finish - start
                   write(stderr, *)my_task_id(),' copies_to_vars',  array_length, i,' : ', finish - start

           endif

           call task_sync()

           ! destroy storage 
          call end_ensemble_manager(ens_handle)
         if (my_task_id() == 0 ) then
          if (debug_flag ==1 ) print*, 'destroyed ensemble_storage ', array_length, i
         end if

        end do numRep

    end do numvar

! finialize mpi
 call finalize_mpi_utilities(async=async)

end program timing
