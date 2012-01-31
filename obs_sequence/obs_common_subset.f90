! DART software - Copyright 2004 - 2011 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

program obs_common_subset

! this program expects to be given pairs, or lists of pairs, of obs_seq
! files as input and does an obs-by-obs comparison.  any obs which does not
! match in type, time, location, or (if present) DART QC is discarded.
! in addition to matching, the DART QCs have to be below a given threshold.
! the intent is to run two experiments with identical input observations
! but some difference in namelist settings or other filter differences.
! then take the two output obs_seq.final files and remove any obs which
! were not assimilated in both experiments for any reason.  then comparing
! the output obs_seq files will be comparing exactly the same obs, and only
! those which were assimilated.  differences in the output diagnostics will
! be because of the experiment differences, not differences in the number of
! obs assimilated.

! running this program creates 2 new output files for each pair of input files;
! the names of the input obs_seq files with a suffix appended.  it can take
! list of obs_seq files, either explicitly listed or in a separate file with
! one input filename per line.


use        types_mod, only : r8, missing_r8, metadatalength, obstypelength
use    utilities_mod, only : register_module, initialize_utilities,            &
                             find_namelist_in_file, check_namelist_read,       &
                             error_handler, E_ERR, E_MSG, nmlfileunit,         &
                             do_nml_file, do_nml_term, get_next_filename,      &
                             finalize_utilities, logfileunit
use     location_mod, only : location_type, get_location, write_location,      &
                             operator(/=)
use      obs_def_mod, only : obs_def_type, get_obs_def_time, get_obs_kind,     &
                             get_obs_def_location, read_obs_def
use     obs_kind_mod, only : max_obs_kinds, get_obs_kind_name
use time_manager_mod, only : time_type, print_date, print_time, set_time,      &
                             set_calendar_type, get_calendar_type,             &
                             operator(/=), operator(>), NO_CALENDAR
use obs_sequence_mod, only : obs_sequence_type, obs_type, write_obs_seq,       &
                             init_obs, assignment(=), get_obs_def, get_obs_key,&
                             init_obs_sequence, static_init_obs_sequence,      &
                             read_obs_seq_header, read_obs_seq, get_num_obs,   &
                             get_first_obs, get_next_obs, get_num_key_range,   &
                             insert_obs_in_seq, get_num_copies, get_num_qc,    &
                             get_copy_meta_data, get_qc_meta_data, get_qc,     &
                             set_copy_meta_data, set_qc_meta_data,             &
                             destroy_obs, destroy_obs_sequence

implicit none

! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = "$URL$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

type(obs_sequence_type) :: seq_in1, seq_in2, seq_out1, seq_out2
type(obs_sequence_type) :: seq_in3, seq_in4, seq_out3, seq_out4
type(obs_type)          :: obs_in1, next_obs_in1, obs_in2, next_obs_in2
type(obs_type)          :: obs_in3, next_obs_in3, obs_in4, next_obs_in4
type(obs_type)          :: obs_out1, prev_obs_out1, obs_out2, prev_obs_out2
type(obs_type)          :: obs_out3, prev_obs_out3, obs_out4, prev_obs_out4
logical                 :: is_this_last1, is_this_last2
logical                 :: is_this_last3, is_this_last4
logical                 :: comparing3 = .true.
logical                 :: comparing4 = .true.
logical                 :: wanted     = .false.
integer                 :: size_seq_in1, num_copies_in1, num_qc_in1
integer                 :: size_seq_in2, num_copies_in2, num_qc_in2
integer                 :: size_seq_in3, num_copies_in3, num_qc_in3
integer                 :: size_seq_in4, num_copies_in4, num_qc_in4
integer                 :: num_copies_in, num_qc_in, size_seq_out
integer                 :: num_inserted, iunit, io, i, j
integer                 :: max_num_obs1, max_num_obs2, file_id
integer                 :: max_num_obs3, max_num_obs4
integer                 :: num_rejected_badqc, num_rejected_diffqc
integer                 :: num_rejected_other
character(len = 129)    :: read_format
logical                 :: pre_I_format, cal
character(len = 256)    :: msgstring, msgstring1, msgstring2
character(len = 164)    :: filename_out1, filename_out2
character(len = 164)    :: filename_out3, filename_out4

character(len = metadatalength) :: dart_qc_meta_data = 'DART quality control'
character(len = metadatalength) :: meta_data
integer                 :: qc_index
integer                 :: num_input_files = 0

! could go into namelist if you wanted more control
integer, parameter      :: print_every = 200
integer, parameter      :: qc_threshold = 1   ! ok if <= this

!----------------------------------------------------------------
! Namelist input with default values

! max_num_input_files : maximum number of input sequence files to be processed
! lazy, pick big numbers.  make them bigger if too small.
integer, parameter               :: max_num_input_files = 1000
integer, parameter               :: max_obs_input_types = 500

character(len = 129) :: filename_seq1(max_num_input_files) = ''
character(len = 129) :: filename_seq2(max_num_input_files) = ''
character(len = 129) :: filename_seq3(max_num_input_files) = ''
character(len = 129) :: filename_seq4(max_num_input_files) = ''
character(len = 129) :: filename_seq_list1  = ''
character(len = 129) :: filename_seq_list2  = ''
character(len = 129) :: filename_seq_list3  = ''
character(len = 129) :: filename_seq_list4  = ''
character(len = 32)  :: filename_out_suffix = '.common'

logical              :: print_only    = .false.
character(len=32)    :: calendar      = 'Gregorian'

namelist /obs_common_subset_nml/ &
         filename_seq1, filename_seq_list1, &
         filename_seq2, filename_seq_list2, &
         filename_seq3, filename_seq_list3, &
         filename_seq4, filename_seq_list4, &
         filename_out_suffix, print_only, calendar

!----------------------------------------------------------------
! Start of the program:
!
! Process each input observation sequence file in turn, optionally
! selecting observations to insert into the output sequence file.
!----------------------------------------------------------------

call setup()

! Read the namelist entry
call find_namelist_in_file("input.nml", "obs_common_subset_nml", iunit)
read(iunit, nml = obs_common_subset_nml, iostat = io)
call check_namelist_read(iunit, io, "obs_common_subset_nml")

! Record the namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=obs_common_subset_nml)
if (do_nml_term()) write(     *     , nml=obs_common_subset_nml)

! the logic here is slightly different than the obs_sequence_tool.
! the user gives us 4 lists of obs_seq files; either in the namelist
! or as the name of a file which contains the list, one per line.
! either way, the lists must be the same length.
! as in the other tools, it is an error to specify both an explicit
! list and the name of file for input.

call handle_filenames(filename_seq1, filename_seq_list1, &
                      filename_seq2, filename_seq_list2, &
                      filename_seq3, filename_seq_list3, &
                      filename_seq4, filename_seq_list4, &
                      num_input_files)

! the default is a gregorian calendar.  if you are using a different type
! set it in the namelist.  this only controls how it prints out the first
! and last timestamps in the obs_seq files.
call set_calendar_type(calendar)

! set a logial to see if we have a calendar or not
cal = (get_calendar_type() /= NO_CALENDAR)

! end of namelist processing and setup


! single pass algorithm (unlike other obs tools).  process the files in
! pairs, making sure the metadata matches.  the number of obs does NOT
! have to match.  iterate them pairwise and keep any matching obs which
! share the same QC values, as long as the QCs are <= threshold.
! (QC tests are skipped if DART QC isn't in the metadata.)

! count of input files was set in the handle_filenames() routine above.
do i = 1, num_input_files

   if  ((len(filename_seq1(i)) == 0) .or. (filename_seq1(i) == "")) then
      write(msgstring, *) 'entry ', i, 'for sequence1 is empty or null'
      call error_handler(E_ERR,'obs_common_subset', msgstring, &
         source,revision,revdate)  ! shouldn't happen
   endif
   if  ((len(filename_seq2(i)) == 0) .or. (filename_seq2(i) == "")) then
      write(msgstring, *) 'entry ', i, 'for sequence2 is empty or null'
      call error_handler(E_ERR,'obs_common_subset', msgstring, &
         source,revision,revdate)  ! shouldn't happen
   endif

   if  ((len(filename_seq3(i)) == 0) .or. (filename_seq3(i) == "")) then
      comparing3 = .false.
      write(msgstring, *) 'entry ', i, 'for sequence3 is empty or null'
      call error_handler(E_MSG,'obs_common_subset', msgstring, &
         source,revision,revdate)  ! may happen
   endif

   if  ((len(filename_seq4(i)) == 0) .or. (filename_seq4(i) == "")) then
      comparing4 = .false.
      write(msgstring, *) 'entry ', i, 'for sequence4 is empty or null'
      call error_handler(E_MSG,'obs_common_subset', msgstring, &
         source,revision,revdate)  ! may happen
   endif

   ! read in the next pair of files.

   call read_obs_seq_header(filename_seq1(i), num_copies_in1, num_qc_in1, &
      size_seq_in1, max_num_obs1, file_id, read_format, pre_I_format, &
      close_the_file = .true.)

   if (max_num_obs1 == 0) then
      write(msgstring,*) 'No obs in input sequence file ', trim(filename_seq1(i))
      call error_handler(E_MSG,'obs_common_subset',msgstring)
      cycle
   endif

   call read_obs_seq_header(filename_seq2(i), num_copies_in2, num_qc_in2, &
      size_seq_in2, max_num_obs2, file_id, read_format, pre_I_format, &
      close_the_file = .true.)

   if (max_num_obs2 == 0) then
      write(msgstring,*) 'No obs in input sequence file ', trim(filename_seq2(i))
      call error_handler(E_MSG,'obs_common_subset',msgstring)
      cycle
   endif

   if (comparing3) then
      call read_obs_seq_header(filename_seq3(i), num_copies_in3, num_qc_in3, &
         size_seq_in3, max_num_obs3, file_id, read_format, pre_I_format, &
         close_the_file = .true.)

      if (max_num_obs3 == 0) then
         write(msgstring,*) 'No obs in input sequence file ', trim(filename_seq3(i))
         call error_handler(E_MSG,'obs_common_subset',msgstring)
         cycle
      endif
   endif

   if (comparing4) then
      call read_obs_seq_header(filename_seq4(i), num_copies_in4, num_qc_in4, &
         size_seq_in4, max_num_obs4, file_id, read_format, pre_I_format, &
         close_the_file = .true.)

      if (max_num_obs4 == 0) then
         write(msgstring,*) 'No obs in input sequence file ', trim(filename_seq4(i))
         call error_handler(E_MSG,'obs_common_subset',msgstring)
         cycle
      endif
   endif

   write(*, *) 'Starting to process input sequence files: '
                   write(*,*)  trim(filename_seq1(i))
                   write(*,*)  trim(filename_seq2(i))
   if (comparing3) write(*,*)  trim(filename_seq3(i))
   if (comparing4) write(*,*)  trim(filename_seq4(i))

   write(logfileunit, *) 'Starting to process input sequence files: '
                   write(logfileunit,*)  trim(filename_seq1(i))
                   write(logfileunit,*)  trim(filename_seq2(i))
   if (comparing3) write(logfileunit,*)  trim(filename_seq3(i))
   if (comparing4) write(logfileunit,*)  trim(filename_seq4(i))

                   call read_obs_seq(filename_seq1(i), 0, 0, 0, seq_in1)
                   call read_obs_seq(filename_seq2(i), 0, 0, 0, seq_in2)
   if (comparing3) call read_obs_seq(filename_seq3(i), 0, 0, 0, seq_in3)
   if (comparing4) call read_obs_seq(filename_seq4(i), 0, 0, 0, seq_in4)

   ! make sure the files have compatible metadata.  this errors out here if not.
   call compare_metadata(seq_in1, seq_in2, filename_seq1(i), filename_seq2(i))

   ! if we get here, the number of copies and qcs must match.  use a single
   ! variable for both.
   num_copies_in = num_copies_in1
   num_qc_in     = num_qc_in1

   ! sanity check - ensure the linked list times are in increasing time order
                   call validate_obs_seq_time(seq_in1, filename_seq1(i))
                   call validate_obs_seq_time(seq_in2, filename_seq2(i))
   if (comparing3) call validate_obs_seq_time(seq_in3, filename_seq3(i))
   if (comparing4) call validate_obs_seq_time(seq_in4, filename_seq4(i))

   ! the output can have no more than the minimum number of input obs.
   size_seq_out = min(size_seq_in1, size_seq_in2)
   if (comparing3) size_seq_out = min(size_seq_out, size_seq_in3)
   if (comparing4) size_seq_out = min(size_seq_out, size_seq_in4)

   ! find which DART qc copy has the DART quality control.  set to -1 if not
   ! present, which will skip the threshold test.   and at this point we know
   ! the two sequences are the same so just pick one and query it.
   qc_index = -1
   do j=1, get_num_qc(seq_in1)
      if(index(get_qc_meta_data(seq_in1,j), dart_qc_meta_data) > 0) then
         qc_index = j
         exit
      endif
   enddo

   ! blank line, start of actually creating output file
   call error_handler(E_MSG,' ',' ')

   ! Initialize individual observation variables
   call init_obs(      obs_in1, num_copies_in, num_qc_in)
   call init_obs( next_obs_in1, num_copies_in, num_qc_in)
   call init_obs(     obs_out1, num_copies_in, num_qc_in)
   call init_obs(prev_obs_out1, num_copies_in, num_qc_in)

   call init_obs(      obs_in2, num_copies_in, num_qc_in)
   call init_obs( next_obs_in2, num_copies_in, num_qc_in)
   call init_obs(     obs_out2, num_copies_in, num_qc_in)
   call init_obs(prev_obs_out2, num_copies_in, num_qc_in)

   if (comparing3) then
      call init_obs(      obs_in3, num_copies_in, num_qc_in)
      call init_obs( next_obs_in3, num_copies_in, num_qc_in)
      call init_obs(     obs_out3, num_copies_in, num_qc_in)
      call init_obs(prev_obs_out3, num_copies_in, num_qc_in)
   endif

   if (comparing4) then
      call init_obs(      obs_in4, num_copies_in, num_qc_in)
      call init_obs( next_obs_in4, num_copies_in, num_qc_in)
      call init_obs(     obs_out4, num_copies_in, num_qc_in)
      call init_obs(prev_obs_out4, num_copies_in, num_qc_in)
   endif

   ! create the output sequences here
                   call init_obs_sequence(seq_out1, num_copies_in, num_qc_in, size_seq_out)
                   call init_obs_sequence(seq_out2, num_copies_in, num_qc_in, size_seq_out)
   if (comparing3) call init_obs_sequence(seq_out3, num_copies_in, num_qc_in, size_seq_out)
   if (comparing4) call init_obs_sequence(seq_out4, num_copies_in, num_qc_in, size_seq_out)

   do j=1, num_copies_in
      meta_data = get_copy_meta_data(seq_in1, j)
                      call set_copy_meta_data(seq_out1, j, meta_data)
                      call set_copy_meta_data(seq_out2, j, meta_data)
      if (comparing3) call set_copy_meta_data(seq_out3, j, meta_data)
      if (comparing4) call set_copy_meta_data(seq_out4, j, meta_data)
   enddo
   do j=1, num_qc_in
      meta_data = get_qc_meta_data(seq_in1, j)
                      call set_qc_meta_data(seq_out1, j, meta_data)
                      call set_qc_meta_data(seq_out2, j, meta_data)
      if (comparing3) call set_qc_meta_data(seq_out3, j, meta_data)
      if (comparing4) call set_qc_meta_data(seq_out4, j, meta_data)
   enddo

   ! not sure what this should do - print after it has generated the common set?
   if (print_only)                  call print_obs_seq(seq_in1, filename_seq1(i))
   if (print_only)                  call print_obs_seq(seq_in2, filename_seq2(i))
   if (print_only .and. comparing3) call print_obs_seq(seq_in3, filename_seq3(i))
   if (print_only .and. comparing4) call print_obs_seq(seq_in4, filename_seq4(i))

   !-------------------------------------------------------------
   ! Start to insert obs from sequence_in into sequence_out
   !
   ! NOTE: insert_obs_in_seq CHANGES the obs passed in.
   !       Must pass a copy of incoming obs to insert_obs_in_seq.
   !--------------------------------------------------------------
   num_inserted = 0
   num_rejected_badqc = 0
   num_rejected_diffqc = 0
   num_rejected_other = 0

   if ( get_first_obs(seq_in1, obs_in1) .and. get_first_obs(seq_in2, obs_in2) )  then

      if (comparing3) then
         if (.not. get_first_obs(seq_in3, obs_in3)) then
            write(msgstring, *)'no first observation in ',trim(filename_seq3(i))
            call error_handler(E_ERR,'obs_common_subset',msgstring,source,revision,revdate)
         endif
      endif

      if (comparing4) then
         if (.not. get_first_obs(seq_in4, obs_in4)) then
            write(msgstring, *)'no first observation in ',trim(filename_seq4(i))
            call error_handler(E_ERR,'obs_common_subset',msgstring,source,revision,revdate)
         endif
      endif

      is_this_last1 = .false.
      is_this_last2 = .false.
      is_this_last3 = .false.
      is_this_last4 = .false.
                      next_obs_in1 = obs_in1
                      next_obs_in2 = obs_in2
      if (comparing3) next_obs_in3 = obs_in3
      if (comparing4) next_obs_in4 = obs_in4

      ObsLoop : do while ( .not. is_this_last1 .and. .not. is_this_last2)

         ! essentially record position in seq_out
                         obs_in1 = next_obs_in1
                         obs_in2 = next_obs_in2
         if (comparing3) obs_in3 = next_obs_in3
         if (comparing4) obs_in4 = next_obs_in4

         if (comparing4) then
            wanted = all_good(obs_in1, obs_in2, qc_index, qc_threshold, obs_in3, obs_in4)
         elseif (comparing3) then
            wanted = all_good(obs_in1, obs_in2, qc_index, qc_threshold, obs_in3)
         else
            wanted = all_good(obs_in1, obs_in2, qc_index, qc_threshold)
         endif

         if ( wanted ) then
                            obs_out1 = obs_in1
                            obs_out2 = obs_in2
            if (comparing3) obs_out3 = obs_in3
            if (comparing4) obs_out4 = obs_in4

            ! Since the stride through the observation sequence file is always
            ! guaranteed to be in temporally-ascending order, we can use the
            ! 'previous' observation as the starting point to search for the
            ! correct insertion point.  This speeds up the insert code a lot.

            if (num_inserted > 0) then
                               call insert_obs_in_seq(seq_out1, obs_out1, prev_obs_out1)
                               call insert_obs_in_seq(seq_out2, obs_out2, prev_obs_out2)
               if (comparing3) call insert_obs_in_seq(seq_out3, obs_out3, prev_obs_out3)
               if (comparing4) call insert_obs_in_seq(seq_out4, obs_out4, prev_obs_out4)
            else
                               call insert_obs_in_seq(seq_out1, obs_out1)
                               call insert_obs_in_seq(seq_out2, obs_out2)
               if (comparing3) call insert_obs_in_seq(seq_out3, obs_out3)
               if (comparing4) call insert_obs_in_seq(seq_out4, obs_out4)
            endif

            ! update position in seq for next insert
                            prev_obs_out1 = obs_out1
                            prev_obs_out2 = obs_out2
            if (comparing3) prev_obs_out3 = obs_out3
            if (comparing4) prev_obs_out4 = obs_out4
            num_inserted = num_inserted + 1

            if (print_every > 0) then
               if (mod(num_inserted,print_every) == 0) then
                  print*, 'inserted number ',num_inserted,' of ',size_seq_out
               endif
            endif

         endif

                         call get_next_obs(seq_in1, obs_in1, next_obs_in1, is_this_last1)
                         call get_next_obs(seq_in2, obs_in2, next_obs_in2, is_this_last2)
         if (comparing3) call get_next_obs(seq_in3, obs_in3, next_obs_in3, is_this_last3)
         if (comparing4) call get_next_obs(seq_in4, obs_in4, next_obs_in4, is_this_last4)

      enddo ObsLoop

   else
      write(msgstring, *)'no first observation in ',trim(filename_seq1(i))
      write(msgstring1,*)' or in ',trim(filename_seq2(i))
      call error_handler(E_MSG,'obs_common_subset', msgstring, text2=msgstring1)
   endif

   if (.not. print_only) then
      print*, '---------  Obs seq file pair #  :         ', i
      print*, 'Number of obs in sequence 1     :         ', size_seq_in1
      print*, 'Number of obs in sequence 2     :         ', size_seq_in2
      print*, 'Number of obs rejected diff qc  :         ', num_rejected_diffqc
      print*, 'Number of obs rejected bad qc   :         ', num_rejected_badqc
      print*, 'Number of obs rejected other    :         ', num_rejected_other
      print*, 'Number of obs copied to output  :         ', num_inserted
      print*, '---------------------------------------------------------'
   endif

                   call destroy_obs_sequence(seq_in1)
                   call destroy_obs_sequence(seq_in2)
   if (comparing3) call destroy_obs_sequence(seq_in3)
   if (comparing4) call destroy_obs_sequence(seq_in4)

                   filename_out1 = trim(filename_seq1(i))//trim(filename_out_suffix)
                   filename_out2 = trim(filename_seq2(i))//trim(filename_out_suffix)
   if (comparing3) filename_out3 = trim(filename_seq3(i))//trim(filename_out_suffix)
   if (comparing4) filename_out4 = trim(filename_seq4(i))//trim(filename_out_suffix)

   write(msgstring, *) 'Starting to write output sequence files'
   call error_handler(E_MSG,'obs_common_subset',msgstring)

                   print*, 'Number of obs in the output seq file :', get_num_key_range(seq_out1)
                   print*, 'and                                  :', get_num_key_range(seq_out2)
   if (comparing3) print*, 'and                                  :', get_num_key_range(seq_out3)
   if (comparing4) print*, 'and                                  :', get_num_key_range(seq_out4)

                   call print_obs_seq(seq_out1, filename_out1)
                   call print_obs_seq(seq_out2, filename_out2)
   if (comparing3) call print_obs_seq(seq_out3, filename_out3)
   if (comparing4) call print_obs_seq(seq_out4, filename_out4)
   if (.not. print_only) then
                      call write_obs_seq(seq_out1, filename_out1)
                      call write_obs_seq(seq_out2, filename_out2)
      if (comparing3) call write_obs_seq(seq_out3, filename_out3)
      if (comparing4) call write_obs_seq(seq_out4, filename_out4)
   else
      write(msgstring,*) 'Output sequence files not created; print_only in namelist is .true.'
      call error_handler(E_MSG,'', msgstring)
   endif

   ! clean up

   call destroy_obs_sequence(seq_out1)
   call destroy_obs(         obs_out1)
   call destroy_obs(          obs_in1)
   call destroy_obs(     next_obs_in1)

   call destroy_obs_sequence(seq_out2)
   call destroy_obs(         obs_out2)
   call destroy_obs(          obs_in2)
   call destroy_obs(     next_obs_in2)

   if (comparing3) then
      call destroy_obs_sequence(seq_out3)
      call destroy_obs(         obs_out3)
      call destroy_obs(          obs_in3)
      call destroy_obs(     next_obs_in3)
   endif

   if (comparing4) then
      call destroy_obs_sequence(seq_out4)
      call destroy_obs(         obs_out4)
      call destroy_obs(          obs_in4)
      call destroy_obs(     next_obs_in4)
   endif

   !call destroy_obs(prev_obs_out1)  ! these are copies of something already
   !call destroy_obs(prev_obs_out2)  ! deleted; duplicate dels if called.

enddo

call shutdown()

!---------------------------------------------------------------------
! end of main program.
!---------------------------------------------------------------------


contains


!---------------------------------------------------------------------
subroutine setup()

! Initialize modules used that require it
call initialize_utilities('obs_common_subset')
call register_module(source,revision,revdate)
call static_init_obs_sequence()

end subroutine setup

!---------------------------------------------------------------------
subroutine shutdown()

call finalize_utilities('obs_common_subset')

end subroutine shutdown

!---------------------------------------------------------------------
subroutine handle_filenames(filename_seq1, filename_seq_list1, &
                            filename_seq2, filename_seq_list2, &
                            filename_seq3, filename_seq_list3, &
                            filename_seq4, filename_seq_list4, &
                            num_input_files)

! sort out the input lists, set the length as a return in num_input_files,
! and make sure what's specified is consistent.

character(len=*), intent(inout) :: filename_seq1(:), filename_seq2(:)
character(len=*), intent(inout) :: filename_seq3(:), filename_seq4(:)
character(len=*), intent(in)    :: filename_seq_list1, filename_seq_list2
character(len=*), intent(in)    :: filename_seq_list3, filename_seq_list4
integer,          intent(out)   :: num_input_files

integer :: indx
logical :: from_file1, from_file2, from_file3, from_file4
character(len=32) :: source1, source2, source3, source4

! if the user specifies neither filename_seq nor filename_seq_list, error
! if the user specifies both, error.
! if the two list lengths aren't equal, error.
! if the user gives one or more filelist names, we make sure the length(s) are
!  not more than maxfiles and read it/them into the explicit list and continue.
! set num_input_files to the number of pairs in the lists

if (filename_seq1(1) == '' .and. filename_seq_list1 == '') then
   call error_handler(E_ERR,'handle_filenames',            &
                      'no filenames specified as input 1',  &
                      source,revision,revdate)
endif
if (filename_seq2(1) == '' .and. filename_seq_list2 == '') then
   call error_handler(E_ERR,'handle_filenames',            &
                      'no filenames specified as input 2',  &
                      source,revision,revdate)
endif
if (filename_seq3(1) == '' .and. filename_seq_list3 == '') then
   comparing3 = .false.
endif
if (filename_seq4(1) == '' .and. filename_seq_list4 == '') then
   comparing4 = .false.
endif

if (comparing4 .and. .not. comparing3) then
   call error_handler(E_ERR,'handle_filenames', &
       'cannot specify filename4 without specifying filename3', &
       source,revision,revdate, &
       text2='it just makes the logic too convoluted.')
endif

! make sure the namelist specifies one or the other but not both
if ((filename_seq1(1) /= '' .and. filename_seq_list1 /= '') .or. &
    (filename_seq2(1) /= '' .and. filename_seq_list2 /= '') .or. &
    (filename_seq3(1) /= '' .and. filename_seq_list3 /= '') .or. &
    (filename_seq4(1) /= '' .and. filename_seq_list4 /= '')) then
   call error_handler(E_ERR,'handle_filenames', &
       'cannot specify both filename_seq and filename_seq_list', &
       source,revision,revdate)
endif

! if they have specified a file which contains a list, read it into
! the filename_seq array and set the count.
if (filename_seq_list1 /= '') then
   source1 = 'filename_seq_list1'
   from_file1 = .true.
else
   source1 = 'filename_seq1'
   from_file1 = .false.
endif

if (filename_seq_list2 /= '') then
   source2 = 'filename_seq_list2'
   from_file2 = .true.
else
   source2 = 'filename_seq2'
   from_file2 = .false.
endif

if ((filename_seq_list3 /= '') .and. comparing3) then
   source3 = 'filename_seq_list3'
   from_file3 = .true.
else
   source3 = 'filename_seq3'
   from_file3 = .false.
endif

if ((filename_seq_list4 /= '') .and. comparing4) then
   source4 = 'filename_seq_list4'
   from_file4 = .true.
else
   source4 = 'filename_seq4'
   from_file4 = .false.
endif

                write(*,*)'filename_seq1(1) ',trim(filename_seq1(1)),' ', trim(source1), ' ',from_file1
                write(*,*)'filename_seq2(1) ',trim(filename_seq2(1)),' ', trim(source2), ' ',from_file2
if (comparing3) write(*,*)'filename_seq3(1) ',trim(filename_seq3(1)),' ', trim(source3), ' ',from_file3
if (comparing4) write(*,*)'filename_seq4(1) ',trim(filename_seq4(1)),' ', trim(source4), ' ',from_file4

! the point of this loop is to count up how many pairs of input seq files we have.
do indx = 1, max_num_input_files
   if (from_file1) &
      filename_seq1(indx) = get_next_filename(filename_seq_list1, indx)
   if (from_file2) &
      filename_seq2(indx) = get_next_filename(filename_seq_list2, indx)
   if (from_file3 .and. comparing3 ) &
      filename_seq3(indx) = get_next_filename(filename_seq_list3, indx)
   if (from_file4 .and. comparing4 ) &
      filename_seq4(indx) = get_next_filename(filename_seq_list4, indx)

   ! a pair of empty names ends the list and we return with the count.
   ! (unless both lists are empty and then we're unhappy)
   if ((filename_seq1(indx) == '') .and. (filename_seq2(indx) == '')) then
      if (indx == 1) then
         call error_handler(E_ERR,'handle_filenames', &
             trim(source1)//' contains no filenames', &
             source,revision,revdate)
      endif

      ! if this isn't the first entry, we are at the end of the list.
      ! set the return count before going.
      num_input_files = indx - 1
      return

   ! catch the cases where the lists aren't the same length, 2 longer than 1
   else if (filename_seq1(indx) == '') then
         call error_handler(E_ERR,'handle_filenames', &
             trim(source2)//' contains more filenames than '//trim(source1), &
             source,revision,revdate)

   ! catch the other case where the lists aren't the same length, 1 longer than 2
   else if (filename_seq2(indx) == '') then
         call error_handler(E_ERR,'handle_filenames', &
             trim(source1)//' contains more filenames than '//trim(source2), &
             source,revision,revdate)
   endif

   ! if both lists had filenames, we didn't match any of the if tests and
   ! we just loop, counting up the number of names.

enddo

write(msgstring, *) 'cannot specify more than ',max_num_input_files,' files'
call error_handler(E_ERR,'handle_filenames', msgstring, source,revision,revdate)

end subroutine handle_filenames

!---------------------------------------------------------------------
subroutine compare_metadata(seq1, seq2, fname1, fname2)

!
! This subroutine compares the metadata for two different observation
! sequences and terminates the program if they are not conformable.
! In order to be merged, the two observation sequences must have the same
! number of qc values, the same number of copies ...
!

 type(obs_sequence_type), intent(IN) :: seq1, seq2
 character(len=*), optional :: fname1, fname2

integer :: num_copies1, num_qc1
integer :: num_copies2, num_qc2
integer :: num_copies , num_qc, i
character(len=metadatalength) :: str1, str2
character(len=255) :: msgstring3

num_copies1 = get_num_copies(seq1)
num_qc1     = get_num_qc(    seq1)

num_copies2 = get_num_copies(seq2)
num_qc2     = get_num_qc(    seq2)

num_copies  = num_copies2
num_qc      = num_qc2

! get this ready in case we have to use it in some later error string.
if (present(fname1) .and. present(fname2)) then
   write(msgstring3,*)'Sequence files ', trim(fname1), ' and ', trim(fname2), &
                      ' are not compatible'
else
  msgstring3 = 'Sequence files cannot be processed because they are not compatible'
endif

if ( num_copies1 /= num_copies2 ) then
   write(msgstring2,*)'Different numbers of data copies found: ', &
                      num_copies1, ' vs ', num_copies2
   num_copies = -1
endif
if ( num_qc1 /= num_qc2 ) then
   write(msgstring2,*)'Different different numbers of QCs found: ', &
                      num_qc1, ' vs ', num_qc2
   num_qc = -1
endif
if ( num_copies < 0 .or. num_qc < 0 ) then
   call error_handler(E_ERR, 'compare_metadata', msgstring3, &
                              source, revision, revdate, text2=msgstring2)
endif

! make sure the metadata matches exactly for both the copies and qcs
! in both files.
CopyMetaData : do i=1, num_copies
   str1 = get_copy_meta_data(seq1,i)
   str2 = get_copy_meta_data(seq2,i)

   ! if they match, write out an informational message and continue.
   ! if they don't match, it's a fatal error.
   if( str1 == str2 ) then
      write(msgstring,*)'copy metadata ',trim(str1), ' in both.'
      call error_handler(E_MSG, 'compare_metadata', msgstring)
   else
      write(msgstring1,*)'copy metadata value mismatch. seq1: ', trim(str1)
      write(msgstring2,*)'copy metadata value mismatch. seq2: ', trim(str2)
      call error_handler(E_ERR, 'compare_metadata', msgstring3, &
              source, revision, revdate, text2=msgstring1, text3=msgstring2)
   endif

enddo CopyMetaData

QCMetaData : do i=1, num_qc
   str1 = get_qc_meta_data(seq1,i)
   str2 = get_qc_meta_data(seq2,i)

   ! if they match, write out an informational message and continue.
   ! if they don't match, it's a fatal error.
   if( str1 == str2 ) then
      write(msgstring,*)'  qc metadata ',trim(str1), ' in both.'
      call error_handler(E_MSG, 'compare_metadata', msgstring)
   else
      write(msgstring1,*)'qc metadata value mismatch. seq1: ', trim(str1)
      write(msgstring2,*)'qc metadata value mismatch. seq2: ', trim(str2)
      call error_handler(E_ERR, 'compare_metadata', msgstring3, &
              source, revision, revdate, text2=msgstring1, text3=msgstring2)
   endif

enddo QCMetaData

end subroutine compare_metadata

!---------------------------------------------------------------------
subroutine print_obs_seq(seq_in, filename)

! you can get more info by running the obs_diag program, but this
! prints out a quick table of obs types and counts, overall start and
! stop times, and metadata strings and counts.

type(obs_sequence_type), intent(in) :: seq_in
character(len=*), intent(in)        :: filename

type(obs_type)          :: obs, next_obs
type(obs_def_type)      :: this_obs_def
logical                 :: is_there_one, is_this_last
integer                 :: size_seq_in
integer                 :: i
integer                 :: this_obs_kind
! max_obs_kinds is a public from obs_kind_mod.f90 and really is
! counting the max number of types, not kinds
integer                 :: type_count(max_obs_kinds), identity_count


! Initialize input obs_types
do i = 1, max_obs_kinds
   type_count(i) = 0
enddo
identity_count = 0

! make sure there are obs left to process before going on.
! num_obs should be ok since we just constructed this seq so it should
! have no unlinked obs.  if it might for some reason, use this instead:
! size_seq_in = get_num_key_range(seq_in)     !current size of seq_in

size_seq_in = get_num_obs(seq_in)
if (size_seq_in == 0) then
   msgstring = 'Obs_seq file '//trim(filename)//' is empty.'
   call error_handler(E_MSG,'print_obs_seq',msgstring)
   return
endif

! Initialize individual observation variables
call init_obs(     obs, get_num_copies(seq_in), get_num_qc(seq_in))
call init_obs(next_obs, get_num_copies(seq_in), get_num_qc(seq_in))

! blank line
call error_handler(E_MSG,'',' ')

write(msgstring,*) 'Processing sequence file ', trim(filename)
call error_handler(E_MSG,'',msgstring)

call print_metadata(seq_in, filename)

!-------------------------------------------------------------
! Start to process obs from seq_in
!--------------------------------------------------------------
is_there_one = get_first_obs(seq_in, obs)

if ( .not. is_there_one )  then
   write(msgstring,*)'no first observation in ',trim(filename)
   call error_handler(E_MSG,'print_obs_seq', msgstring)
endif

! process it here
is_this_last = .false.

call get_obs_def(obs, this_obs_def)
call print_time(get_obs_def_time(this_obs_def), ' First timestamp: ')
! does not work with NO_CALENDAR
if (cal) call print_date(get_obs_def_time(this_obs_def), '   calendar Date: ')

ObsLoop : do while ( .not. is_this_last)

   call get_obs_def(obs, this_obs_def)
   this_obs_kind = get_obs_kind(this_obs_def)
   if (this_obs_kind < 0) then
      identity_count = identity_count + 1
   else
      type_count(this_obs_kind) = type_count(this_obs_kind) + 1
   endif
!   print *, 'obs kind index = ', this_obs_kind
!   if(this_obs_kind > 0)print *, 'obs name = ', get_obs_kind_name(this_obs_kind)

   call get_next_obs(seq_in, obs, next_obs, is_this_last)
   if (.not. is_this_last) then
      obs = next_obs
   else
      call print_time(get_obs_def_time(this_obs_def), '  Last timestamp: ')
      if (cal) call print_date(get_obs_def_time(this_obs_def), '   calendar Date: ')
   endif

enddo ObsLoop


write(msgstring, *) 'Number of obs processed  :          ', size_seq_in
call error_handler(E_MSG, '', msgstring)
write(msgstring, *) '---------------------------------------------------------'
call error_handler(E_MSG, '', msgstring)
do i = 1, max_obs_kinds
   if (type_count(i) > 0) then
      write(msgstring, '(a32,i8,a)') trim(get_obs_kind_name(i)), &
                                     type_count(i), ' obs'
      call error_handler(E_MSG, '', msgstring)
   endif
enddo
if (identity_count > 0) then
   write(msgstring, '(a32,i8,a)') 'Identity observations', &
                                  identity_count, ' obs'
   call error_handler(E_MSG, '', msgstring)
endif

! another blank line
call error_handler(E_MSG, '', ' ')

! Time to clean up

call destroy_obs(     obs)
call destroy_obs(next_obs)

end subroutine print_obs_seq

!---------------------------------------------------------------------
subroutine validate_obs_seq_time(seq, filename)

! this eventually belongs in the obs_seq_mod code, but for now
! try it out here.  we just fixed a hole in the interactive create
! routine which would silently let you create out-of-time-order
! linked lists, which gave no errors but didn't assimilate the
! right obs at the right time when running filter.   this runs
! through the times in the entire sequence, ensuring they are
! monotonically increasing in time.  this should help catch any
! bad files which were created with older versions of code.

type(obs_sequence_type), intent(in) :: seq
character(len=*),        intent(in) :: filename

type(obs_type)          :: obs, next_obs
type(obs_def_type)      :: this_obs_def
logical                 :: is_there_one, is_this_last
integer                 :: size_seq, obs_count
integer                 :: key
type(time_type)         :: last_time, this_time


! make sure there are obs left to process before going on.
size_seq = get_num_obs(seq)
if (size_seq == 0) then
   msgstring = 'Obs_seq file '//trim(filename)//' is empty.'
   call error_handler(E_MSG,'validate_obs_seq_time',msgstring)
   return
endif

! Initialize individual observation variables
call init_obs(     obs, get_num_copies(seq), get_num_qc(seq))
call init_obs(next_obs, get_num_copies(seq), get_num_qc(seq))

obs_count = 0

!-------------------------------------------------------------
! Start to process obs from seq
!--------------------------------------------------------------
is_there_one = get_first_obs(seq, obs)

! we already tested for 0 obs above, so there should be a first obs here.
if ( .not. is_there_one )  then
   write(msgstring,*)'no first obs in sequence ' // trim(filename)
   call error_handler(E_ERR,'validate_obs_seq_time', &
                      msgstring, source, revision, revdate)
   return
endif

is_this_last = .false.
last_time = set_time(0, 0)
ObsLoop : do while ( .not. is_this_last)

   call get_obs_def(obs, this_obs_def)
   this_time = get_obs_def_time(this_obs_def)

   if (last_time > this_time) then
      ! bad time order of observations in linked list
      call print_time(last_time, ' previous timestamp: ')
      if (cal) call print_date(last_time, '      calendar date: ')
      call print_time(this_time, '     next timestamp: ')
      if (cal) call print_date(this_time, '      calendar date: ')

      key = get_obs_key(obs)
      write(msgstring1,*)'obs number ', key, ' has earlier time than previous obs'
      write(msgstring2,*)'observations must be in increasing time order, file ' // trim(filename)
      call error_handler(E_ERR,'validate_obs_seq_time', msgstring2, &
                         source, revision, revdate, &
                         text2=msgstring1)
   endif

   last_time = this_time
   obs_count = obs_count + 1

   call get_next_obs(seq, obs, next_obs, is_this_last)
   if (.not. is_this_last) obs = next_obs

enddo ObsLoop

! clean up
call destroy_obs(     obs)
call destroy_obs(next_obs)

! technically not a time validation, but easy to check.  obs_count should never
! be larger than size_seq - that's a fatal error.  obs_count < size_seq would
! suggest there are obs in the file that aren't part of the linked list.
! this does not necessarily indicate a fatal error but it's not a common
! situation and might indicate someone should check on the file.
if (obs_count /= size_seq) then
   write(msgstring,*) 'input sequence ', trim(filename)
   call error_handler(E_MSG,'validate_obs_seq_time', msgstring)

   write(msgstring,*) 'total obs in file: ', size_seq, '  obs in linked list: ', obs_count
   if (obs_count > size_seq) then
      ! this is a fatal error
      write(msgstring1,*) 'linked list obs_count > total size_seq, should not happen'
      call error_handler(E_ERR,'validate_obs_seq_time', msgstring, &
                         source, revision, revdate, &
                         text2=msgstring1)
   else
      ! just warning msg
      write(msgstring1,*) 'only observations in linked list will be processed'
      call error_handler(E_MSG,'validate_obs_seq_time', msgstring, &
                         source, revision, revdate, text2=msgstring1)
   endif
endif

end subroutine validate_obs_seq_time

!---------------------------------------------------------------------
subroutine print_metadata(seq1, fname1)

!
! print out the metadata strings, trimmed
!

type(obs_sequence_type), intent(in) :: seq1
character(len=*), optional :: fname1

integer :: num_copies , num_qc, i
character(len=metadatalength) :: str1
character(len=255) :: msgstring3

num_copies = get_num_copies(seq1)
num_qc     = get_num_qc(    seq1)

if ( num_copies < 0 .or. num_qc < 0 ) then
   write(msgstring3,*)' illegal copy or obs count in file '//trim(fname1)
   call error_handler(E_ERR, 'print_metadata', msgstring3, &
                      source, revision, revdate)
endif

MetaDataLoop : do i=1, num_copies
   str1 = get_copy_meta_data(seq1,i)

   write(msgstring,*)'Data Metadata: ',trim(str1)
   call error_handler(E_MSG, '', msgstring)

enddo MetaDataLoop

QCMetaData : do i=1, num_qc
   str1 = get_qc_meta_data(seq1,i)

   write(msgstring,*)'  QC Metadata: ', trim(str1)
   call error_handler(E_MSG, '', msgstring)

enddo QCMetaData

end subroutine print_metadata


! compare location, time, type, QC, plus make sure QC <= threshold
!---------------------------------------------------------------------
function all_good(obs_in1, obs_in2, qc_index, qc_threshold, obs_in3, obs_in4)
type(obs_type),           intent(in) :: obs_in1, obs_in2
integer,                  intent(in) :: qc_index, qc_threshold
type(obs_type), optional, intent(in) :: obs_in3, obs_in4
logical                              :: all_good

type(obs_def_type)  :: test1_obs_def,  test2_obs_def,  test3_obs_def,  test4_obs_def
integer             :: test1_obs_type, test2_obs_type, test3_obs_type, test4_obs_type
type(time_type)     :: test1_obs_time, test2_obs_time, test3_obs_time, test4_obs_time
type(location_type) :: test1_obs_loc,  test2_obs_loc,  test3_obs_loc,  test4_obs_loc
integer             :: test1_qc,       test2_qc,       test3_qc,       test4_qc
real(r8)            :: temp(1)

! assume failure so we can return as soon as we know they don't match.
all_good = .false.

call get_obs_def(obs_in1, test1_obs_def)
test1_obs_loc  = get_obs_def_location(test1_obs_def)
test1_obs_time = get_obs_def_time(test1_obs_def)
test1_obs_type = get_obs_kind(test1_obs_def)
if (qc_index < 0) then
   test1_qc = 0
else
   call get_qc(obs_in1, temp, qc_index)
   test1_qc = nint(temp(1))
endif

call get_obs_def(obs_in2, test2_obs_def)
test2_obs_loc  = get_obs_def_location(test2_obs_def)
test2_obs_time = get_obs_def_time(test2_obs_def)
test2_obs_type = get_obs_kind(test2_obs_def)
if (qc_index < 0) then
   test2_qc = 0
else
   call get_qc(obs_in2, temp, qc_index)
   test2_qc = nint(temp(1))
endif

if (present(obs_in3)) then
   call get_obs_def(obs_in3, test3_obs_def)
   test3_obs_loc  = get_obs_def_location(test3_obs_def)
   test3_obs_time = get_obs_def_time(test3_obs_def)
   test3_obs_type = get_obs_kind(test3_obs_def)
   if (qc_index < 0) then
      test3_qc = 0
   else
      call get_qc(obs_in3, temp, qc_index)
      test3_qc = nint(temp(1))
   endif
else
   test3_qc       = test1_qc
   test3_obs_type = test1_obs_type
   test3_obs_time = test1_obs_time
   test3_obs_loc  = test1_obs_loc
endif

if (present(obs_in4)) then
   call get_obs_def(obs_in4, test4_obs_def)
   test4_obs_loc  = get_obs_def_location(test4_obs_def)
   test4_obs_time = get_obs_def_time(test4_obs_def)
   test4_obs_type = get_obs_kind(test4_obs_def)
   if (qc_index < 0) then
      test4_qc = 0
   else
      call get_qc(obs_in4, temp, qc_index)
      test4_qc = nint(temp(1))
   endif
else
   test4_qc       = test1_qc
   test4_obs_type = test1_obs_type
   test4_obs_time = test1_obs_time
   test4_obs_loc  = test1_obs_loc
endif

! first compare the integer values; it's quicker than location.
! in this case leave the QC until last so we can get statistics on why
! an obs is rejected.  if a lot of obs fail the other tests, then the input
! files aren't as identical as they're supposed to be.

if ((test1_obs_type /= test2_obs_type) .or. &
    (test1_obs_type /= test3_obs_type) .or. &
    (test1_obs_type /= test4_obs_type) .or. &
    (test2_obs_type /= test3_obs_type) .or. &
    (test2_obs_type /= test4_obs_type) .or. &
    (test3_obs_type /= test4_obs_type)) then
   num_rejected_other = num_rejected_other + 1
   return
endif

if ((test1_obs_time /= test2_obs_time) .or. &
    (test1_obs_time /= test3_obs_time) .or. &
    (test1_obs_time /= test4_obs_time) .or. &
    (test2_obs_time /= test3_obs_time) .or. &
    (test2_obs_time /= test4_obs_time) .or. &
    (test3_obs_time /= test4_obs_time)) then
   num_rejected_other = num_rejected_other + 1
   return
endif

if ((test1_obs_loc /= test2_obs_loc) .or. &
    (test1_obs_loc /= test3_obs_loc) .or. &
    (test1_obs_loc /= test4_obs_loc) .or. &
    (test2_obs_loc /= test3_obs_loc) .or. &
    (test2_obs_loc /= test4_obs_loc) .or. &
    (test3_obs_loc /= test4_obs_loc)) then
   num_rejected_other = num_rejected_other + 1
   return
endif

if ((test1_qc /= test2_qc) .or. &
    (test1_qc /= test3_qc) .or. &
    (test1_qc /= test4_qc) .or. &
    (test2_qc /= test3_qc) .or. &
    (test2_qc /= test4_qc) .or. &
    (test3_qc /= test4_qc)) then
   num_rejected_diffqc = num_rejected_diffqc + 1
   return
endif

! this assumes we have already tested for both qcs being equal, so it can
! just test one and know whether they are both over the threshold or not.
if (test1_qc > qc_threshold) then
   num_rejected_badqc = num_rejected_badqc + 1
   return
endif

! all match - good return.
all_good = .true.

end function all_good

!---------------------------------------------------------------------
! currently unused:

function horiz_location_equal(loc1,loc2)

! function to compare only the lat & lon and ignore the vert location.

type(location_type), intent(in) :: loc1, loc2
logical                         :: horiz_location_equal

real(r8) :: l1(3), l2(3)

l1 = get_location(loc1)
l2 = get_location(loc2)

horiz_location_equal = .false.

if ( abs(l1(1)  - l2(1) ) > epsilon(l1(1) ) ) return
if ( abs(l1(2)  - l2(2) ) > epsilon(l1(2) ) ) return

horiz_location_equal = .true.

end function horiz_location_equal

!---------------------------------------------------------------------
end program obs_common_subset

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

