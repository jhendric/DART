! Data Assimilation Research Testbed -- DART
! Copyright 2004-2007, Data Assimilation Research Section
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

program wakeup_filter

! spread out on the processors in the same order as the executable
! model and filter programs, and echo into the fifo (named pipe)
! a message to wake up the sleeping filter program.
!

! <next few lines automatically updated by version control software, do not edit>
! $Revision: 1.1 $
! $Date: 2007/03/06 22:15:47 $
! $Id: wakeup_filter.f90,v 1.1 2007/03/06 22:15:47 thoar Exp $


use mpi_utilities_mod, only : initialize_mpi_utilities, finalize_mpi_utilities, &
                              restart_task




call initialize_mpi_utilities("Wakeup_Filter")

call restart_task()

call finalize_mpi_utilities()


end program wakeup_filter