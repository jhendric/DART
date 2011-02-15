#!/bin/tcsh 
#
# SCRIPT:	job_setup.csh
# AUTHOR:	T. R. Whitcomb
#           Naval Research Laboratory
#
# Sets up the configuration (e.g. Linux modules) for a script run
# in a resource manager - this file is not executed, but is sourced
# by various run scripts.
######

# These are based on the ACESGrid setup at MIT
if ( -f /etc/profile.d/modules.csh ) then
    source /etc/profile.d/modules.csh
endif
module load mpich/pgi
module load mpiexec