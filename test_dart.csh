#!/bin/csh
#
# Data Assimilation Research Testbed -- DART
# Copyright 2004, 2005, Data Assimilation Initiative, University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
#
# <next three lines automatically updated by CVS, do not edit>
# $Id$
# $Source$
# $Name$

set SNAME = $0
set clobber

switch ( $#argv )
   case 0:
      # supplying no arguments -- echo usage not
      breaksw
   default:
      echo " "
      echo "usage: $SNAME:t"
      echo " "
      echo "This script compiles 'filter' for a wide range of models and then does"
      echo "relatively extensive tests of the L96 programs with a variety of options."
      echo " "
      echo "This must be run from the top-level 'DART' directory."
      echo " "
      echo "This is a pretty verbose process, so if you are logging the output,"
      echo "make sure you have plenty of space:"
      echo " "
      echo "./$SNAME:t |& tee DART_test.log"
      echo " "
      echo "can easily result in a 750 Kb log file"
      exit 1
      breaksw
endsw

if ( ! -d models/lorenz_96 ) then
   echo "models/lorenz_96 does not exist. $SNAME:t must be run from the top-level"
   echo "DART directory -- please try again."
   exit 2
else
   set DARTHOME = `pwd`
endif

echo "The top-level DART directory (DARTHOME) is $DARTHOME"

#----------------------------------------------------------------------
# Compile 'filter' for a wide range of models.
#----------------------------------------------------------------------

@ makenum  = 1
@ modelnum = 101
foreach MODEL ( 9var MITgcm_annulus bgrid_solo cam wrf pe2lyr \
            lorenz_04 lorenz_63 lorenz_84 lorenz_96 lorenz_96_2scale )

    echo "-------------------------------------------"
    echo "Compiling $MODEL at "`date`
    echo ""

    cd ${DARTHOME}/models/${MODEL}/work
    rm -f *.o *.mod filter Makefile input.nml.filter_default 
    csh mkmf_filter   || exit $makenum
    make              || exit $modelnum
    rm -f *.o *.mod filter Makefile input.nml.filter_default

   @ makenum  = $makenum  + 1
   @ modelnum = $modelnum + 1
end

#----------------------------------------------------------------------
# Lots of tests for L96
#----------------------------------------------------------------------

echo "-------------------------------------------"
echo "Testing lorenz_96 (L96) at "`date`
echo ""

cd ${DARTHOME}/models/lorenz_96/work

# Make sure that all .o, .mod and executables are gone
rm -rf *.o *.mod assim_region create_fixed_network_seq create_obs_seq filter
rm -rf integrate_model perfect_model_obs

# Begin by compiling all programs; need to stop if an error is detected
csh mkmf_assim_region
make                              || exit
csh mkmf_create_fixed_network_seq || exit
make                              || exit
csh mkmf_create_obs_sequence      || exit
make                              || exit
csh mkmf_filter                   || exit
make                              || exit
csh mkmf_integrate_model          || exit
make                              || exit
csh mkmf_perfect_model_obs        || exit
make                              || exit

# Setup appropriate namelists for a 1000-step test run
# Begin by modifying perfect model namelist and getting rid of all other stuff
echo ':0'                              > vi_script
echo '/start_from_restart'            >> vi_script
echo ':s/false/true/'                 >> vi_script
echo '/output_restart'                >> vi_script
echo ':s/false/true/'                 >> vi_script
echo '/ensemble_manager_nml'          >> vi_script
echo ':.,$ delete'                    >> vi_script
echo ':wq'                            >> vi_script
vi -s vi_script input.nml.perfect_model_obs_default

# Prepend this to the filter namelist and create input.nml
cat input.nml.perfect_model_obs_default input.nml.filter_default > input.nml

# Need to modify rest of input.nml for test run
echo ':0'                              > vi_script
echo '/ens_size'                      >> vi_script
echo ':s/20/80/'                      >> vi_script
echo '/start_from_restart'            >> vi_script
echo ':s/false/true/'                 >> vi_script
echo '/output_restart'                >> vi_script
echo ':s/false/true/'                 >> vi_script
echo '/num_output_state_members'      >> vi_script
echo ':s/0/20/'                       >> vi_script
echo '/num_groups'                    >> vi_script
echo ':s/1/4/'                        >> vi_script
echo '/cutoff'                        >> vi_script
echo ':s/0.2/1000000.0/'              >> vi_script
echo '/cov_inflate'                   >> vi_script
echo ':s/-1.0/1.05/'                  >> vi_script
echo '/read_binary_restart_files'     >> vi_script
echo ':s/true/false/'                 >> vi_script
echo ':wq'                            >> vi_script
vi -s vi_script input.nml

# Create an obs_sequence file
rm -rf obs_seq.in obs_seq.out obs_seq.final
./create_obs_sequence < ../random_obs.input || exit
echo 'set_def.out'		> temp_input
echo '1'                       >> temp_input
echo '1000'                    >> temp_input
echo '0 0'                     >> temp_input
echo '0 3600'                  >> temp_input
echo 'obs_seq.in'              >> temp_input

./create_fixed_network_seq < temp_input      || exit

# Run the perfect model and the filter
./perfect_model_obs  || exit
./filter             || exit

# Need to do visual matlab inspection of this output for now
# plot_total_err
matlab -nojvm

#-----------------------------------------------------------------------
echo "-----------------------------------------------------------------"
echo "Prepare to set up a sequence of 10 hour runs for testing"
echo "In first case, just do 10 days and output filter restarts in both"
echo "the single file and multiple file format for later testing"
echo "-----------------------------------------------------------------"
#-----------------------------------------------------------------------

# Create an obs_sequence file for 10 hour tests
rm -rf obs_seq.in obs_seq.out obs_seq.final
./create_obs_sequence < ../random_obs.input  || exit
echo 'set_def.out'              > temp_input
echo '1'                       >> temp_input
echo '10'                      >> temp_input
echo '0 0'                     >> temp_input
echo '0 3600'                  >> temp_input
echo 'obs_seq.in'              >> temp_input
echo ':wq'                     >> temp_input
./create_fixed_network_seq < temp_input      || exit

# Need to start from binary at end of previous long run
mv perfect_restart  perfect_ics.spun_up
mv filter_restart   filter_ics.spun_up

# To get reproducing across a variety of options, need num_domains 2 or greater
# And binary restart files
# Start from previous restart and create new restart
echo ':0'                                  > vi_script
echo '/num_domains'                       >> vi_script
echo ':s/1/3/'                            >> vi_script
echo ':0'                                 >> vi_script
echo '/read_binary_restart_files'         >> vi_script
echo ':s/false/true/'                     >> vi_script
echo ':0'                                 >> vi_script
echo '/restart_in_file_name'              >> vi_script
echo ':s/perfect_ics/perfect_ics.spun_up' >> vi_script
echo '/filter_nml'                        >> vi_script
echo '/restart_in_file_name'              >> vi_script
echo ':s/filter_ics/filter_ics.spun_up'   >> vi_script
echo ':wq'                                >> vi_script
vi -s vi_script input.nml

# Run the perfect model and the filter to produce single restart file
./perfect_model_obs   || exit
./filter              || exit
mv perfect_restart  perfect_ics.10hour
mv filter_restart   filter_ics.10hour

# Now run the filter again to produce multiple restart files
echo ':0'                                 > vi_script
echo '/single_restart_file_out'          >> vi_script
echo ':s/true/false/'                    >> vi_script
echo ':wq'                               >> vi_script
vi -s vi_script input.nml
# Run the filter, the files will end up in filter_restart.01
./filter  || exit


# Now do a second 10 hour run from the end of the first 10 hour
# Also change the filter back to only produce single restart file
echo ':0'                                          > vi_script
echo '/restart_in_file_name'                      >> vi_script
echo ':s/perfect_ics.spun_up/perfect_ics.10hour'  >> vi_script
echo '/filter_nml'                                >> vi_script
echo '/restart_in_file_name'                      >> vi_script
echo ':s/filter_ics.spun_up/filter_ics.10hour'    >> vi_script
echo '/single_restart_file_out'                   >> vi_script
echo ':s/false/true/'                             >> vi_script
echo ':wq'                                        >> vi_script
vi -s vi_script input.nml
./perfect_model_obs  || exit
./filter             || exit


# Do some tests on perfect_model_obs options first
# Set up baseline output file
mv perfect_restart  perfect_restart.baseline
mv obs_seq.out      obs_seq.out.baseline
mv True_State.nc    True_State.nc.baseline

# Test storing the ensemble on disk instead of in core
echo ':0'                              > vi_script
echo '/ensemble_manager'              >> vi_script
echo '/in_core'                       >> vi_script
echo ':s/true/false/'                 >> vi_script
echo ':wq'                            >> vi_script
vi -s vi_script input.nml

./perfect_model_obs  || exit
mv perfect_restart  perfect_restart.out_of_core
mv obs_seq.out      obs_seq.out.out_of_core
mv True_State.nc    True_State.nc.out_of_core
diff perfect_restart.out_of_core  perfect_restart.baseline || exit
diff obs_seq.out.out_of_core      obs_seq.out.baseline     || exit
#diff True_State.nc.out_of_core    True_State.nc.baseline   ||exit

# Test the two async options
# Need to get the scripts (problem here because script names are machine dependent)
cp -p ../shell_scripts/*.csh .
cp -p ${DARTHOME}/shell_scripts/advance_ens_fisher.csh     advance_ens.csh
cp -p ${DARTHOME}/shell_scripts/assim_filter_fisher.csh   assim_filter.csh
cp -p ${DARTHOME}/shell_scripts/filter_server_fisher.csh filter_server.csh

# Change to async 2 and go back to in_core for ensemble
echo ':0'                             > vi_script
echo '/perfect_model_obs_nml'        >> vi_script
echo '/async'                        >> vi_script
echo ':s/0/2/'                       >> vi_script
echo '/in_core'                      >> vi_script
echo ':s/false/true/'                >> vi_script
echo ':wq'                           >> vi_script
vi -s vi_script input.nml

./perfect_model_obs  || exit
mv perfect_restart  perfect_restart.2
mv obs_seq.out      obs_seq.out.2
mv True_State.nc    True_State.nc.2
diff perfect_restart.2  perfect_restart.baseline  || exit
diff obs_seq.out.2      obs_seq.out.baseline      || exit
#diff True_State.nc.2    True_State.nc.baseline    || exit

# Now try option 3 with a filter_server
echo ':0'                             > vi_script
echo '/perfect_model_obs_nml'        >> vi_script
echo '/async'                        >> vi_script
echo ':s/2/3/'                       >> vi_script
echo ':wq'                           >> vi_script
vi -s vi_script input.nml

./filter_server.csh &
./perfect_model_obs  || exit
mv perfect_restart  perfect_restart.3
mv obs_seq.out      obs_seq.out.3
cp True_State.nc    True_State.nc.3
diff perfect_restart.3  perfect_restart.baseline  || exit
diff obs_seq.out.3      obs_seq.out.baseline      || exit
#diff True_State.nc.3    True_State.nc.baseline    || exit

# Next, start checking filter options for this case

# Begin by checking single versus multiple restarts
cp obs_seq.out.baseline obs_seq.out
./filter  || exit
mv obs_seq.final        obs_seq.final.baseline
mv filter_restart       filter_restart.baseline
mv assim_tools_restart  assim_tools_restart.baseline
cp Prior_Diag.nc        Prior_Diag.nc.baseline
cp Posterior_Diag.nc    Posterior_Diag.nc.baseline

matlab -nojvm

# Now do a run with multiple input files, NOTE: they are filter_restart.00??
echo ':0'                                     > vi_script
echo '/filter_nml'                           >> vi_script
echo '/restart_in_file_name'                 >> vi_script
echo ':s/filter_ics.10hour/filter_restart/'  >> vi_script
echo '/single_restart_file_in'               >> vi_script
echo ':s/true/false/'                        >> vi_script
echo ':wq'                                   >> vi_script
vi -s vi_script input.nml

./filter  || exit
mv obs_seq.final        obs_seq.final.in_files
mv filter_restart       filter_restart.in_files
mv assim_tools_restart  assim_tools_restart.in_files
cp Prior_Diag.nc        Prior_Diag.nc.in_files
cp Posterior_Diag.nc    Posterior_Diag.nc.in_files

diff obs_seq.final.in_files        obs_seq.final.baseline       || exit
diff filter_restart.in_files       filter_restart.baseline      || exit
diff assim_tools_restart.in_files  assim_tools_restart.baseline || exit
#diff Prior_Diag.nc.in_files        Prior_Diag.nc.baseline
#diff Posterior_Diag.nc.in_files    Posterior_Diag.nc.baseline


# Next switch the number of domains from 3 to 5
echo ':0'                                      > vi_script
echo '/num_domains'                           >> vi_script
echo ':s/3/5/'                                >> vi_script
echo ':wq'                                    >> vi_script
vi -s vi_script input.nml

./filter  || exit

diff obs_seq.final        obs_seq.final.baseline        || exit
diff filter_restart       filter_restart.baseline       || exit
#diff assim_tools_restart  assim_tools_restart.baseline  || exit
#diff Prior_Diag.nc        Prior_Diag.nc.baseline
#diff Posterior_Diag.nc    Posterior_Diag.nc.baseline

# Next switch the number of domains back to 3; try parallel option 2
echo ':0'                                      > vi_script
echo '/do_parallel'                           >> vi_script
echo ':s/0/2/'                                >> vi_script
echo '/num_domains'                           >> vi_script
echo ':s/5/3/'                                >> vi_script
echo ':wq'                                    >> vi_script
vi -s vi_script input.nml

./filter || exit

diff obs_seq.final        obs_seq.final.baseline        || exit
diff filter_restart       filter_restart.baseline       || exit
diff assim_tools_restart  assim_tools_restart.baseline  || exit
#diff Prior_Diag.nc        Prior_Diag.nc.baseline        || exit
#diff Posterior_Diag.nc    Posterior_Diag.nc.baseline    || exit


# Try parallel option 3
echo ':0'                                      > vi_script
echo '/do_parallel'                           >> vi_script
echo ':s/2/3/'                                >> vi_script
echo ':wq'                                    >> vi_script
vi -s vi_script input.nml

csh ./filter_server.csh &
./filter  || exit

diff obs_seq.final        obs_seq.final.baseline        || exit
diff filter_restart       filter_restart.baseline       || exit
diff assim_tools_restart  assim_tools_restart.baseline  || exit
#diff Prior_Diag.nc        Prior_Diag.nc.baseline
#diff Posterior_Diag.nc    Posterior_Diag.nc.baseline



# Go back to parallel option 0 and proceed
echo ':0'                                      > vi_script
echo '/do_parallel'                           >> vi_script
echo ':s/3/0/'                                >> vi_script
echo ':wq'                                    >> vi_script
vi -s vi_script input.nml

./filter  || exit

diff obs_seq.final        obs_seq.final.baseline        || exit
diff filter_restart       filter_restart.baseline       || exit
diff assim_tools_restart  assim_tools_restart.baseline  || exit
#diff Prior_Diag.nc        Prior_Diag.nc.baseline
#diff Posterior_Diag.nc    Posterior_Diag.nc.baseline

echo ""
echo "Testing complete  at "`date`
echo "-------------------------------------------"

