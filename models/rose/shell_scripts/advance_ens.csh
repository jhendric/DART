#!/bin/tcsh

#   Multi-processor jobs must be submitted as batch, under PBS
### Job name
#PBS -N dart_rose
### Declare job non-rerunable
#PBS -r n
### Output files
#PBS -e dart_rose.err
#PBS -o dart_rose.log
### Queue name (small, medium, long, verylong)
#PBS -q small
#PBS -l nodes=1

### This job's working directory; must cd to it, or it will run in /home...
if ($?PBS_O_WORKDIR) then
   cd $PBS_O_WORKDIR
else
   setenv PBS_O_WORKDIR `pwd`
endif

### Output to confirm job characteristics
if ($?PBS_JOBNAME) then
   echo Running $PBS_JOBNAME on host `hostname`
else
   echo "Running on host "`hostname`
endif
echo Time is `date`
echo Directory is `pwd`
echo This job runs on the following processors:

### Define number of processors; # of lines in PBS_NODEFILE
if ($?PBS_NODEFILE) then
   setenv NPROCS `wc -l < $PBS_NODEFILE`
else
   setenv PBS_NODEFILE nodefile
   rm -f $PBS_NODEFILE
   set NPROCS = 10
   set inode = 1
   while($inode <= $NPROCS)
      echo node$inode >> $PBS_NODEFILE
      @ inode ++
   end
endif
cat "$PBS_NODEFILE"
echo This job has allocated $NPROCS nodes

# First line of filter_control should have number of model states to be integrated
set nensmbl = `head -1 filter_control`

# figure # batches of CAM runs to do, from # ensemble members and # processors
@ nbatch = $nensmbl / $NPROCS
if ($nensmbl % $NPROCS != 0 ) @ nbatch++
echo $nbatch batches will be executed

# Create a directory for each member to run in for namelists
set element = 0
set batch = 1
while($batch <= $nbatch)
   foreach node ( `cat $PBS_NODEFILE` )
      @ element++
      if ($element > $nensmbl) goto all_elements_done

      rsh $node "csh $PBS_O_WORKDIR/advance_model.csh $PBS_O_WORKDIR $element /scratch/local/tmp$user$$$element " &

      echo rsh $node \
         "csh $PBS_O_WORKDIR/advance_model.csh $PBS_O_WORKDIR $element /scratch/local/tmp$user$$$element " &

   end
# Another way to monitor progress.  batchflag has other info to start,
# so this echo can be removed and scripts will still work.
   echo waiting to finish batch $batch  >> $PBS_O_WORKDIR/batchflag
   wait
   @ batch++
end
all_elements_done:

# Wait for all *background* processes to finish up
wait

# signal to sync_submit.csh to continue
rm -f $PBS_O_WORKDIR/batchflag
