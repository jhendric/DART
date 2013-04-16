#!/bin/csh

# scaling of transpose on Yellowstone
 
#BSUB -J timing 	
#BSUB -o timing.%J.log
#BSUB -e timing.%J.err
#BSUB -q regular 
#BSUB -n 2 
#BSUB -R "span[ptile=16]"
#BSUB -P P86850054 
#BSUB -W 0:40

# glade is current running dir.

# These options are set in Yellowstone by default
# enable 64k pagesizes - this needs doing only once
#ldedit -bdatapsize=64K -bstackpsize=64K -btextpsize=64K ./filter

# pin tasks to processors
setenv TARGET_CPU_LIST "-1"

mpirun.lsf ./timing

exit 0

