#!/bin/csh

# timing transpose

#BSUB -J timing
#BSUB -o timing.%J.log
#BSUB -q debug 
#BSUB -n 2 
#BSUB -P 86850054 
#BSUB -W 0:03

# enable 64k pagesizes - this needs doing only once
ldedit -bdatapsize=64K -bstackpsize=64K -btextpsize=64K ./timing

# pin tasks to processors
setenv TARGET_CPU_LIST "-1"

mpirun.lsf /usr/local/bin/launch ./timing

exit 0

