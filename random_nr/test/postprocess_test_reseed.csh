#!/bin/csh
#
# DART software - Copyright 2004 - 2011 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# $Id$

foreach SEEDSEP ( 2 60 3600 21600 43200 86400 )

   cd `printf seeds_%05d $SEEDSEP`
   set fname = `printf ran_test_r8_%05d_full.out $SEEDSEP`

   # focus on just the following trials 

   grep 'seed'             $fname >! seeds
   grep 'uniform     1000' $fname >! u1000
   grep 'gaussian    1000' $fname >! g1000
   grep 'gaussian 1000000' $fname >! g1000000

   # strip out the string and leave the numbers to ingest into matlab

   sed s/"trial "//           seeds    >! bob
   sed s/"new seed: "//       bob      >! seeds
   sed s/"uniform     1000"// u1000    >! unif_1000
   sed s/"gaussian    1000"// g1000    >! gauss_1000
   sed s/"gaussian 1000000"// g1000000 >! gauss_1000000

   rm u1000 g1000 g1000000 bob

   cd ..

end

exit $status

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

