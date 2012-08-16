#!/bin/csh
#
# DART software - Copyright 2004 - 2011 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# $Id$
#
# This is an example script for how to stage the files in CENTRALDIR
# in preparation for a 'perfect model' or OSSE.
#
#==============================================================================
# Set the commands so we can avoid problems with aliases, etc.
#==============================================================================

set   MOVE = '/usr/local/bin/mv -fv'
set   COPY = '/usr/local/bin/cp -fv --preserve=timestamps'
set   LINK = '/usr/local/bin/ln -fvs'
set REMOVE = '/usr/local/bin/rm -fr'

set   MOVE = '/bin/mv -fv'
set   COPY = '/bin/cp -fvp'
set   LINK = '/bin/ln -fvs'
set REMOVE = '/bin/rm -fr'

#==============================================================================
# Stage all the required files in CENTRALDIR
#==============================================================================

set CENTRALDIR = `pwd`

set NOAHDIR = /Users/thoar/svn/DART/devel/models/noah/src/hrldas-v3.3
set DARTDIR = /Users/thoar/svn/DART/devel/models/noah

${COPY} ${NOAHDIR}/Run/wrfinput.template    wrfinput  || exit 1
${COPY} ${NOAHDIR}/Run/namelist.hrldas             .  || exit 1
${COPY} ${NOAHDIR}/Run/Noah_hrldas_beta            .  || exit 1
${COPY} ${NOAHDIR}/Run/SOILPARM.TBL                .  || exit 1
${COPY} ${NOAHDIR}/Run/VEGPARM.TBL                 .  || exit 1
${COPY} ${NOAHDIR}/Run/GENPARM.TBL                 .  || exit 1
${COPY} ${NOAHDIR}/Run/URBPARM.TBL                 .  || exit 1

${COPY} ${DARTDIR}/work/obs_seq.in                 .  || exit 2
${COPY} ${DARTDIR}/work/input.nml                  .  || exit 2
${COPY} ${DARTDIR}/work/perfect_model_obs          .  || exit 2
${COPY} ${DARTDIR}/work/dart_to_noah               .  || exit 2
${COPY} ${DARTDIR}/work/noah_to_dart               .  || exit 2
${COPY} ${DARTDIR}/shell_scripts/run_pmo.csh       .  || exit 2
${COPY} ${DARTDIR}/shell_scripts/advance_model.csh .  || exit 2

# need a single noah restart file to be used as THE TRUTH.
# the input.nml:model_nml noah_netcdf_filename = 'restart.nc'
# the assimilate.csh scripts wants an ensemble member node

ln -sv ${NOAHDIR}/Run/hourly_output/RESTART.2004010107_DOMAIN1 restart.nc
ln -sv restart.nc restart.0001.nc

./noah_to_dart                || exit 3
${MOVE} dart_ics perfect_ics  || exit 4

echo
echo "CENTRALDIR is ${CENTRALDIR}"
echo "Configure     ${CENTRALDIR}/input.nml"
echo "Configure     ${CENTRALDIR}/namelist.hrldas"
echo "Configure     ${CENTRALDIR}/wrfinput"
echo "execute       ${CENTRALDIR}/run_pmo.csh"
echo

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

