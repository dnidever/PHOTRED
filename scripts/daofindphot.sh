#!/bin/sh
#
export image=${1}
#
# Was the aperture file input?
if [ $# -lt 2 ]; then
  export apertures="photo.opt"
else
  export apertures=${2}  
fi
#
#  If required files do not exist, don't run this script.
#
if [ ! -s ${image}.fits ]; then
   echo "ERROR: input image [ ${image}.fits ] not found."
   exit 1
fi
if [ ! -s ${image}.opt ]; then
   echo "ERROR: ${image}.opt required to run DAOPHOT."
   exit 1
fi
if [ ! -s ${apertures} ]; then
   echo "ERROR: apertures file ${apertures} not found."
   exit 1
fi
#  Delete files to avoid errors in DAOPHOT.
#
rm ${image}.log      >& /dev/null
rm ${image}.coo      >& /dev/null
rm ${image}.ap       >& /dev/null
#
#  If daophot.opt does NOT exist, create it
#
if [ ! -s daophot.opt ]; then
  cp ${image}.opt daophot.opt
fi
#
daophot << END_DAOPHOT >> ${image}.log
OPTIONS
${image}.opt

ATTACH ${image}.fits
FIND
1,1
${image}.coo
y
PHOTOMETRY
${apertures}

${image}.coo
${image}.ap
EXIT
END_DAOPHOT
