#!/bin/sh
export image=${1}
rm ${image}.log      >& /dev/null
rm ${image}.coo      >& /dev/null
rm ${image}.ap       >& /dev/null
rm ${image}a.coo     >& /dev/null
rm ${image}a.ap      >& /dev/null
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
photo.opt

${image}.coo
${image}.ap
EXIT
END_DAOPHOT

# Get the magnitude limit
maglim=`grep "Estimated magnitude limit" ${image}.log | awk '{print $6-1.0}'`

# Pick PSF stars
daophot << END_DAOPHOT >> ${image}.log
OPTIONS
${image}.opt

ATTACH ${image}.fits
PICKPSF
${image}.ap
100,${maglim}
${image}a.coo
PHOTOMETRY
photo.opt

${image}a.coo
${image}a.ap
EXIT
END_DAOPHOT
