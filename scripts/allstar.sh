#!/bin/sh
#
#############################################################################
#
# ALLSTAR.SH - Script to run ALLSTAR after a PSF already exists
#
# USAGE - ./getpsf.sh <image name w/o extension>
#
# HISTORY - 01/07/2005 : Created by altering daostep1.sh daostep2.sh.
#
#           02/03/2005 : Modified - run another pass of goodpsf after
#                        subtracting neighbor stars.
#
#           02/03/2005 : Added the goodpsf.f F77 program to the end of 
#                        this script file and commented out.
#
#           02/07/2005 : Added the preapcor.sh script at the end
#
#           03/31/2005 : Modified to use LSTFILTER and added lstfilter.f
#                        F77 program to the end of this script.
#
#           04/19/2005 : Added WHILE loop in 2nd step for initial filtering.
#
#
# Script made by Sangmo Tony Sohn
#
#############################################################################
#
daophot="/net/astro/bin/daophot"
allstar="/net/astro/bin/allstar"
#goodpsf="/net/halo/bin/goodpsf"
#lstfilter="/net/halo/bin/lstfilter"
export image=${1}
#
#  If required files do not exist, don't run this script.
#
if [ ! -s ${image}.fits ]; then
   echo "ERROR: input image [ ${image}.fits ] not found."
   exit 1
fi
if [ ! -s ${image}.psf ]; then
   echo "ERROR: input image [ ${image}.psf ] not found."
   exit 1
fi
if [ ! -s ${image}.opt ]; then
   echo "ERROR: ${image}.opt required to run DAOPHOT."
   exit 1
fi
if [ ! -s apcor.opt ]; then
   echo "ERROR: apcor.opt required to run DAOPHOT."
   exit 1
fi
if [ ! -s ${image}.als.opt ]; then
   echo "ERROR: ${image}.als.opt required to run ALLSTAR."
   exit 1
fi
#
echo "Running ALLSTAR on ${image}.fits :"
#
###########################################
###                                     ###
### Perform PSF photometry with ALLSTAR ###
###                                     ###
###########################################
#
# If allstar.opt does NOT exist, create it
#
if [ ! -s allstar.opt ]; then
  cp ${image}.als.opt allstar.opt
fi
#
#rm allstar.inp    >& /dev/null
rm ${image}.als.inp >& /dev/null
rm ${image}.als   >& /dev/null
rm ${image}s.fits >& /dev/null
#
echo "" >> ${image}.log
echo "===================" >> ${image}.log
echo "==  Run ALLSTAR. ==" >> ${image}.log
echo "===================" >> ${image}.log
echo "" >> ${image}.log
#
cat ${image}.als.opt > ${image}.als.inp
echo $image'.fits'  >> ${image}.als.inp
echo $image'.psf'   >> ${image}.als.inp
echo $image'.ap'    >> ${image}.als.inp
echo $image'.als'   >> ${image}.als.inp
echo $image's.fits' >> ${image}.als.inp
allstar < ${image}.als.inp >> ${image}.log
#rm allstar.inp >& /dev/null
#
