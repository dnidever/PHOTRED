#!/bin/sh
#
#############################################################################
#
# DAOPHOT.SH - Script to perform automated PSF photometry using DAOPHOT II
#
# USAGE - ./daophot.sh <image name w/o extension>
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
if [ ! -s ${image}.opt ]; then
   echo "ERROR: ${image}.opt required to run DAOPHOT."
   exit 1
fi
if [ ! -s photo.opt ]; then
   echo "ERROR: photo.opt required to run DAOPHOT."
   exit 1
fi
if [ ! -s ${image}.als.opt ]; then
   echo "ERROR: ${image}.als.opt required to run ALLSTAR."
   exit 1
fi
if [ ! -s ${image}.psf ]; then
   echo "ERROR: ${image}.psf required to run ALLSTAR."
   exit 1
fi
#
echo "Starting photometry on ${image}.fits :"
#
###############################################################
###                                                         ###
### 1st step : select PSF candidates and obtain initial PSF ###
###                                                         ###
###############################################################
#
echo "   Running 1st step ...."
#
#  Delete files to avoid errors in DAOPHOT.
#
rm ${image}.log      >& /dev/null
rm ${image}.coo      >& /dev/null
rm ${image}.ap      >& /dev/null
rm ${image}.psf.tmp      >& /dev/null
# Temporarily move the psf file out of the way
#   otherwise we can't do aperture photometry the normal way
mv ${image}.psf ${image}.psf.tmp     >& /dev/null
#
echo "===================================" >> ${image}.log
echo "== Step 1 : Detection            ==" >> ${image}.log
echo "===================================" >> ${image}.log
echo "" >> ${image}.log
#
#  Copy the {image}.opt file to daophot.opt for safety.
#
#rm daophot.opt >& /dev/null
#cp ${image}.opt daophot.opt
#
#  If daophot.opt does NOT exist, create it
#
if [ ! -s daophot.opt ]; then
  cp ${image}.opt daophot.opt
fi
#
#  Run DAOPHOT to perform FIND, PHOTOMETRY, and PICKPSF on frame.
#  You may change the number of stars in PICKPSF.
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


${image}.coo
${image}.ap
EXIT
END_DAOPHOT


######################################################
###                                                ###
### 2nd step : perform PSF photometry with ALLSTAR ###
###                                                ###
######################################################
#
echo "   Running 2nd step ...."
# Move the psf file back
mv ${image}.psf.tmp ${image}.psf     >& /dev/null
#
# If allstar.opt does NOT exist, create it
#
if [ ! -s allstar.opt ]; then
  cp ${image}.als.opt allstar.opt
fi
#
#rm allstar.inp    >& /dev/null
rm ${image}.als.inp    >& /dev/null
rm ${image}.als   >& /dev/null
rm ${image}s.fits >& /dev/null
#
echo "" >> ${image}.log
echo "===========================" >> ${image}.log
echo "== Step 2 : Run ALLSTAR. ==" >> ${image}.log
echo "===========================" >> ${image}.log
echo "" >> ${image}.log
#
cat ${image}.als.opt > ${image}.als.inp
echo $image'.fits'  >> ${image}.als.inp
echo $image'.psf'   >> ${image}.als.inp
echo $image'.ap'    >> ${image}.als.inp
echo $image'.als'   >> ${image}.als.inp
echo $image's.fits' >> ${image}.als.inp
echo exit >> ${image}.als.inp
allstar < ${image}.als.inp >> ${image}.log
#rm allstar.inp >& /dev/null

echo ""
