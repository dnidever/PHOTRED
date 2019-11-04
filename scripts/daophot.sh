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
export baseworkdir=${2}
#
#  If required files do not exist, don't run this script.
#
if [ ! -s ${image}.fits ] && [ ! -s ${image}.fits.fz ]; then
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
if [ ! -s apcor.opt ]; then
   echo "ERROR: apcor.opt required to run DAOPHOT."
   exit 1
fi
if [ ! -s ${image}.als.opt ]; then
   echo "ERROR: ${image}.als.opt required to run ALLSTAR."
   exit 1
fi
if [ ! -s goodpsf.pro ]; then
   echo "ERROR: GOODPSF program required to filter bad PSF stars."
   exit 1
fi
if [ ! -s lstfilter.py ]; then
   echo "ERROR: LSTFILTER program required to filter bad PSF stars."
   exit 1
fi
###################################################
# Working in temporary directory, copy files there
if [ -n ${baseworkdir} ]; then
   export curdir=`pwd`
   # Create base working directory 
   if [ ! -s ${baseworkdir} ] && [ ! -d ${baseworkdir} ]; then
      mkdir ${baseworkdir}
   fi
   # Create temporary directory
   export workdir=`mktemp -d --tmpdir=${baseworkdir} dao.XXXXXX`
   echo "Working in temporary directory ${workdir}"
   # Copy the files that we need
   #  fits, photo.opt, apcor.opt, opt, als.opt, goodpsf.pro, lstfilter.py
   cp -f photo.opt apcor.opt goodpsf.pro lstfilter.py ${workdir} >& /dev/null
   cp -f ${image}.opt ${image}.als.opt ${workdir} >& /dev/null
   if [ -s ${image}.fits ]; then
     cp -f ${image}.fits ${workdir} >& /dev/null
   else
     cp -f ${image}.fits.fz ${workdir} >& /dev/null
   fi
   # Copy the cmn files
   if [ -s ${image}.cmn.lst ]; then
     cp -f ${image}.cmn.ap ${image}.cmn.coo ${image}.cmn.log ${image}.lst ${workdir} >& /dev/null
   fi
   # Go there
   cd ${workdir}
fi
#
# Using fpack compressed file
export fpackfile=0
export fullfile=${image}.fits
if [ ! -s ${image}.fits ] && [ -s ${image}.fits.fz ]; then
   echo "Temporarily uncompressing ${image}.fits.fz"
   funpack ${image}.fits.fz
   export fpackfile=1
   export fullfile=${image}.fits.fz
fi
#
echo "Starting photometry on ${fullfile} :"
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
rm ${image}.ap       >& /dev/null
rm ${image}.lst      >& /dev/null
rm ${image}.lst1     >& /dev/null
rm ${image}.psf      >& /dev/null
rm ${image}.nei      >& /dev/null
rm ${image}.psf.log  >& /dev/null
#
echo "===================================" >> ${image}.log
echo "== Step 1 : Construct PSF ver. 1 ==" >> ${image}.log
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

# Get the magnitude limit
maglim=`grep "Estimated magnitude limit" ${image}.log | awk '{print $6-1.0}'`
#maglim=`grep "Estimated magnitude limit" ${image}.log | awk '{print $6-3.0}'`
echo "Using magnitude limit = ${maglim}" >> ${image}.log

# Make psfini.coo file.  This is the list to start with for PSF stars.
cp ${image}.ap ${image}.psfini.ap >& /dev/null
# Remove "bad" sources using the CONFIRMED CELESTIAL SOURCES list if it exists
# The .cmn.lst file is created by PHOTRED_COMMONSOURCES.PRO
if [ -s ${image}.cmn.lst ]; then
  echo "srcfilter,'${image}.ap','${image}.cmn.lst','${image}.psfini.ap'" | idl >& /dev/null
fi

# Pick PSF stars
daophot << END_DAOPHOT >> ${image}.log
OPTIONS
${image}.opt

ATTACH ${image}.fits
PICKPSF
${image}.psfini.ap
100,${maglim}
${image}.lst
EXIT
END_DAOPHOT
#
#  NOT ENOUGH PSF stars, increase the magnitude limit
#
psfstars=`wc -l ${image}.lst | awk '{print $1-3}'`
echo "Number of PSF stars = ${psfstars}" >> ${image}.log

# Increase Magnitude limit
if [ ${psfstars} -lt 12 ]; then
maglim2=`echo $maglim | awk '{print $1+3.0}'`
echo "Not enough PSF stars. Increasing the magnitude limit to ${maglim2}" >> ${image}.log
rm ${image}.lst      >& /dev/null

# Pick PSF stars AGAIN
daophot << END_DAOPHOT >> ${image}.log
OPTIONS
${image}.opt

ATTACH ${image}.fits
PICKPSF
${image}.psfini.ap
100,${maglim2}
${image}.lst
EXIT
END_DAOPHOT

fi
# End Increase Magnitude Limit
#
#  Filter bad sources from ${image}.lst and copy it as ${image}.lst1 using LSTFILTER
#
echo " " >> ${image}.log 
echo "LSTFILTER ${image}.lst" >> ${image}.log
## Always run the LOCAL version of "lstfilter.py"
./lstfilter.py ${image}.lst ${image}.lst1  >> ${image}.log
#
#  Run DAOPHOT to construct a first version of PSF.
#
daophot << END_DAOPHOT >> ${image}.psf.log
OPTIONS
${image}.opt

ATTACH ${image}.fits
PSF
${image}.ap
${image}.lst1
${image}.psf

EXIT
END_DAOPHOT
# The last blank line is to overwrite the .nei file
#
# If the image names has extra dots "." in it (i.e. F1.obj1140_4.fits) then PSF writes the
# nei file to F1.nei.  Copy this file to the right filename.  If they are the same
# it will just give an error, but still be okay.
# The .nei file actually isn't used for anything, but I'm not sure if this is necessary.
# COMMENTING THIS OUT SINCE IT'S NOT NECESSARY.
#neifile=`grep "File with PSF stars and neighbors" ${image}.psf.log | awk '{print $8}'`
#cp -f ${neifile} ${image}.nei >& /dev/null
#
#
#################################################################################
###                                                                           ###
### 2nd step : filter out bad PSF stars and construct a second version of PSF ###
###                                                                           ###
#################################################################################
#
echo "   Running 2nd step ...."
#
rm ${image}.lst2     >& /dev/null
rm ${image}.lst1.chi ${image}.lst2.chi >& /dev/null
#
echo "" >> ${image}.log
echo "===================================" >> ${image}.log
echo "== Step 2 : Construct PSF ver. 2 ==" >> ${image}.log
echo "===================================" >> ${image}.log
echo "" >> ${image}.log
#
# This gets the lines from the log file that has the chi for each star
# and gives "saturated" and "defective" stars a chi=1.00
# and saves it to image.lst1.chi
ln1=`grep -n "Profile errors:" ${image}.psf.log | sed 's/:/ /g' | awk '(NR==1){print $1+2}'`
ln2=`grep -n "File with PSF stars and neighbors =" ${image}.psf.log | sed 's/:/ /g' | awk '(NR==1){print $1-4}'`
nl=$((ln2-ln1+1))
tail -n +${ln1} ${image}.psf.log | head -n ${nl} | sed -e 's/saturated/ 1.000   /g' -e 's/defective/ 1.000   /g' > ${image}.lst1.chi
#
cp ${image}.lst1 ${image}.lst2
cp ${image}.lst1.chi ${image}.lst2.chi
nbad=100
#
#######
#
#  The WHILE loop for rejecting stars with either '*' or '?' next to them.
#
#  This re-runs DAOPHOT PSF and keeps rejecting bad stars, until only good stars
#  are left.  
#
while [ $nbad -gt 0 ]; do 
echo "goodpsf,'${image}.lst2','${image}.lst2.chi','${image}.lstt'" | idl >& /dev/null
#goodpsf << END_GOODPSF >& /dev/null
#${image}.lst2
#${image}.lst2.chi
#${image}.lstt
#END_GOODPSF
#
rm ${image}.psf.log
rm ${image}.lst2 ${image}.lst2.chi
mv ${image}.lstt ${image}.lst2
#
daophot << END_DAOPHOT >> ${image}.psf.log
OPTIONS
${image}.opt

ATTACH ${image}.fits
PSF
${image}.ap
${image}.lst2
${image}.psf


EXIT
END_DAOPHOT
#
ln1=`grep -n "Profile errors:" ${image}.psf.log | sed 's/:/ /g' | awk '(NR==1){print $1+2}'`
ln2=`grep -n "File with PSF stars and neighbors =" ${image}.psf.log | sed 's/:/ /g' | awk '(NR==1){print $1-4}'`
nl=$((ln2-ln1+1))
tail -n +${ln1} ${image}.psf.log | head -n ${nl} | sed -e 's/saturated/ 1.000   /g' -e 's/defective/ 1.000   /g' > ${image}.lst2.chi
nbad=`grep -e '?' -e '*' ${image}.lst2.chi | wc -l`
done
#
# End of WHILE loop
#
#######
#
####################################################################
###                                                              ###
### 3rd step : subtract neighbor stars and construct a final PSF ###
###                                                              ###
####################################################################
#
echo "   Running 3rd step ...."
#
rm ${image}.grp      >& /dev/null
rm ${image}.nst      >& /dev/null
rm ${image}a.fits    >& /dev/null
rm ${image}.lst2.chi >& /dev/null
rm ${image}.lst3     >& /dev/null
rm ${image}.lst3.chi >& /dev/null
rm ${image}.plst     >& /dev/null
rm ${image}.psf.log  >& /dev/null
rm ${image}.plst.chi >& /dev/null
rm ${image}.psf.log  >& /dev/null
#
echo "" >> ${image}.log
echo "=============================================" >> ${image}.log
echo "== Step 3 : Construct final version of PSF ==" >> ${image}.log
echo "=============================================" >> ${image}.log
echo "" >> ${image}.log
#
#  Run DAOPHOT to subtract neighbor stars.
#
daophot << END_DAOPHOT >> ${image}.log
OPTIONS
${image}.opt

ATTACH ${image}.fits
GROUP
${image}.nei
${image}.psf
5.
${image}.grp
NSTAR
${image}.psf
${image}.grp
${image}.nst
SUBSTAR
${image}.psf
${image}.nst
y
${image}.lst2
${image}a.fits

EXIT
END_DAOPHOT
# Sometimes there's a segmentation fault on NSTAR
if [ ! -s ${image}a.fits ]; then
   echo "Possible segmentation fault in NSTAR.  Rerunning SUBSTAR to create a.fits"
daophot << END_DAOPHOT >> ${image}.log
OPTIONS
${image}.opt

ATTACH ${image}.fits
SUBSTAR
${image}.psf
${image}.nst
y
${image}.lst2
${image}a.fits

EXIT
END_DAOPHOT
fi
#
#  Run DAOPHOT to create ${image}.psf.log
#
daophot << END_DAOPHOT >> ${image}.psf.log
OPTIONS
${image}.opt

ATTACH ${image}a.fits
PSF
${image}.ap
${image}.lst2
${image}.psf


EXIT
END_DAOPHOT
#
ln1=`grep -n "Profile errors:" ${image}.psf.log | sed 's/:/ /g' | awk '(NR==1){print $1+2}'`
ln2=`grep -n "File with PSF stars and neighbors =" ${image}.psf.log | sed 's/:/ /g' | awk '(NR==1){print $1-4}'`
nl=$((ln2-ln1+1))
tail -n +${ln1} ${image}.psf.log | head -n ${nl} | sed -e 's/saturated/ 1.000   /g' -e 's/defective/ 1.000   /g' > ${image}.lst2.chi
#
cp ${image}.lst2 ${image}.lst3
cp ${image}.lst2.chi ${image}.lst3.chi
nbad=100
#
#######
#
#  The WHILE loop for rejecting stars with either '*' or '?' next to them.
#
while [ $nbad -gt 0 ]; do 
echo "goodpsf,'${image}.lst3','${image}.lst3.chi','${image}.plst'" | idl >& /dev/null
#goodpsf << END_GOODPSF >& /dev/null
#${image}.lst3
#${image}.lst3.chi
#${image}.plst
#END_GOODPSF
#
rm ${image}.psf.log
rm ${image}.lst3 ${image}.lst3.chi
mv ${image}.plst ${image}.lst3
#
daophot << END_DAOPHOT >> ${image}.psf.log
OPTIONS
${image}.opt

ATTACH ${image}a.fits
PSF
${image}.ap
${image}.lst3
${image}.psf


EXIT
END_DAOPHOT
#
ln1=`grep -n "Profile errors:" ${image}.psf.log | sed 's/:/ /g' | awk '(NR==1){print $1+2}'`
ln2=`grep -n "File with PSF stars and neighbors =" ${image}.psf.log | sed 's/:/ /g' | awk '(NR==1){print $1-4}'`
nl=$((ln2-ln1+1))
tail -n +${ln1} ${image}.psf.log | head -n ${nl} | sed -e 's/saturated/ 1.000   /g' -e 's/defective/ 1.000   /g' > ${image}.lst3.chi
nbad=`grep -e '?' -e '*' ${image}.lst3.chi | wc -l`
done
#
# End of WHILE loop
#
#######
#
mv ${image}.lst3 ${image}.plst
mv ${image}.lst3.chi ${image}.plst.chi
#rm ${image}.lst1 ${image}.lst2 ${image}.lst1.chi ${image}.lst2.chi
#
######################################################
###                                                ###
### 4th step : perform PSF photometry with ALLSTAR ###
###                                                ###
######################################################
#
echo "   Running 4th step ...."
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
echo "== Step 4 : Run ALLSTAR. ==" >> ${image}.log
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
#
##################################################
###                                            ###
### 5th step : get readied for running DAOGROW ###
###                                            ###
##################################################
#
echo "   Running 5th step ...."
rm ${image}a.log   >& /dev/null
rm ${image}a.ap    >& /dev/null
rm ${image}a.als   >& /dev/null
rm ${image}as.fits >& /dev/null
rm ${image}a.als.inp >& /dev/null
#rm allstar.opt     >& /dev/null
#
#
# If allstar.opt does NOT exist, create it
#
if [ ! -s allstar.opt ]; then
  cp ${image}.als.opt allstar.opt
fi
#
echo "" >> ${image}.log
echo "===============================================" >> ${image}.log
echo "== Step 5 : Prepare for aperture correction. ==" >> ${image}.log
echo "===============================================" >> ${image}.log
echo "" >> ${image}.log
#
#  Run DAOPHOT on the image with neighbor-subtracted PSF stars.
#
daophot << END_DAOPHOT >> ${image}a.log
OPTIONS
${image}.opt

ATTACH ${image}a.fits
PHOTOMETRY
apcor.opt

${image}.plst
${image}a.ap
exit

exit
END_DAOPHOT
#cat ${image}.als.opt > allstar.opt
cat ${image}.als.opt > ${image}a.als.inp
echo ${image}a.fits >> ${image}a.als.inp
echo ${image}.psf   >> ${image}a.als.inp
echo ${image}a.ap   >> ${image}a.als.inp
echo ${image}a.als  >> ${image}a.als.inp
echo ${image}as.fits  >> ${image}a.als.inp
echo exit  >> ${image}a.als.inp
allstar < ${image}a.als.inp >> ${image}a.log
#rm allstar.inp
rm ${image}as.fits >& /dev/null
# Delete temporarily uncompressed fits file
if [ ${fpackfile} == 1 ]; then
   echo "Removing temporarily uncompressed ${image}.fits file"
   rm ${image}.fits >& /dev/null
fi
# Fpack compress s.fits file
if [ -s ${image}s.fits ]; then
  # Not sure of the best way to test if fpack is available
  #  just try it and see if it works
  rm ${image}s.fits.fz >& /dev/null
  fpack ${image}s.fits >& /dev/null
  # It worked!
  if [ -s ${image}s.fits.fz ]; then
    echo "Fpack compressing ${image}s.fits"
    rm ${image}s.fits >& /dev/null
  fi
fi
###################################################
# Working in temporary directory, copy files back
if [ -n ${workdir} ]; then
   echo "Copying files back to original directory"
   # Delete some files
   rm photo.opt apcor.opt goodpsf.pro lstfilter  >& /dev/null
   rm ${image}.fits ${image}.fits.fz daophot.opt allstar.opt >& /dev/null
   rm ${image}.opt ${image}.als.opt >& /dev/null
   rm ${image}.cmn.coo ${image}.cmn.ap ${image}.cmn.lst ${image}.cmn.log >& /dev/null
   # Copy files back
   cp -f * ${curdir} >& /dev/null
   # Delete all temporary files
   rm -f * >& /dev/null
   # cd back
   cd ${curdir}
   # Delete temporary directory
   #  leave the base working directory
   rmdir ${workdir}
fi
echo ""

