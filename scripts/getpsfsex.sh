#!/bin/sh
#
#############################################################################
#
# GETPSF.SH - Script to perform automated PSF photometry using DAOPHOT II
#                and Source Extractor
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
# Generate the SExtractor config file
echo "photred_mksexconfig,'${image}.fits','${image}.sex','${image}.cat'" | idl >& /dev/null
#
# Run SExtractor
sex ${image}.fits -c ${image}.sex >& /dev/null
#
# Convert SExtractor 
echo "sex2daophot,'${image}.cat','${image}.fits','${image}.coo'" | idl >& /dev/null
#
#  Run DAOPHOT to perform FIND, PHOTOMETRY, and PICKPSF on frame.
#  You may change the number of stars in PICKPSF.
#
daophot << END_DAOPHOT >> ${image}.log
OPTIONS
${image}.opt

ATTACH ${image}.fits
PHOTOMETRY


${image}.coo
${image}.ap
EXIT
END_DAOPHOT

# Get the magnitude limit using the 95% percentile mag - 1.5 mag
nl=`wc -l ${image}.coo | awk '{print int(0.95*$1)}'`
maglim=`awk '{if (NR>3) print $4}' ${image}.coo | sort | head -${nl} | tail -1 | awk '{print $1-1.5}'`
#maglim=`grep "Estimated magnitude limit" ${image}.log | awk '{print $6-1.0}'`

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
# Always run the LOCAL version of "lstfilter"
#./lstfilter << END_LSTFILTER >> ${image}.log
#${image}.lst
#${image}.coo
#${image}.lst1
#END_LSTFILTER
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
neifile=`grep "File with PSF stars and neighbors" ${image}.psf.log | awk '{print $8}'`
cp -f ${neifile} ${image}.nei >& /dev/null
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
rm ${image}.als.inp >& /dev/null
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
allstar < ${image}.als.inp >> ${image}.log
#rm allstar.inp >& /dev/null
#
