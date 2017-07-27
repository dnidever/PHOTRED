#!/bin/bash

#################################################
## 
## This script performs the preprocess 
## stage to prepare the data for 
## crowding script
##
## AUTHOR: Antonio Dorta <adorta@iac.es>
## DATE: 2017-03-02
##
## Operations:
## 1) Read filenames to get the number of chips
## 2) Split the input stars file so that every
##    chip has a similar number of FULL lines
## 3) Create a directory per chip to store all 
##    needed files
## 4) Create a new MCH file with ".alf" files
##    instead of ".als". A link to this file
##    using underscores is made to avoid
##    problems with daomaster
## 5) Link/copy/move all specified files to
##    the right chip directory
## 
#################################################

# New extension for MCH with .alf 
ext=".alf.mch"
# Files to be copied
EXT_TO_COPY=(.fits .opt .alf .mch .psf .als .opt .als.opt .ap .raw .mag .log .weights .scale .zero _comb.psf _comb.opt _comb.als.opt _shift.mch .phot)
#EXT_TO_COPY=(alf)
FILES_TO_COPY=(apcor.lst extinction fields)

mkdir -p logs
INLIST=logs/ADDSTAR.inlist

#FILE_FIELDS=addstar_fields
#FILE_CHIPS=addstar_chips
starsDir=STARS
starsFn=input_stars

ARGS=( \
"1: whether link/copy/move data files from source to destination (or 'inlist' to only create input list)" \
"2: path to data files" \
"3: create (1) or not (anything else) a directory per field" \
"4: create (1) or not (anything else) a directory per chip" \
"5: path to input stars file" \
"6: split (1) or not (anything else) the input stars file among all chips" \
"7: file with info of chips" \
"8 (Optional): path to add as prefix if needed when creating links (only if arg6 is 'link')" \
)
# "8 (Optional): destination (current directory if none). It will be created if it does NOT exist" \


if [[ $# -ge 1 ]] && [[ $1 == 'inlist' ]]
then 
  # If there is only one argument and it is "inlist", only the list will be created
  find `pwd` -regextype sed -regex ".*F[0-9]*-[0-9]*_[0-9]*\.mch" | sort > $INLIST
  echo "Input list for ADDSTAR stage has been created in $INLIST. Skipping file transferring..."
  echo ""
  exit
elif [[ $# -eq 8 ]] 
then
  prefix=$7
elif [[ $# -eq 7 ]] 
then
  # If there is no destination directory, use the current one
  prefix=""
else
  echo -n "ERROR -> Syntax: $0 "
  count=1
  for args in "${ARGS[@]}"
  do
    echo -n "arg$count "
    (( count++ ))
  done
	echo ""
	echo "ARGUMENTS:"
	for arg in "${ARGS[@]}"
  do
    echo  " * arg$arg"
    (( count++ ))
  done

  >&2 echo "Syntax error when calling $0"	
  exit 1
fi

# Check that all needed arguments were specified
transferMethod=$1
orig=$2
sepFields=$3
sepChips=$4
starsOrig=$5
starsSplit=$6
chipsFile=$7
dest="."
 

# Determine whether copy|move|link the data files
case $transferMethod in
  "copy") cmd='cp -f'
  ;;
  "move") cmd='mv -f'
  ;;
  "link") cmd='ln -sf'
  ;;
*) cmd='cp -f'
  ;;
esac
if [[ $transferMethod == "move" ]]
then
  starsCmd='cp -f'
else
  starsCmd=$cmd
fi

# Check that star file exists
if [ ! -f $starsOrig ]; then
  >&2 echo "File $starsOrig NOT found"
  exit 2
fi

# Create DESTINATION
mkdir -p $dest

# Get CHIPS number based on filenames (skip all those missing chips)

FIELDS=(`find $orig -name "F*-*_??.fits" -printf "%f\n" | sed "s/^F//" | sed "s/-.*//" | sort -u`)
TOTAL_FIELDS=${#FIELDS[@]}

CHIPS=(`find $orig -name "F*-*_??.fits" -printf "%f\n" | sed "s/.*_//" | sed "s/\..*//" | sort -u`)
TOTAL_CHIPS=${#CHIPS[@]}

echo FIELDS: $FIELDS
echo CHIPS: $CHIPS

if [[ $TOTAL_FIELDS -eq 0 ]]
then
  >&2 echo "ERROR: No fields have been detected, check path to data files ($orig)"
  exit 3
fi
echo "Detected $TOTAL_FIELDS field(s)"


if [[ $TOTAL_CHIPS -eq 0 ]]
then
  >&2 echo "ERROR: No chips have been detected, check path to data files ($orig)"
  exit 4
fi
echo "Detected $TOTAL_CHIPS chip(s)"



starsDest="$dest/$starsDir"
mkdir -p $starsDest
if [[ $starsSplit -eq 1 ]]
then
  # SPLIT input file of stars into as many files as chips we have (with similar number of FULL lines!!)
  echo "Splitting $starsOrig in $TOTAL_CHIPS with similar number of lines (${starsFn}_XX.txt)"
  split -n l/${TOTAL_CHIPS} --additional-suffix=".stars" $starsOrig

  i=0
  for f in *.stars
  do
	  chip=${CHIPS[$i]}
    #mkdir -p $dest/chip_$chip
    #mv $f $dest/chip_$chip/${starsFn}_$chip.txt
    mv $f $starsDest/${starsFn}_${chip}.txt
    i=$(( i+1 ))
  done

else
  cp $starsOrig $starsDest/${starsFn}.txt
  starsFile=${starsFn}.txt
fi




# Check if origin directory is absolute or relative
if [[ "$orig" != /* && "$transferMethod" == "link"  ]]
then
  # Relative path, convert to absolute
	orig=`readlink -fv $orig`
  orig="$prefix$orig"
	echo "NEW ORIG: $orig"
fi


# Store current directory and move to destination
pushd .
cd $dest

# Create file with fields and chips
#rm -f $FILE_FIELDS
#for fld in ${FIELDS[@]}
#do
#	echo "F$fld" >> $FILE_FIELDS
#done

#rm -f $FILE_CHIPS
#for chp in ${CHIPS[@]}
#do
#	echo $chp >> $FILE_CHIPS
#done


# Create MCH with .alf files instead of .als
# Create symbolic links FX_...mch to FX-...mch 
#(daomaster does NOT work with '-')
echo "Creating MCH files"
for f in $orig/*\-*_??.mch
do 
	fn="`basename $f .mch`$ext"
	sed "s/\.als/\.alf/g" $f > $fn
	ln -s "$fn" "${fn//-/_}"
done

relDir="."
if [[ $sepFields -eq 1 ]]
then
  relDir="${relDir}/.."
fi
if [[ $sepChips -eq 1 ]]
then
  relDir="${relDir}/.."
fi

# Transfer common files
for file in ${FILES_TO_COPY[@]}
do
  $cmd $orig/${file} .
done


 
# Create a directory per EACH cheap and transfer related files to it
for fld in ${FIELDS[@]}
do
  for chp in ${CHIPS[@]}
  do
		destDir='.'
    if [[ $sepFields -eq 1 ]]
    then
      destDir="${destDir}/F${fld}"
    fi
    if [[ $sepChips -eq 1 ]]
    then
      destDir="${destDir}/chip${chp}"
    fi
    if [[ "$destDir" != "." ]]
    then
      mkdir -p $destDir
    fi
     

    echo " *** Create directory $destdir for Field F${fld} and Chip ${chp}. Use $cmd to transfer files"
    mv F${fld}*_${chp}.* $destDir

    # Transfer group of files by extension
	  for ext in ${EXT_TO_COPY[@]}
    do
      $cmd $orig/F${fld}-*_${chp}${ext} $destDir
    done



    # Get FI and PS values from OPT files of each chip
    grep "^FI" $orig/F${fld}*_${chp}.opt | awk '{print $NF}' > $destDir/F${fld}-chipFI_${chp}.dat
    grep "^PS" $orig/F${fld}*_${chp}.opt | awk '{print $NF}' > $destDir/F${fld}-chipPS_${chp}.dat

    # Transfer stars file

    cd $destDir
    if [[ $starsSplit -eq 1 ]]
    then
      starsFile="*${chp}.txt" 
    fi

    # Transfer stars
    $starsCmd $relDir/$starsDest/$starsFile F${fld}-stars_${chp}.txt

    # Link to common files
    for file in ${FILES_TO_COPY[@]}
    do
      ln -s $relDir/$file .
    done

    # Go back to base directory
    cd $relDir
      
  done
done

# Transfer CHIPS info file
echo "$starsCmd $orig/${chipsFile} ."
$starsCmd $orig/${chipsFile} .

# Transfer some other needed files
for file in ${FILES_TO_COPY[@]}
do
  $starsCmd $orig/$file .
done

# Create the INLIST file with all FX-NNNNNNN_YY.mch files
find `pwd` -regextype sed -regex ".*F[0-9]*-[0-9]*_[0-9]*\.mch" | sort > $INLIST
echo "Input list for ADDSTAR stage has been created in $INLIST"
 

# Return to initial directory
popd 
exit 0

#DONE!!!
