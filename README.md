# PHOTRED
Automated DAOPHOT-based PSF photometry pipeline

PHOTRED is an automated DAOPHOT-based PSF photometry pipeline.  It's very generic but does require reduced flat images.
PHOTRED performs WCS fitting, aperture photometry, single-image PSF photometry (ALLSTAR), source matching across multiple images,
forced PSF photometry across multiple exposures using a master source list created from a deep stack of all exposures
(ALLFRAME), aperture correction, calibration using photometric transformation equations, and dereddening.
STDRED is also included and can be used to derive the transformation equations from exposures of standard star fields.

PHOTRED does require the following that need to be installed separately:
- IRAF
- IDL
- The IDL Astronomy Uers's Library
- The stand-alone fortran version of DAOPHOT/ALLFRAME
- SExtractor

**This README file is now the main documentation.**
**Old** PDF manuals for PHOTRED and STDRED are included in the doc/ directory.  

IDL "sav" files of all the main PHOTRED programs are in the sav/ directory.  These can be used with the IDL virtual
machine for people who don't have an IDL license.

# Contents

* [Installation Instructions](#Installation_Instructions)
  * [1. Checkout PHOTREDL](#1_Checkout_PHOTRED)
  * [2. Download IDL Astro User's Library](#2_Download_Astro_Library)
  * [3. Make sure IDL/IRAF are available](#3_IDL_IRAF_Available)
  * [4. Make sure DAOPHOT/ALLFRAME is installed](#4_DAOPHOT_Installed)
  * [5. Schlegel Maps](#5_Schlegel_Maps)
  * [6. Setup Your IRAF Login File](#6_Setup_IRAF)
* [Running Instructions](#Running_Instructions)
  * [1. Data](#1_Data)
  * [2. Transformation Equations](#2_Transformation_Equations)
  * [3. Setup File](#3_Setup_File)
  * [4. Make sure your Imager is in the "imagers" file](#4_Check_Imager_File)
  * [5. Check "filters" file](#5_Check_Filters_File)
  * [6. Check "extinction" file](#6_Check_Extinction_File)
  * [7. Run PHOTRED_RENAME](#7_Run_PHOTRED_RENAME)
  * [8. Run PHOTRED](#8_Run_PHOTRED)		
* [Adding New Imagers](#Adding_New_Imagers)
* [Basic Explanation](#Basic_Explanation)
* [STAGES](#Stages)
  * [RENAME](#Stage_RENAME)
  * [SPLIT](#Stage_SPLIT)
  * [WCS](#Stage_WCS)
  * [DAOPHOT](#Stage_DAOPHOT)
  * [MATCH](#Stage_MATCH)
  * [ALLFRAME](#Stage_ALLFRAME)
  * [APCOR](#Stage_APCOR)
  * [ASTROM](#Stage_ASTROM)
  * [CALIB](#Stage_CALIB)
  * [COMBINE](#Stage_COMBINE)
  * [DEREDDEN](#Stage_DEREDDEN)
  * [SAVE](#Stage_SAVE)
  * [HTML](#Stage_HTML)
	

# <a name="Installation_Instructions"></a>Installation Instructions

## <a name="1_Checkout_PHOTRED"></a> 1. Checkout PHOTRED repository

Checkout this PHOTRED github repository.  Update your IDL_PATH in your
.cshrc, .tcshrc or .login file so that the "PHOTRED/pro/" directory is
included.  I would put PHOTRED near the beginning of the list so that
it's programs have priority if there happens to be other programs with
the same names.

You could also copy the PHOTRED/pro directory directly into your
"~/idl/" directory (maybe in a "photred/" subdirectory).  It's
probably better to keep full repository in a separate directly so that
it's easy to pull in updates.

There is one fortran code (lstfilter.f) that needs to be compiled.
Compile it like this:

```
gfortran lstfilter.f -o lstfilter
```

## <a name="2_Download_Astro_Library"></a> 2. Download IDL Astro User's Library

You will also need the IDL Astro User's Library.  Most major astronomy
institutions already have this installed but it's good to check that
you have an updated copy.  You can download a tar file with all the
programs ("astron.dir.tar.gz") from their ftp site:

http://idlastro.gsfc.nasa.gov/ftp/

**This is a bit outdated and I don't think it works anymore**.
There can be a problems if you have older copies of the IDL Astro
User's programs in your IDL directory or programs with the same
name. I've included a program called "checkidlastro.sh" in the
photred_idl.tar file that will print out the names of programs in your
~/idl/ directory (and subdirectories) that have the same name as
programs in the IDL Astro User's Library (run it by typing
"./checkidlastro.sh" in your ~/idl/ directory; it needs the file
"idlastro.lst" which is also included). If it finds anything I would
erase the offending program or rename it (e.g. proname.pro ->
proname.pro.bak).

## <a name="3_IDL_IRAF_Available"></a> 3. Make sure IDL/IRAF are available

PHOTRED needs IDL and IRAF to run. If you don't have IDL then you
might consider buying a license from Exelis Visual Information
Solutions (http://exelisvis.com/). Student licenses are available and
full licenses have come down in price. If these are too expensive,
then you might consider installing the free version of IDL, the GNU
Data Language (http://gnudatalanguage.sourceforge.net). I have not
tested PHOTRED on GDL, but I think GDL should have all of the
functionality needed to run PHOTRED. However, some PHOTRED programs
might need to be tweaked to work with GDL.   You can also download
the IDL Virtual Machine for free and run IDL "sav" files. I have
made IDL sav files for the PHOTRED pipeline so they can be run with
the IDL Virtual Machine.  They are in the "sav/" directory of this
repository.  To run a program type
```
idl -vm=progname.sav
```

IRAF is freely downloadable from http://iraf.net.

## <a name="4_DAOPHOT_Installed"></a> 4. Make sure DAOPHOT/ALLFRAME is installed

Make sure that DAOPHOT/ALLSTAR/ALLFRAME and SExtractor are installed. Type:

```
which daophot
which allstar
which daomaster
which daogrow
which allframe
which sex
```

They should all return the name of the program. If you get an error
then you will need to install that program.

Currently daophot and allstar are hardwired to /net/astro/bin/ in the
daophot.sh and getpsf.sh scripts. If this is not where your version of
daophot/allstar are located then change those lines or comment them
out.

The newest version of ALLFRAME had an error in it that caused it to
print out zero's for chi and sharp. The fixed version on the UVa Astro
machines is called "/net/halo/bin/allframe.2004.fixed". This is
currently hardcoded in photred_allframe.pro and allframe.pro. Make
sure that this program exists.

## <a name="5_Schlegel_Maps"></a> 5. Schlegel Maps

PHOTRED uses the Schlegel maps to deredden the photometry. The FITS
files can be downloaded from here. On the UVa Astro machines the files
are located here: /net/grass/catalogs/reddening/

You need to update your ".cshrc" file to set the Unix environment
variable "DUST_DIR" to the directory where the Schlegel dust maps are
located (actually one directory up in the directory tree). At UVa
Astro this is the line you should add to your ".cshrc" file:

```
setenv DUST_DIR /net/grass/catalogs/reddening/
```

Now test that it works (you must already have installed the PHOTRED
IDL files):

```
% idl
IDL>print,dust_getval(10,10)
 ￼￼0.472163
```

If you get an error here, then there is a problem. Check that all the
files and required programs are there.

## <a name="6_Setup_IRAF"></a> 6. Setup Your IRAF Login File

PHOTRED calls some IRAF programs. In order for this to work you need
to edit your IRAF "login.cl" file so that it doesn't print out any
messages to the screen. Comment out the 9 lines following "# Set the
terminal type." near the top of the file, and possibly also the 4
lines following "# Delete any old MTIO lock (magtape position) files."
(these give problems on the Pleione cluster). If that still doesn't
work, create a blank file called ".hushiraf" in your iraf directory.

```
touch .hushiraf
```

Start IRAF by typing "cl" in your IRAF directory and see what
happens. Nothing should be written to the screen except "cl>" or maybe
"ecl>". If it's still printing other things to the screen, then you'll
need to comment out more lines from the "login.cl" file.

If you have a non-standard IRAF installation and you need to specify
an absolute path (and possibly additional optional flags), you'll
need to make an alias for "cl" (that would be in your shell start-up file),
e.g.
```
alias cl "/absolute/path/to/cl/command/cl -plusoptions"
```

# <a name="Running_Instructions"></a>Running Instructions

## <a name="1_Data"></a> 1. Data

Start by putting all the final flat frames in their nights
directory. Process them with IRAF's CCDRED, MSCRED or Armin Rest's
SuperMacho pipeline to get final flat images. PHOTRED only runs on one
night's worth of data. If you have 5 nights of data then you will need
to run PHOTRED 5 times (they can be running simultaneously in separate
directories).

FIX BAD PIXELS!!!. Make sure to fix any bad pixels or columns in the
images because otherwise they might get detected as "sources" and mess
up WCS, DAOPHOT, MATCH, and ALLFRAME.

## <a name="2_Transformation_Equations"></a> 2. Transformation Equations

In order to do calibration step (CALIB) PHOTRED needs the photometry
transformation equations which means doing the standard star
reduction. There is a separate pipeline for the standard star
reduction called STDRED that you should run first to get the
transformation equations for the entire run. You can run all of the
PHOTRED stages before CALIB (which include the cpu intensive DAOPHOT
and ALLFRAME stages) without the transformation equations, and then do
the final stages once you have transformation equations.

The transformation file needs to be in this format. Each filter gets
two lines. The first line should contain: filter, color, nightly
zero-point term, airmass term, color term, airmass*color term, color^2
term. The second line gives the uncertainties in the five
terms. Normally the last two terms (airmass*color and color^2 are
0.0000).

n1.trans

```
M M-T  -0.9990 0.1402 -0.1345 0.0000 0.0000
       1.094E-02 5.037E-03 2.010E-03 0.0000 0.0000
T M-T   -0.0061 0.0489 0.0266 0.0000 0.0000
       6.782E-03 3.387E-03 1.374E-03 0.0000 0.0000
D M-D   1.3251 0.1403 -0.0147 0.0000 0.0000
       1.001E-02 5.472E-03 2.653E-02 0.0000 0.0000
```

## <a name="3_Setup_File"></a> 3. Setup File

PHOTRED needs a "photred.setup" file to run. This file specifies a few
important parameters. Here's an example of a "photred.setup" file. The
various parameters are described below.

```
##### REQUIRED #####
scriptsdir  /idl/idllocal/photred/PHOTRED/scripts/
irafdir     /home/smash/
telescope   Blanco
instrument  DECAM
observatory CTIO
nmulti      10
filtref     g,i,r,z,u
trans       nocalib.trans
##### OPTIONAL #####
sepfielddir  1
keepmef      0
redo         0
#skipwcs     0
#wcsup       N
#wcsleft     E
#pixscale    0.50
#wcsrefname  USNO-B1
#searchdist  60
#wcsrmslim   1.0
#hyperthread  1
daopsfva      2
daofitradfwhm 1.0
psfcomsrc     1
psfcomglobal  1
#psfcomgauss  1
#mchmaxshift  50.0
finditer      2
#alfdetprog  sextractor
#alfnocmbimscale 0
#alfexclude  F1,F3
ddo51radoffset  1
keepinstr   1
avgmag      1
avgonlymag  0
todered     M,T,D,M-T,M-D
#toextadd    M,T,D,M-T,M-D
#cmd2cdaxes  M,M-T,M-D
stdfile	     /home/smash/SMASH/standards.txt
##### STAGES #####
rename
split
wcs
daophot
match
allframe
apcor
astrom
calib
combine
deredden
save
html
```

Parameter  |  Description
------------ | -------------
**REQUIRED** |  **REQUIRED**
scriptsdir |  The absolute path to the directory that contains the PHOTRED scripts (i.e. daophot.sh, etc.)
irafdir | The absolute path to your IRAF directory (that contains your login.cl file)
telescope | The name of the telescope (e.g. Blanco, Swope)
instrument | The name of the instrument (e.g. MOSAIC)
observatory | (OPTIONAL) The name of the observatory. This is needed if the header does not contain the AIRMASS and it needs to be calculated from the date, ra/dec and observatory location.
nmulti |  The number of processors PHOTRED should use for the DAOPHOT and ALLFRAME stages (only relevant if on the Pleione cluster).
filtref |  The shortname of the filter (specified in the "filters" file) to be used as the reference frame (e.g. the M frame). This can be a comma-delimited priority-ordered list of filters (e.g. g,i,r,z,u). If there are multiple observations in this filter then the longest exposure in this filter will be used.
trans | The name of the file that contains the photometric transformation equations.
keepmef | OPTIONAL. Multi-extension files (MEF) are split by PHOTRED. Do you want PHOTRED to keep the MEF files: YES=1, NO=0 (i.e. erase them).
**OPTIONAL**  |  **OPTIONAL**
sepfielddir |  Put each field in a separate directory (this is now the default option), otherwise everything will go in the main directory and can slow down processing because a very large (~100,000) number of fileds.
keepmef  |  Keep the original multi-exension files (MEF).
redo | PHOTRED will NOT reprocess files that have already been processed unless "redo" is set. This can also be set as a keyword on the command line (i.e. IDL>photred,/redo).
skipwcs | Set this if your images already have correct WCS in their headers and you don't want the WCS to be refit in the WCS stage.
wcsup | What cardinal direction (i.e. N, S, E or W) is "up" in the image? This is only used for non-standard setups.
wcsleft | What cardinal direction (i.e. N, S, E or W) is "left" in the image? This is only used for non-standard setups.
pixscale | The plate scale in arcseconds/pixel. This is ONLY used for non-"standard" imagers (i.e. not MOSAIC, IMACS, LBC or Swope) where the pixel scale cannot be determined from the image headers.
wcsrefname | The name of the WCS reference catalog to use. The two options are 'USNO-B1' and '2MASS-PSC'. USNO-B1 is the default. The astrometric accuracy of the 2MASS catalog is better (~0.170 arcsec) than USNO-B1 (~0.270 arcsec), but it does not go as deep (R~18) as USNO-B1 (R~20). So if you have deep images then definitely use USNO-B1, but if you have moderately deep images then 2MASS-PSC is probably better (and faster).
searchdist | This sets the search distance (in arcmin) for WCS fitting (PHOTRED_WCS). Normally this is not needed. The default is 2*image size > 60 arcmin (i.e. whichever is greater). This is normally sufficient. If the WCS isn't fitting correctly then try setting "searchdist" to a larger value.
wcsrmslim | This is the maximum RMS (in arcseconds) allowed for an acceptable WCS fit. The default is 1.0 arcseconds. Normally the RMS values are ~0.2-0.3 arcseconds.
hyperthread | This allows multiple jobs to be running (daophot and allframe only) on a computer (such as halo or stream) that has multiple processors. It's similar to running it on a cluster.
wcscaterrlim | This is an error limit on the detected sources in an image used to match to an astrometric reference catalog.
daopsfva | The DAOPHOT .opt VA option, which is for the spatially-varying analytical PSF (0: constant, 1: linear, 2: quadratic).  The default internal value is 2.
daofitradfwhm | The DAOPHOT .opt file PSF Fitting Radius to use in terms of the FWHM.  1.0 x FWHM is the default option.  Smaller values are often better in crowded regions.
psfcomsrc | If this is set to 1 then only sources detected in other frames of the same field are used as PSF stars. This makes sure that your PSF stars are not contaminated by cosmic rays or other junk. HIGHLY RECOMMENDED.
psfcomglobal | This finds PSF stars in other exposures of the same field in a "global" manner that works better with large dithers.  HIGHLY RECOMMENDED.
psfcomgauss | This fits Gaussians to each potential PSF star and requires the Gaussian parameters to be "reasonable".
mchmaxshift | This sets a maximum constraint on the X/Y shifts found in the MATCH stage. WARNING, only use this if you ABSOLUTELY known that your frames are well aligned already and do not have large dithers.
finditer | The number of times to iteratively find sources in ALLFRAME (allfprep). The default is 2.
alfdetprog | The program to use for source detection in the ALLFRAME stage (allfprep). The options are "sextractor" and "daophot". The default is "sextractor". SExtractor is generally better at finding faint sources and returns a stellaricity probability value which is very useful. HOWEVER, SExtractor fails in VERY crowded regions. It's best to use DAOPHOT for very crowded images.
alfnocmbimscale | Do not "scale" the combined images in the ALLFRAME prep stage.
alfexclude  |  Comma-delimited list of fields (e.g., F1, F3, F5) to exclude from ALLFRAME processing.
alfusecmn  |  Use the reference image common sources file to pick PSF stars for the combined image.
alftiletype  |  The type of combination method to use.  The old method is "ORIG" while the new and improved method is "WCS".
ddo51radoffset  |  There is a photometric offset in the DDO51 filter that depends on the radial distance from the center of the field. Currently this is only observed in the CTIO+MOSAIC data. Setting this parameter will remove this offset (done in CALIB). If you use this make sure to also use it in STDRED.
keepinstr | CALIB should keep the instrumental magnitudes in the final output file.
avgmag | CALIB should calculate average magnitudes in filters that were observed multiple times. The individual magnitudes are also kept.
avgonlymag | Same as "avgmag" but only keeps the average magnitudes.
todered | The magnitudes and colors the DEREDDEN stage should deredden. The short names of the filters and colors should be used (i.e. M, T, M-T). The dereddened magnitudes and colors will have the same names but with a "0" (zero) after it (i.e. M0 for M). The dash is removed for the color names.
toextadd | The extinction and reddening values to add. The short names of the filters and colors should be used (i.e. M, T, M-T). The extinctions will have the same name as the magnitude but with an "A" at the beginning (i.e. AM for M), and the reddening values will have the same name as the color names but with the dash removed and an "E" at the beginning (i.e. EMT for M-T).
cmd2cdaxes | The magnitudes and colors to use for the CMD and 2CD plots in PHOTRED_HTML.PRO. It should be in this format: MAG, MAG1-MAG2, MAG3-MAG4, i.e. M,M-T,M-D. For this example the CMD would be M vs. M-T and the 2CD would be M-D vs. M-T. The CMD and 2CD will be connected (CMD on top, 2CD on bottom) and will share the x-axis which be the first color. The magnitude will also be used for other plots (i.e. error vs. magnitude, chi vs. magnitude, etc).
stdfile  |  List of standard star fields to exclude from processing and move to the "standards/" directory.

Then add all the stage names that you want to process from the
following list: rename, split, wcs, daophot, match, allframe, apcor,
calib, astrom, combine, deredden, and save

The number of ALLFRAME iterations can be set in the "allframe.opt" file, the MA parameter. The default is 50.

## <a name="4_Check_Imager_File"></a> 4. Make sure your Imager is in the "imagers" file

There is an "imagers" file in your scripts directory. If the imager
you are using is not in the list then add it at the end. You need the
telescope name, instrument names, observaotry name (as found in the
OBSERVATORY.PRO IDL program), the number of amplifiers, and the
separator character (only necessary for multi-amplifier imagers). The
TELESCOPE and INSTRUMENT names in your "photred.setup" file needs to
match the ones in the "imagers" file.

```
# TELESCOPE   INSTRUMENT  OBSERVATORY  NAMPS  SEPARATOR
blanco        decam       CTIO         62     _
blanco        mosaic      CTIO         16     _
mayall        mosaic      KPNO         8      _
baade         imacs       LCO          8      c
lbt           lbc         MGIO         4      _
swope         ccd         LCO          1
rrrt          sbig        FMO          1      -
2p2           wfi         LASILLA      8      _
```

The case is not important. The separator character is used to separate
the amplifier number from the main frame name for multi-amplifier
images (i.e. "obj1034_1.fits" is the filename of the first amplifier
FITS file of the "obj1034" image). The separator character is almost
always the underscore "_". The only exception (for now) is IMACS which
uses a "c". If your images are in the multi-extension format (MEF)
then they will be split up into amplifier files in the SPLIT stage of
PHOTRED with the IRAF task MSCSPLIT which uses the underscore as the
separator character.

## <a name="5_Check_Filters_File"></a> 5. Check "filters" file

PHOTRED uses short names for filters, and these are stored in the
"filters" file. One is provided in the scripts tar file. This is what
it looks like:

```
'M Washington c6007'    M
'M Washington k1007'    M
'M'                     M
'I c6028'               T
'I Nearly-Mould k1005'  T
'T'                     T
'T2'                    T
'DDO 51 c6008'          D
'D51 DDO c6008'         D
'D51 DDO 51 c6008'      D
'D51 DDO 51 k1008'      D
'D'                     D
'D51'                   D
'DDO51'                 D
'DDO-51'                D
'DDO_51'                D
'DDO 51'                D
'51'                    D
'T1'                    T1
'C'                     C
'B Harris c6002'        B
'B-BESSEL'              B
'B'                     B
'V'                     V
'R-BESSEL'              R
'R'                     R
'u DECam c0006 3500.0 1000.0'         u
'g DECam SDSS c0001 4720.0 1520.0'    g
'r DECam SDSS c0002 6415.0 1480.0'    r
'i DECam SDSS c0003 7835.0 1470.0'    i
'z DECam SDSS c0004 9260.0 1520.0'    z
'Y DECam c0005 10095.0 1130.0'        Y
'u'                     u
'g'                     g
'r'                     r
'i'                     i
'z'                     z
'Y'                     Y
```

The first column is the text that is found in the FILTER keyword in
the FITS header. The second column is the shortname. These can be
repeated because different observatories have different names for the
same filter. These shortname will be used in the transformation file
and the "extinction" file.

Make sure that your filters appears in the list (leading/trailing
spaces are not important. If they don't then PHOTRED will make a new
entry in the "filters" file for your filter and make a new
shortname. This is NOT desirable because the shortname probaby won't
match what you have in your transformation file or in the "extinction"
file.

Here's a way to double-check if your filter is in the "filters"
file. Copy the "filters" file from your scriptsdir to the directory
where your data is. Then run PHOTRED_GETFILTER on one FITS file per
filter. Change the shortname to a single capital letter if possible.

```
IDL>print,photred_getfilter('ccd1001.fits')
D
```

PHOTRED_GETFILTER will print out the shortname of the filter. It will
tell you if it didn't find the filter in the "filters" file and what
new entry it added. It's preferable to have the most up to date
"filters" file in the scriptsdir directory so that it can be used for
the next run. PHOTRED uses the "filters" file in the main directory
(where "photred.setup" and the data are located) if there is one,
otherwise it will copy the "filters" file from the scriptsdir
directory.

## <a name="6_Check_Extinction_File"></a> 6. Check "extinction" file

If you want your photometry dereddened then PHOTRED needs to know what
extinction value to use for each filter. This is stored in the
"extinction" file. The first column has the filter shortname and the
second column A(filter)/E(B-V) for that filter. The reddenings for
colors can be derived from these values, so no entries for colors are
needed.

```
# Filter    A(filter)/E(B-V)
M           3.43
T           1.83
D           3.37
u           4.239
g           3.303
r           2.285
i           1.698
z           1.263
```

Make sure your filter and its appropriate extinction value appears in
this file. Try to update the "extinction" As with the "filters" file
PHOTRED uses the "extinction" file in the data directory if there is
one, otherwise it will copy the "extinction" file from the scriptsdir
directory.


## <a name="7_Run_PHOTRED_RENAME"></a> 7. Run PHOTRED_RENAME

Okay, now you're ready to PHOTRED. PHOTRED has 13 stages and there is
a separate IDL program for each stage (e.g. PHOTRED_DAOPHOT). Each
stage can be run on it's own. The PHOTRED program is actually just a
giant wrapper for the 13 stages (and the ones specified in the
"photred.setup" file) run in the correct order. If you ever want to
just run ONE stage then it's probably easier to just run that stage at
the command line instead of editing the "photred.setup" file and
running photred.

It's preferable to run the very first stage, PHOTRED_RENAME, by itself
from the command line and double-check the results. This stage
prepends a string to each FITS filename that indicates what field it
is (i.e. obj1101.fits -> F1-obj1101.fits), and it creates a "fields"
file that has the actual field name for each field "shortname". This
is information that PHOTRED_MATCH needs to match up the various ALS
files for each field/chip group. It's very important that the files
are renamed properly. So check closely the text that is
output. PHOTRED_RENAME uses the first "word" in the OBJECT FITS
keyword as the field name. Any
zero/flat/twilight/sky/pointing/focus/test frames and standard star
frames (SA98, SA110, SA114 and NGC3680) are moved to a "calib/"
directory.

It's a good idea to run PHOTRED_RENAME in "testing" mode, so you can
see how it will rename files without it actually doing anything. Just
type "photred_rename,/testing". The first thing PHOTRED_RENAME does is
check that all of the FITS header parameters can be found (readnoise,
gain, ut-time, filter, exposure time, ra, dec, date, and airmass). If
any of these cannot be found then it will spit out errors. Watch for
these! Check that the actual field names (not the "shortnames") in
"fields" are correct. This information is not used until PHOTRED_SAVE
renames the final photometry files. **Make sure to change the "fields"
file ONLY AFTER running PHOTRED_RENAME in the NORMAL mode.** In testing
mode PHOTRED_RENAME will write the field information to
"fields.testing" instead of "fields".

If the files were not renamed properly, rename them by hand and update
the "fields" file. Also, update the "logs/RENAME.outlist" file. It
might be easiest to delete the "logs/RENAME.outlist" file and remake
it by typing
```
% ls F*.fits >> logs/RENAME.outlist
```

## <a name="8_Run_PHOTRED"></a> 8. Run PHOTRED

Now run PHOTRED. Start idl, type "photred" and you're off!!  You can
also run PHOTRED in the background. Make a batch file called
"photred.batch" that has a single line with "photred". You can then
run this batch file with idlbatch or idlbatchn (a "niced"
version). The "idlbatch" programs will run the IDL job in the
background and create a log file called "photred.batch.log". Make sure
to put the "idlbatch" programs in your ~/bin/ directory and that
~/bin/ is also in your path (check your .cshrc file).

Double check (in the logfile) that the WCS is being fit correctly. The
Total RMS should be around 0.2-0.4 arcsec. If it's not working
properly check that the images have RA/DEC or CRVAL1/CRVAL2 in them
and that the pixel scale is correct. If it still isn't working then
you can try setting "searchdist" larger (the default is 2*image size >
60 arcmin). You can also set the maximum acceptable RMS with
"wcsrmslim" in the "photred.setup" file.

You can use PHOTRED_SUMMARY (run it in the same directory) to get an
update on the progress of PHOTRED. Most of the stages update their
respective lists and so the number of files in the
inlist/outlist/success/failure lists will tell you what's going
on. However, the DAOPHOT and ALLFRAME stages update their lists after
all files have finished, so for these two stages PHOTRED_SUMMARY
checks to see which files have the expected output (i.e. .als or .mag
files). This information is listed in the "COMPLETED" column. If you
are redoing some files then these numbers won't be accurate.

## <a name="Adding_New_Imagers"> Adding new imagers

PHOTRED currently works on data from KPNO+MOSAIC, CTIO+MOSAIC
(Blanco), Swope CCD, IMACS, LBT Camera (LBC), and DECam. New imagers
can be added, but there are a couple of things that need to be
double-checked:
* Make sure to add your imager to the "imagers" file in the scripts directory.
* Make sure that PHOTRED_GETUTTIME.PRO, PHOTRED_GETFILTER,
PHOTRED_GETEXPTIME.PRO, PHOTRED_GETGAIN.PRO, PHOTRED_GETRDNOISE.PRO,
PHOTRED_GETDATE.PRO and PHOTRED_GETAIRMASS.PRO return the proper
values. PHOTRED_RENAME.PRO checks that all of the appropriate keywords
are in the FITS headers. Run PHOTRED_RENAME in /testing mode and see
if you get any errors. If there are errors then the above programs
might need to be modified to deal with the new data type.
* Add the necessary filters to the "filters" file
* Make sure that PHOTRED_WCS.PRO can properly process the
images. PIXSCALE might need to be specified in the "photred.setup"
file.

## <a name="Basic_Explanation"></a> Basic Explanation

Photred is meant to be run on one night's data at a time. The standard
star reduction should already have been done and the transformation
equation put in the directory.

The pipeline is split into stages and the files are "shuttled" from
stage to stage via lists of files. Each stage has an INLIST and
OUTLIST. The INLIST is the list of files to process, and the OUTLIST
is the list of files output. Normally the INLIST of files is moved
over from the OUTLIST of the previous stage. The INLIST files that are
successfully processed are removed from the INLIST file, and are added
to the SUCCESS list. INLIST files that are NOT successfully process
are left in the INLIST file and are added to the FAILURE list.

Each stage has several log files associated with it:
- **INLIST** The list of files to process. These are normally moved over
from the OUTLIST of the previous stage.
- **OUTLIST** The files successfully output from the stage. These might be
in a different format from the INLIST files.
- **SUCCESS** The files in INLIST that were successfully processed
- **FAILURE** The files in INLIST that were NOT successfully processed LOG A
running log of what the stage has done

The main interface and logistical work of PHOTRED is done in IDL. Most
of the heavy processing is done by Peter Stetson's photometric codes
DAOPHOT, ALLSTAR, DAOMATCH, DAOMASTER, ALLFRAME and DAOGROW. Some IRAF
tasks, such as MSCCMATCH, are also used.

PHOTRED has the built-in capability to run multiple DAOPHOT or
ALLFRAME processes at the same time on Pleione (i.e. our beowulf
cluster). You have to be logged in to Pleione and be able to access
the data in order to take advantage of the multi-plexing
capability. The maximum number of processes that will be allowed to be
in the PBS queue at any given time is set by the NMULTI parameter in
the photred.setup file (NMULTI=8 is a good value to use).

# <a name="Stages"></a> STAGES

##  <a name="Stage_RENAME"></a> RENAME

### Basic Explanation

This program renames object fits files so that it includes their field
information. For example, ccd1001.fits gets renamed to
F1.ccd1001.fits. Any calibration frames (zero, dflat, sflat, etc.) get
moved to the "calib/" directory without getting renamed.

### Lists

It will create the inlist from all fits files in the directory Outlist
will be of all files that aren't zero, dflat, sflat, etc. and were
successfully put in a "field". It also creates the "fields"
file. Files are renamed: F1.obj1023.fits, F2.obj1045.fits, etc.

**INLIST** (fits) - Creates it itself from fits files in directory

Single-Chip | Split Multi-chip | Multi-chip (MEF)
------------ | ------------- | -------------
zero1001.fits | zero1001_1.fits | zero1001.fits
ccd1012.fits | ccd1012_2.fits | ccd1012.fits
ccd1024.fits | ccd1024c3.fits | ccd1024.fits
ccd1053.fits | ccd1053c5.fits | ccd1053.fits

**OUTLIST** (fits) - The renamed object files

Single-Chip | Split Multi-chip | Multi-chip (MEF)
------------ | ------------- | -------------
F1.ccd1012.fits | F1.ccd1012_2.fits | F1.ccd1012.fits
F2.ccd1024.fits | F2.ccd1024c3.fits | F2.ccd1024.fits
F3.ccd1053.fits | F3.ccd1053c5.fits | F3.ccd1053.fits


##  <a name="Stage_SPLIT"></a> SPLIT

### Basic Explanation

This splits multiple-extension files (MEF) into separate files for
each amp/chip. Non-MEF files are not affected.

### Lists

All the files in rename.outlist are put in split.inlist. All non-MEF
(single chip) files automatically go to split.outlist. All
successfully split files (not original MEF files) are put in
split.outlist

**INLIST** (fits) - Copied from rename.outlist

Single-Chip | Split Multi-chip | Multi-chip (MEF)
------------ | ------------- | -------------
F1.ccd1012.fits | F1.ccd1012_2.fits | F1.ccd1012.fits
F2.ccd1024.fits | F2.ccd1024c3.fits | F2.ccd1024.fits
F3.ccd1053.fits | F3.ccd1053c5.fits | F3.ccd1053.fits

**OUTLIST** (fits) - All split files that are split okay, or single-chip files

Single-Chip |  Split Multi-chip | Multi-chip (MEF)
------------ | ------------- | -------------
F1.ccd1012.fits | F1.ccd1012_2.fits | F1.ccd1012_1.fits, F1.ccd1012_2.fits, ...
F2.ccd1024.fits | F2.ccd1024c3.fits | F2.ccd1024_1.fits, F2.ccd1024_2.fits, ...
F3.ccd1053.fits | F3.ccd1053c5.fits | F3.ccd1053_1.fits, F3.ccd1053_2.fits, ...

From now on the "Split Multi-chip" and "Multi-chip" files will "look"
the same, since now the MEF files have been split.


## <a name="Stage_WCS"></a> WCS

### Basic Explanation

This program gets the correct WCS for images.

### Lists

It takes all files from split.outlist and puts them into
wcs.inlist. All files that succeeded get put in wcs.outlist.

**INLIST** (fits) - Moved from split.outlist

Single-Chip | Split Multi-chip | Multi-chip (MEF)
------------ | ------------- | -------------
F1.ccd1012.fits | F1.ccd1012_2.fits | F1.ccd1012.fits
F2.ccd1024.fits | F2.ccd1024c3.fits | F2.ccd1024.fits
F3.ccd1053.fits | F3.ccd1053c5.fits | F3.ccd1053.fits

**OUTLIST** (fits) - All object files that are given a proper wcs.

Single-Chip |  Split Multi-chip | Multi-chip (MEF)
------------ | ------------- | -------------
F1.ccd1012.fits | F1.ccd1012_2.fits | F1.ccd1012.fits
F2.ccd1024.fits | F2.ccd1024c3.fits | F2.ccd1024.fits
F3.ccd1053.fits | F3.ccd1053c5.fits | F3.ccd1053.fits


## <a name="Stage_DAOPHOT"></a> DAOPHOT

### Basic Explanation

This program runs DAOPHOT on all the images (using Tony's script).

### Lists

If split.outlist exists then they are put into daophot.inlist. If
split.outlist does NOT exist then the wcs.outlist is taken
instead. All fits files that successfully run through daophot get put
into daophot.outlist

**INLIST** (fits) - Moved from split.outlist

Single-Chip | Split Multi-chip
------------ | -------------
F1.ccd1012.fits | F1.ccd1012_2.fits
F2.ccd1024.fits | F2.ccd1024c3.fits
F3.ccd1053.fits | F3.ccd1053c5.fits

**OUTLIST** (als) - All files in inlist that are succcessfully processed

Single-Chip | Split Multi-chip
------------ | -------------
F1.ccd1012.als | F1.ccd1012_2.als
F2.ccd1024.als | F2.ccd1024c3.als
F3.ccd1053.als | F3.ccd1053c5.als


## <a name="Stage_MATCH"></a> MATCH

### Basic Explanation

This program runs DAOMATCH and DAOMASTER on the ALS files which combines the photometry from the various filters.

### Lists

It takes all the files from daophot.outlist and puts them into match.inlist. It uses the fields (from the filenames) to figure out which files to together and should be matched. All of the MCH files to in the match.outlist.

**INLIST** (als) - Moved from daophot.outlist

Single-Chip | Split Multi-chip
------------ | -------------
F1.ccd1012.als | F1.ccd1012_2.als
F2.ccd1024.als | F2.ccd1024c3.als
F3.ccd1053.als | F3.ccd1053c5.als

**OUTLIST** (mch) - All files in inlist that were successfully matched (ALL??)

Single-Chip | Split Multi-chip
------------ | -------------
F1.ccd1012.mch | F1.ccd1012_2.mch
F2.ccd1024.mch | F2.ccd1024c3.mch
F3.ccd1053.mch | F3.ccd1053c5.mch

There won't be as many MCH files in the outlist as ALS files in the
inlist. If there are 3 frames per field then there will be 3x as many
ALS files as MCH files. The MCH files will have the names of the
"reference" frame.

If any of the ALS files don't match then the entire "set" of images fails.


## <a name="Stage_ALLFRAME"></a> ALLFRAME

### Basic Explanation

This program runs ALLFRAME on all of the MCH files. ALLFRAME does PSF fitting on images of all filters/bands similtaneously.

### Lists

list of all mch files Takes the list of MCH files from the
match.outlist and puts them in allframe.inlist. All fields that
succeed, their MAG files get put into allframe.outlist.

**INLIST** (mch) - Moved from match.outlist

Single-Chip | Split Multi-chip
------------ | -------------
F1.ccd1012.mch | F1.ccd1012_2.mch
F2.ccd1024.mch | F2.ccd1024c3.mch
F3.ccd1053.mch | F3.ccd1053c5.mch

**OUTLIST** (mag) - Every MCH file that was processed successfully by
  allframe and has a MAG file.

Single-Chip | Split Multi-chip
------------ | -------------
F1.ccd1012.mag | F1.ccd1012_2.mag
F2.ccd1024.mag | F2.ccd1024c3.mag
F3.ccd1053.mag | F3.ccd1053c5.mag


## <a name="Stage_APCOR"></a> APCOR

### Basic Explanation

This programs find the aperture correction for all the files using DAOGROW.

### Lists

Takes all of the fits files from daophot.success(!!) and puts them
into apcor.inlist. All fits files that have an aperture correction in
the final apcor.lst get put into the apcor.outlist.

**INLIST** (fits) - COPIED from daophot.success

Single-Chip | Split Multi-chip
------------ | -------------
F1.ccd1012.fits | F1.ccd1012_2.fits
F2.ccd1024.fits | F2.ccd1024c3.fits
F3.ccd1053.fits | F3.ccd1053c5.fits

**OUTLIST** (fits) - Every FITS file that was successfully given an
aperture correction in the final apcor.lst file

Single-Chip | Split Multi-chip
------------ | -------------
F1.ccd1012.fits | F1.ccd1012_2.fits
F2.ccd1024.fits | F2.ccd1024c3.fits
F3.ccd1053.fits | F3.ccd1053c5.fits

**SUCCESS** (fits) - Same as outlist.


## <a name="Stage_ASTROM"></a> ASTROM

### Basic Explanation

This program gets coordinates for all stars from the WCS in the
reference image.

### Lists

This takes all of the mag files from ALLFRAME.outlist or mch files
from MATCH.outlist. All files that are successfully given coordinates
are put into the ASTROM.outlist (with .ast endings). There will be a
separate file for each field chip/amp.

**INLIST** (mag/mch) - Moved from ALLFRAME.outlist or if that does not
exist, then from MATCH.outlist

Single-Chip | Split Multi-chip
------------ | -------------
F1.ccd1012.mag | F1.ccd1012_2.mag
F2.ccd1024.mag | F2.ccd1024c3.mag
F3.ccd1053.mag | F3.ccd1053c5.mag

**OUTLIST** (ast) - Every MAG/MCH file that was processed successfully

Single-Chip | Split Multi-chip
------------ | -------------
F1.ccd1012.ast | F1.ccd1012_2.ast
F2.ccd1024.ast | F2.ccd1024c3.ast
F3.ccd1053.ast | F3.ccd1053c5.ast


## <a name="Stage_CALIB"></a> CALIB

### Basic Explanation

This uses the transformation equations to convert the instrumental
magnitudes to calibrated magnitudes.

### Lists

Get the list of all ast files from ASTROM.outlist. All files that are
successfully calibrated are put into the calib.outlist (with .phot
endings). There will be a .phot file for each field chip/amp.

**INLIST** (ast) - Moved from ASTROM.outlist

Single-Chip | Split Multi-chip
------------ | -------------
F1.ccd1012.ast | F1.ccd1012_2.ast
F2.ccd1024.ast | F2.ccd1024c3.ast
F3.ccd1053.ast | F3.ccd1053c5.ast

**OUTLIST** (phot) - Every AST file that was successfully calibrated

Single-Chip | Split Multi-chip
------------ | -------------
F1.ccd1012.phot | F1.ccd1012_2.phot
F2.ccd1024.phot | F2.ccd1024c3.phot
F3.ccd1053.phot | F3.ccd1053c5.phot


## <a name="Stage_COMBINE"></a> COMBINE

### Basic Explanation

This combines all the photometry from the various chips/amps for multi- chip data.

### Lists

This combines all chips/amps of a multi-chip frame. It takes all the
files in CALIB.outlist and puts them into COMBINE.inlist. All files
that are successfully combined are put into the COMBINE.outlist. There
will be a separate output file for each field.

**INLIST** (phot) - Moved from CALIB.outlist

Single-Chip | Split Multi-chip
------------ | -------------
F1.ccd1012.phot | F1.ccd1012_2.phot
F2.ccd1024.phot | F2.ccd1024c3.phot
F3.ccd1053.phot | F3.ccd1053c5.phot

**OUTLIST** (cmb) - All of the chips/amps get combined. For single-chip
data the PHOT files are just copied over.

Single-Chip | Split Multi-chip
------------ | -------------
F1.ccd1012.cmb | F1.ccd1012.cmb
F2.ccd1024.cmb | F2.ccd1024.cmb
F3.ccd1053.cmb | F3.ccd1053.cmb


## <a name="Stage_DEREDDEN"></a> DEREDDEN

### Basic Explanation

This program dereddens the magnitudes (and colors specified in the
setup file) using the Schlegel maps.

### Lists

This takes all the files in COMBINE.outlist and dereddens the
magnitudes and colors. All files that are successful get put into
DEREDDEN.outlist There will be a separate output file for each field
with a .dered ending.

**INLIST** (ast) - Moved from combine.outlist

Single-Chip | Split Multi-chip
------------ | -------------
F1.ccd1012.cmb | F1.ccd1012.cmb
F2.ccd1024.cmb | F2.ccd1024.cmb
F3.ccd1053.cmb | F3.ccd1053.cmb

**OUTLIST** (dered) - All of the files that were successfully dereddened.

Single-Chip | Split Multi-chip
------------ | -------------
F1.ccd1012.dered | F1.ccd1012.dered
F2.ccd1024.dered | F2.ccd1024.dered
F3.ccd1053.dered | F3.ccd1053.dered


## <a name="Stage_SAVE"></a> SAVE

### Basic Explanation

This renames the final photometry files with the field name.

### Lists

Takes all of the files from deredden.outlist and renames them to have
the field names with .final extensions. It also saves IDL save files
of the final photometry structures.

**INLIST** (dered) - Moved from deredden.outlist

Single-Chip | Split Multi-chip
------------ | -------------
F1.ccd1012.dered | F1.ccd1012.dered
F2.ccd1024.dered | F2.ccd1024.dered
F3.ccd1053.dered | F3.ccd1053.dered

**OUTLIST** (final/dat/fits) - Rename the files with their field names. An IDL
save file of the photometry structure is saved to FIELD.dat and FITS binary
table saved to FIELD.fits.gz.

Single-Chip | Split Multi-chip
------------ | -------------
FIELD1.final/dat/fits | FIELD1.final/dat/fits
FIELD2.final/dat/fits | FIELD2.final/dat/fits
FIELD3.final/dat/fits | FIELD3.final/dat/fits


## <a name="Stage_HTML"></a> HTML

### Basic Explanation

This makes plots and summary HTML files for all the fields output by
SAVE. The plots and HTML files will be put in the "html/"
directory. The main index page is called "html/index.html"

### Lists

Takes all of the files from SAVE.outlist and creates plots and HTML summary pages.

**INLIST** (dat) - Moved from SAVE.outlist

Single-Chip | Split Multi-chip
------------ | -------------
F1.ccd1012.dat | F1.ccd1012.dat
F2.ccd1024.dat | F2.ccd1024.dat
F3.ccd1053.dat | F3.ccd1053.dat

**OUTLIST** NO outlist at the moment.
