# PHOTRED
Automated DAOPHOT-based PSF photometry pipeline

PHOTRED is an automated DAOPHOT-based PSF photometry pipeline.  It's very generic but does require reduced flat images.
PHOTRED does WCS fitting, aperture photometry, single-image PSF photometry (ALLSTAR), source matching across multiple images,
forced PSF photometry across multiple exposures using a master source list created from a deep stack of all exposures
(ALLFRAME), aperture correction, calibration using photometric transformation equations, and dereddening.
STDRED is also included and can be used to derive the transformation equations from exposures of standard star fields.

PHOTRED does require the following that need to be installed separately:
-IRAF
-IDL
-The IDL Astronomy Uers's Library
-The stand-alone fortran version of DAOPHOT/ALLFRAME
-SExtractor

PDF manuals for PHOTRED and STDRED are included in the doc/ directory.

IDL "sav" files of all the main PHOTRED programs are in the sav/ directory.  These can be used with the IDL virtual
machine for people who don't have an IDL license.

Installation Instructions

1. Download the PHOTRED IDL programs
Download the PHOTRED IDL programs tar file (last updated 06/02/08). Copy this to your IDL directory (most likely ~/idl/) and unpack it:

```
gunzip photred_idl.tar.gz tar -xvf photred_idl.tar
```

Let it overwrite any older programs by the same name. You need the new
versions!  You will also need the IDL Astro User's Library. The
programs are automatically available on UVa Astronomy. If you don't
have the IDL Astro User's Library then you should download the
programs from their download site (get "astron.tar.gz").

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

2. Download the PHOTRED scripts

Download the PHOTRED scripts tar file (last updated 06/02/08). Make a
directory where these scripts will reside (i.e. ~/photred/), and copy
the tar file to the directory and unpack it:

gunzip photred_scripts.tar.gz tar -xvf photred_scripts.tar

There are two fortran codes (lstfilter.f and makemag.f) that need to
be compiled. The tar file includes compiled versions that were
compiled on a Linux system. If you are planning to run PHOTRED on a
Sun machine or are having problems with the programs recompile them:

```
gfortran lstfilter.f -o lstfilter
```

3. Make sure IDL/IRAF are available

PHOTRED needs IDL and IRAF run. If you don't have IDL then you might
consider buying a license from ITT Visual Information
Solutions. Student licenses are available and full licenses have come
down in price. If these are too expensive, then you might consider
installing the free version of IDL, the GNU Data Language. I have not
tested PHOTRED on GDL, but I think GDL should have all of the
functionality needed to run PHOTRED. However, some PHOTRED programs
might need to be tweaked to work with GDL.  You can also download the
IDL Virtual Machine for free and run IDL "sav"

You can also download the IDL Virtual Machine for free and run IDL
"sav" files. I have made IDL sav files for the PHOTRED pipeline so
they can be run with the IDL Virtual Machine. Download the tar file
(last updated 05/12/08) and put them in your ~/idl/ directory. To run
a program type "idl -vm=progname.sav".  IRAF is freely downloadable
from iraf.net.

4. Make sure DAOPHOT/ALLFRAME is installed

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

5. Schlegel Maps

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
mycomputer % idl IDL>print,dust_getval(10,10)
 ￼￼0.472163
```

If you get an error here, then there is a problem. Check that all the
files and required programs are there.

6. Setup Your IRAF Login File

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

Running Instructions

1. Data

Start by putting all the final flat frames in their nights
directory. Process them with IRAF's CCDRED, MSCRED or Armin Rest's
SuperMacho pipeline to get final flat images. PHOTRED only runs on one
night's worth of data. If you have 5 nights of data then you will need
to run PHOTRED 5 times (they can be running simultaneously in separate
directories).

FIX BAD PIXELS!!!. Make sure to fix any bad pixels or columns in the
images because otherwise they might get detected as "sources" and mess
up WCS, DAOPHOT, MATCH, and ALLFRAME.

2. Transformation Equations

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

3. Setup File

PHOTRED needs a "photred.setup" file to run. This file specifies a few
important parameters. Here's an example of a "photred.setup" file. The
various parameters are described below.

```
##### REQUIRED #####
scriptsdir /net/home/dln5q/daophot/
irafdir /net/home/dln5q/iraf/
telescope Blanco
instrument MOSAIC
observatory CTIO
nmulti 1
filtref M
trans blanco.trans
##### OPTIONAL #####
keepmef
redo
#skipwcs
#wcsup
#wcsleft
#pixscale
#wcsrefname
#searchdist
#wcsrmslim
#hyperthread 1
psfcomsrc 1
#mchmaxshift 10.0
finditer 2
#alfdetprog sextractor
#ddo51radoffset 1
todered M,T,D,M-T,M-D
#toextadd M,T,D,M-T,M-D
keepinstr 1
avgmag 1
avgonlymag 0
#cmd2cdaxes
##### STAGES ##### rename
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
REQUIRED |  
scriptsdir |  The absolute path to the directory that contains the PHOTRED scripts (i.e. daophot.sh, etc.)
irafdir | The absolute path to your IRAF directory (that contains your login.cl file)
telescope | The name of the telescope (e.g. Blanco, Swope)
instrument | The name of the instrument (e.g. MOSAIC)
observatory | (OPTIONAL) The name of the observatory. This is needed if the header does not contain the AIRMASS and it needs to be calculated from the date, ra/dec and observatory location.
nmulti |  The number of processors PHOTRED should use for the DAOPHOT and ALLFRAME stages (only relevant if on the Pleione cluster).
filtref |  The shortname of the filter (specified in the "filters" file) to be used as the reference frame (e.g. the M frame). This can be a comma-delimited priority-ordered list of filters (e.g. g,i,r,z,u). If there are multiple observations in this filter then the longest exposure in this filter will be used.
trans | The name of the file that contains the photometric transformation equations.
keepmef | OPTIONAL. Multi-extension files (MEF) are split by PHOTRED. Do you want PHOTRED to keep the MEF files: YES=1, NO=0 (i.e. erase them).
redo | OPTIONAL. PHOTRED will NOT reprocess files that have already been processed unless "redo" is set. This can also be set as a keyword on the command line (i.e. IDL>photred,/redo).
skipwcs | OPTIONAL. Set this if your images already have correct WCS in their headers and you don't want the WCS to be refit in the WCS stage.
wcsup | OPTIONAL. What cardinal direction (i.e. N, S, E or W) is "up" in the image? This is only used for non-standard setups.
wcsleft | OPTIONAL. What cardinal direction (i.e. N, S, E or W) is "left" in the image? This is only used for non-standard setups.
pixscale | The plate scale in arcseconds/pixel. This is ONLY used for non-"standard" imagers (i.e. not MOSAIC, IMACS, LBC or Swope) where the pixel scale cannot be determined from the image headers.
wcsrefname | OPTIONAL. The name of the WCS reference catalog to use. The two options are 'USNO-B1' and '2MASS-PSC'. USNO-B1 is the default. The astrometric accuracy of the 2MASS catalog is better (~0.170 arcsec) than USNO-B1 (~0.270 arcsec), but it does not go as deep (R~18) as USNO-B1 (R~20). So if you have deep images then definitely use USNO-B1, but if you have moderately deep images then 2MASS-PSC is probably better (and faster).
searchdist | OPTIONAL. This sets the search distance (in arcmin) for WCS fitting (PHOTRED_WCS). Normally this is not needed. The default is 2*image size > 60 arcmin (i.e. whichever is greater). This is normally sufficient. If the WCS isn't fitting correctly then try setting "searchdist" to a larger value.
wcsrmslim | This is the maximum RMS (in arcseconds) allowed for an acceptable WCS fit. The default is 1.0 arcseconds. Normally the RMS values are ~0.2-0.3 arcseconds.
hyperthread | This allows multiple jobs to be running (daophot and allframe only) on a computer (such as halo or stream) that has multiple processors. It's similar to running it on a cluster.
psfcomsrc | If this is set to 1 then only sources detected in other frames of the same field are used as PSF stars. This makes sure that your PSF stars are not contaminated by cosmic rays or other junk. HIGHLY RECOMMENDED.
mchmaxshift | OPTIONAL. This sets a maximum constraint on the X/Y shifts found in the MATCH stage. WARNING, only use this if you ABSOLUTELY known that your frames are well aligned already.
finditer | The number of times to iteratively find sources in ALLFRAME (allfprep). The default is 2.
alfdetprog | OPTIONAL. The program to use for source detection in the ALLFRAME stage (allfprep). The options are "sextractor" and "daophot". The default is "sextractor". SExtractor is generally better at finding faint sources and returns a stellaricity probability value which is very useful. HOWEVER, SExtractor fails in VERY crowded regions. It's best to use DAOPHOT for very crowded images.
ddo51radoffset | OPTIONAL. There is a photometric offset in the DDO51 filter that depends on the radial distance from the center of the field. Currently this is only observed in the CTIO+MOSAIC data. Setting this parameter will remove this offset (done in CALIB). If you use this make sure to also use it in STDRED.
keepinstr | OPTIONAL. CALIB should keep the instrumental magnitudes in the final output file.
avgmag | OPTIONAL. CALIB should calculate average magnitudes in filters that were observed multiple times. The individual magnitudes are also kept.
avgonlymag | OPTIONAL. Same as "avgmag" but only keeps the average magnitudes.
todered | OPTIONAL. The magnitudes and colors the DEREDDEN stage should deredden. The short names of the filters and colors should be used (i.e. M, T, M-T). The dereddened magnitudes and colors will have the same names but with a "0" (zero) after it (i.e. M0 for M). The dash is removed for the color names.
toextadd | OPTIONAL. The extinction and reddening values to add. The short names of the filters and colors should be used (i.e. M, T, M-T). The extinctions will have the same name as the magnitude but with an "A" at the beginning (i.e. AM for M), and the reddening values will have the same name as the color names but with the dash removed and an "E" at the beginning (i.e. EMT for M-T).
cmd2cdaxes | The magnitudes and colors to use for the CMD and 2CD plots in PHOTRED_HTML.PRO. It should be in this format: MAG, MAG1-MAG2, MAG3-MAG4, i.e. M,M-T,M-D. For this example the CMD would be M vs. M-T and the 2CD would be M-D vs. M-T. The CMD and 2CD will be connected (CMD on top, 2CD on bottom) and will share the x-axis which be the first color. The magnitude will also be used for other plots (i.e. error vs. magnitude, chi vs. magnitude, etc).


Then add all the stage names that you want to process from the
following list: rename, split, wcs, daophot, match, allframe, apcor,
calib, astrom, combine, deredden, and save

The number of ALLFRAME iterations can be set in the "allframe.opt" file, the MA parameter. The default is 50.