pro photred,dirs,redo=redo,stp=stp

;+
;
; PHOTRED
;
; This is a photometry pipeline that final reduced images
; and outputs calibrated PSF photometry with astrometry
; using the DAOPHOT suite of programs.
;
; PHOTRED should be run on ONE whole night.  All files
; should be in one directory.
;
; A "photred.setup" parameter should be in the directory with
; the necessary parameters.  here is an example:
;
;   ##### REQUIRED #####
;   scriptsdir  /net/home/dln5q/daophot/
;   irafdir     /net/home/dln5q/iraf/
;   telescope   Blanco
;   instrument  MOSAIC
;   nmulti      1
;   filtref     M
;   trans       blanco.trans
;   ##### OPTIONAL #####
;   keepmef     0
;   redo        0
;   #searchdist 60
;   #wcsrmslim  1.0
;   finditer    2
;   keepinstr   1
;   avgmag      1
;   avgonlymag  0
;   todered     M,T,D,M-T,M-D
;   #toextadd    M,T,D,M-T,M-D
;   ##### STAGES #####
;   rename
;   split
;   wcs
;   daophot
;   match
;   allframe
;   apcor
;   calib
;   astrom
;   combine
;   deredden
;   save
;
;
; There are 12 stages to PHOTRED.  Each one has it's own
; program (PHOTRED_stagename.PRO).  They can be run in order by
; putting them in the "photred.setup" file, or they can be run
; separately from the command line.
;
; RENAME
; SPLIT
; WCS
; DAOPHOT
; MATCH
; ALLFRAME
; APCOR
; CALIB
; ASTROM
; COMBINE
; DEREDDEN
; SAVE
;
;
; RENAME
;   This program renames object fits files so that it includes their
;   field information. For example, ccd1001.fits gets renamed to
;   F1.ccd1001.fits.  Any calibration frames (zero, dflat, sflat,
;   etc.) get moved to the "calib/" directory without getting renamed.
;
;   It will create the inlist from all fits files in the directory
;   Outlist will be of all files that aren't zero, dflat, sflat, etc.
;   and were successfully put in a "field".  It also creates the "fields" file.
;   Files are renamed: F1.obj1023.fits, F2.obj1045.fits, etc.
;   INLIST:  FITS  - Creates it itself from fits files in directory
;           Single-Chip            Multi-chip
;           zero1001.fits          zero1001.fits
;           ccd1012.fits           ccd1012.fits
;           ccd1024.fits           ccd1024.fits
;           ccd1053.fits           ccd1053c1.fits     IMACS
;   OUTLIST: FITS  - All object files that are properly renamed
;           F1.ccd1012.fits        F1.ccd1012.fits
;           F2.ccd1024.fits        F2.ccd1024.fits
;           F3.ccd1053.fits        F3.ccd1053c1.fits  IMACS
;   SUCCESS: FITS  - Same as outlist??
;
; SPLIT
;   This splits multiple-extension files (MEF) into separate files
;   for each amp/chip. Non-MEF files are not affected.
;
;   list of MEF or single chip files
;   All the files in rename.outlist are put in split.inlist.  All non-MEF
;   (single chip) files automatically go to split.outlist.  All
;   successfully split files (not original MEF files) are put in
;   split.outlist
;   INLIST:  FITS  - Moved from wcs.outlist
;           Single-Chip            Multi-chip
;           F1.ccd1012.fits        F1.ccd1012.fits 
;           F2.ccd1024.fits        F2.ccd1024.fits
;           F3.ccd1053.fits        F3.ccd1053c1.fits  IMACS
;   OUTLIST: FITS  - All split files that are split okay, or single-chip files
;           F1.ccd1012.fits        F1.ccd1012_1.fits, F1.ccd1012_2.fits, ...
;           F2.ccd1024.fits        F2.ccd1024_1.fits, F2.ccd1024_2.fits, ...
;           F3.ccd1053.fits        F3.ccd1053c1.fits
;   SUCCESS: FITS  - All inlist files that were successfully split or single-chip files
;           F1.ccd1012.fits        F1.ccd1012.fits
;           F2.ccd1024.fits        F2.ccd1024.fits
;           F3.ccd1053.fits        F3.ccd1053c1.fits  IMACS
;
; WCS
;   This program gets the correct WCS for images.
;
;   It takes all files from split.outlist and puts them into wcs.inlist.
;   All files that succeeded get put in wcs.outlist.
;   INLIST:  FITS  - Moved from rename.outlist
;           Single-Chip            Split Multi-chip
;           F1.ccd1012.fits        F1.ccd1012_1.fits, F1.ccd1012_2.fits, ...
;           F2.ccd1024.fits        F2.ccd1024_1.fits, F2.ccd1024_2.fits, ..
;           F3.ccd1053.fits        F3.ccd1053c1.fits     IMACS
;   OUTLIST: FITS  - All object files that are given a proper wcs.
;           F1.ccd1012.fits        F1.ccd1012_1.fits, F1.ccd1012_2.fits, ...
;           F2.ccd1024.fits        F2.ccd1024_1.fits, F2.ccd1024_2.fits, ..
;           F3.ccd1053.fits        F3.ccd1053c1.fits     IMACS
;   SUCCESS: FITS  - Same as outlist
;
;
; DAOPHOT
;   This program runs DAOPHOT on all the images (using Tony Sohn's script).
;
;   It takes all files from wcs.outlist and puts them into daophot.inlist
;   All fits files that successfully run through daophot get put into
;   daophot.outlist
;   INLIST:  FITS  - Moved from split.outlist
;           Single-Chip            Split Multi-chip
;           F1.ccd1012.fits        F1.ccd1012_1.fits, F1.ccd1012_2.fits, ...
;           F2.ccd1024.fits        F2.ccd1024_1.fits, F2.ccd1024_2.fits, ...
;           F3.ccd1053.fits        F3.ccd1053c1.fits  IMACS
;   OUTLIST: ALS  - All files in inlist that are succcessfully processed
;           F1.ccd1012.als         F1.ccd1012_1.als, F1.ccd1012_2.als, ...
;           F2.ccd1024.als         F2.ccd1024_1.als, F2.ccd1024_2.als, ...
;           F3.ccd1053.als         F3.ccd1053c1.als
;   SUCCESS: FITS  - All fits files that are successfully processed
;           F1.ccd1012.fits        F1.ccd1012_1.fits, F1.ccd1012_2.fits, ...
;           F2.ccd1024.fits        F2.ccd1024_1.fits, F2.ccd1024_2.fits, ...
;           F3.ccd1053.fits        F3.ccd1053c1.fits
;
; MATCH
;   This program runs DAOMATCH and DAOMASTER on the ALS files which
;   combines the photometry from the various filters.
;
;   It takes all the files from daophot.outlist and puts them into
;   match.inlist.  It uses the fields (from the filenames) to figure
;   out which files to together and should be matched.
;   All of the MCH files to in the match.outlist.
;   INLIST:  ALS  - Moved from daophot.outlist
;           Single-Chip            Split Multi-chip
;           F1.ccd1012.als         F1.ccd1012_1.als, F1.ccd1012_2.als, ...
;           F2.ccd1024.als         F2.ccd1024_1.als, F2.ccd1024_2.als, ...
;           F3.ccd1053.als         F3.ccd1053c1.als  IMACS
;   OUTLIST: MCH  - All files in inlist that were successfully matched (ALL??)
;           F1.ccd1012.mch         F1.ccd1012_1.mch, F1.ccd1012_2.mch, ...
;           F2.ccd1024.mch         F2.ccd1024_1.mch, F2.ccd1024_2.mch, ...
;           F3.ccd1053.mch         F3.ccd1053c1.mch
;   SUCCESS: ALS  - All als files that were successfully combined
;           F1.ccd1012.als         F1.ccd1012_1.als, F1.ccd1012_2.als, ...
;           F2.ccd1024.als         F2.ccd1024_1.als, F2.ccd1024_2.als, ...
;           F3.ccd1053.als         F3.ccd1053c1.als
;   There won't be as many MCH files in the outlist as ALS files in the inlist.
;   If there are 3 frames per field then there will be 3x as many ALS files as
;   MCH files.  The MCH files will have the names of the "reference" frame.
;
;   If some ALS files are not matched, but others are should this MCH file
;   still be put in the OUTLIST??
;
; ALLFRAME
;   This program runs ALLFRAME on all of the MCH files. ALLFRAME does
;   PSF fitting on images of all filters/bands similtaneously.
;
;   Takes the list of MCH files from the match.outlist and puts them
;   in allframe.inlist.  All fields that succeed, their MAG files get
;   put into allframe.outlist.
;   INLIST:  MCH  - Moved from match.outlist
;           Single-Chip            Split Multi-chip
;           F1.ccd1012.mch         F1.ccd1012_1.mch, F1.ccd1012_2.mch, ...  
;           F2.ccd1024.mch         F2.ccd1024_1.mch, F2.ccd1024_2.mch, ...  
;           F3.ccd1053.mch         F3.ccd1053c1.mch  IMACS
;   OUTLIST: MAG  - Every MCH file that was processed successfully by allframe
;                   and has a MAG file.
;           F1.ccd1012.mag         F1.ccd1012_1.mag, F1.ccd1012_2.mag, ...
;           F2.ccd1024.mag         F2.ccd1024_1.mag, F2.ccd1024_2.mag, ...
;           F3.ccd1053.mag         F3.ccd1053c1.mag
;   SUCCESS: MCH  - All mch files that were successfully processed
;           F1.ccd1012.mch         F1.ccd1012_1.mch, F1.ccd1012_2.mch, ...
;           F2.ccd1024.mch         F2.ccd1024_1.mch, F2.ccd1024_2.mch, ...
;           F3.ccd1053.mch         F3.ccd1053c1.mch
;
; APCOR
;   This programs find the aperture correction for all the files using DAOGROW.
;
;   Takes all of the fits files from daophot.success(!!) and puts them into
;   apcor.inlist.  All files that have an aperture correction in the final
;   apcor.lst get put into the apcor.outlist.
;   INLIST:  FITS  - COPIED from daophot.success 
;           Single-Chip            Split Multi-chip
;           F1.ccd1012.als         F1.ccd1012_1.als, F1.ccd1012_2.als, ...
;           F2.ccd1024.als         F2.ccd1024_1.als, F2.ccd1024_2.als, ...
;           F3.ccd1053.als         F3.ccd1053c1.als  IMACS
;   OUTLIST: FITS  - Every ALS/FITS file that was successfully given an aperture
;                   correction in the final apcor.lst file
;           F1.ccd1012.als         F1.ccd1012_1.als, F1.ccd1012_2.als, ...
;           F2.ccd1024.als         F2.ccd1024_1.als, F2.ccd1024_2.als, ...
;           F3.ccd1053.als         F3.ccd1053c1.als
;   SUCCESS: FITS  - Same as outlist.
;
; ASTROM
;   This program gets coordinates for all stars from the WCS in the
;   reference image.
;
;   Get the list of all mag files from allframe.outlist or mch files
;   from match.outlist.  All files that are successfully given coordinates
;   are put into ASTROM.outlist (with .ast endings).  There will be separate
;   file for each field chip/amp.
;   INLIST:  MAG/MCH  - Moved from ALLFRAME.outlist or if that does not exist, then MATCH.outlist
;           Single-Chip            Split Multi-chip
;           F1.ccd1012.mag         F1.ccd1012_1.mag, F1.ccd1012_2.mag, ...
;           F2.ccd1024.mag         F2.ccd1024_1.mag, F2.ccd1024_2.mag, ...
;           F3.ccd1053.mag         F3.ccd1053c1.mag  IMACS
;   OUTLIST: AST  - Every MAG/MCH file that was processed successfully
;           F1.ccd1012.ast         F1.ccd1012_1.ast, F1.ccd1012_2.ast, ...
;           F2.ccd1024.ast         F2.ccd1024_1.ast, F2.ccd1024_2.ast, ...
;           F3.ccd1053.ast         F3.ccd1053c1.ast
;   SUCCESS: PHOT  - All MAG/MCH files that were successfully processed
;           F1.ccd1012.mag         F1.ccd1012_1.mag, F1.ccd1012_2.mag, ...
;           F2.ccd1024.mag         F2.ccd1024_1.mag, F2.ccd1024_2.mag, ...
;           F3.ccd1053.mag         F3.ccd1053c1.mag
;
; CALIB
;   This uses the transformation equations to convert the instrumental
;   magnitudes to calibrated magnitudes.
;
;   Get the list of all ast files from ASTROM.outlist
;   All files that are successfully calibrated
;   are put into the CALIB.outlist (with .phot endings).
;   There will be a .phot file for each field chip/amp.
;   INLIST:  AST  - Moved from ASTROM.outlist
;           Single-Chip            Split Multi-chip
;           F1.ccd1012.ast         F1.ccd1012_1.ast, F1.ccd1012_2.ast, ...
;           F2.ccd1024.ast         F2.ccd1024_1.ast, F2.ccd1024_2.ast, ...
;           F3.ccd1053.ast         F3.ccd1053c1.ast  IMACS
;   OUTLIST: PHOT  - Every AST file that was successfully calibrated
;           F1.ccd1012.phot        F1.ccd1012_1.phot, F1.ccd1012_2.phot, ...
;           F2.ccd1024.phot        F2.ccd1024_1.phot, F2.ccd1024_2.phot, ...
;           F3.ccd1053.phot        F3.ccd1053c1.phot
;   SUCCESS: AST  - All ast files that were successfully processed
;           F1.ccd1012.ast         F1.ccd1012_1.ast, F1.ccd1012_2.ast, ...
;           F2.ccd1024.ast         F2.ccd1024_1.ast, F2.ccd1024_2.ast, ...
;           F3.ccd1053.ast         F3.ccd1053c1.ast
;
; COMBINE
;   This combines all the photometry from the various chips/amps for
;   multi-chip data.
;
;   This combines all chips/amps of a multi-chip frame.  It takes all
;   the files in CALIB.outlist and puts them into COMBINE.inlist.
;   All files that are successfully combined are put into the
;   COMBINE.outlist.  There will be a separate output file for each
;   field.
;   INLIST:  PHOT  - Moved from CALIB.outlist
;           Single-Chip            Split Multi-chip
;           F1.ccd1012.phot        F1.ccd1012_1.phot, F1.ccd1012_2.phot, ...
;           F2.ccd1024.phot        F2.ccd1024_1.phot, F2.ccd1024_2.phot, ...
;           F3.ccd1053.phot        F3.ccd1053c1.phot  IMACS
;   OUTLIST: AST  - All of the chips/amps get combined.  For single-chip
;                   data the PHOT files are just copied over.
;           F1.ccd1012.cmb         F1.ccd1012.cmb
;           F2.ccd1024.cmb         F2.ccd1024.cmb
;           F3.ccd1053.cmb         F3.ccd1053.cmb
;   SUCCESS: AST  - All PHOT files that were successfully processed
;           F1.ccd1012.phot        F1.ccd1012_1.phot, F1.ccd1012_2.phot, ...
;           F2.ccd1024.phot        F2.ccd1024_1.phot, F2.ccd1024_2.phot, ...
;           F3.ccd1053.phot        F3.ccd1053c1.phot
;
; DEREDDEN
;   This program dereddens the magnitudes (and colors specified in the
;   setup file) using the Schlegel maps. 
;
;   This takes all the files in combine.outlist and dereddens the magnitudes
;   and colors.  All files that are successful get put into deredden.outlist
;   There will be a separate output file for each field with a .dered
;   ending.
;   INLIST:  CMB  - Moved from COMBINE.outlist
;           Single-Chip            Split Multi-chip
;           F1.ccd1012.cmb         F1.ccd1012.cmb
;           F2.ccd1024.cmb         F2.ccd1024.cmb
;           F3.ccd1053.cmb         F3.ccd1053.cmb  IMACS
;   OUTLIST: DERED  - All of the files that were successfully dereddened.
;           F1.ccd1012.dered       F1.ccd1012.dered
;           F2.ccd1024.dered       F2.ccd1024.dered
;           F3.ccd1053.dered       F3.ccd1053.dered
;   SUCCESS: CMB  - All CMB files that were successfully processed
;           F1.ccd1012.cmb         F1.ccd1012.cmb
;           F2.ccd1024.cmb         F2.ccd1024.cmb
;           F3.ccd1053.cmb         F3.ccd1053.cmb
;
; SAVE
;   This renames the final photometry files with the field name.
;
;   Takes all of the files from deredden.outlist and renames them to have
;   the field names with .final extensions.  It also saves IDL save files of
;   the final photometry structures.
;   INLIST:  DERED  - Moved from deredden.outlist
;           Single-Chip            Split Multi-chip
;           F1.ccd1012.dered       F1.ccd1012.dered
;           F2.ccd1024.dered       F2.ccd1024.dered
;           F3.ccd1053.dered       F3.ccd1053.dered  IMACS
;   OUTLIST: FINAL/DAT - Rename the files with their field names.
;                        An IDL save file of the photometry structure
;                        is saved to FIELD.dat
;           FIELD1.final/dat       FIELD1.final/dat
;           FIELD2.final/dat       FIELD2.final/dat
;           FIELD3.final/dat       FIELD3.final/dat
;   SUCCESS: DERED  - All CMB files that were successfully processed
;           F1.ccd1012.dered       F1.ccd1012.dered
;           F2.ccd1024.dered       F2.ccd1024.dered
;           F3.ccd1053.dered       F3.ccd1053.dered
;
;
; Each stage has several log files associated with it:
;  INLIST   The list of files to process
;  OUTLIST  The successfully outputted files
;  SUCCESS  The files in INLIST that were successfully processed
;  FAILURE  The files in INLIST that were NOT successfully processed
;  LOG      A running log of what the stage has done
;
;
; If one file fails in match/combine then the whole set of images
; fails.  This will make it easy to "redo" them.
;
; Each stage can be run separately from the command line if desired.
;
; A file called "field" will be created by PHOTRED_RENAME with the
; unique field names and "shornames" that are prepended to file
; names.
;  F1   80S296
;  F2   42S285
;  F3   63S026
;  ...
;
;
; By D. Nidever  Jan-May 2008
;-



COMMON photred,setup


; Start the logfile
;------------------
; format is photred.DATETIME.log
jd = systime(/julian)
caldat,jd,month,day,year,hour,minute,second
smonth = strtrim(month,2)
if month lt 10 then smonth = '0'+smonth
sday = strtrim(day,2)
if day lt 10 then sday = '0'+sday
syear = strmid(strtrim(year,2),2,2)
shour = strtrim(hour,2)
if hour lt 10 then shour='0'+shour
sminute = strtrim(minute,2)
if minute lt 10 then sminute='0'+sminute
ssecond = strtrim(round(second),2)
if second lt 10 then ssecond='0'+ssecond
logfile = 'photred.'+smonth+sday+syear+shour+sminute+ssecond+'.log'
JOURNAL,logfile

; Print info
;-----------
host = GETENV('HOST')
print,''
print,'############################################'
print,'Starting PHOTRED   ',systime(0)
print,'Running on ',host
print,'############################################'
print,''

; Check that all of the required programs are available
;------------------------------------------------------
;  Each sub-program will do its own test.
progs = ['photred_rename','photred_split','photred_daophot','photred_wcs','range','add_tag',$
         'photred_apcor','photred_match','photred_calib','photred_combine','photred_astrom',$
         'photred_deredden','photred_save','readline','strsplitter','readpar','check_iraf',$
         'photred_allframe','pbs_daemon','pbs_makescript','pbs_checkstat','allframe','allfprep',$
         'photred_getfilter','photred_getexptime','photred_getuttime','mkdel','apcor','readlist',$
         'printlog','photred_loadsetup','photred_mkopt','photred_getgain','photred_getrdnoise',$
         'matchstars','roi_cut','srcmatch','undefine','array_indices2','dust_getval','first_el',$
         'importascii','mad','maketemp','maxloc','minloc','mktemp','mpfit','photred_getairmass',$
         'photred_getdate','photred_getinput','scale','scale_vector','slope','stringize','daomatch',$
         'getpixscale','head_xyad','head_adxy','loadmch','photcalib','photred_photcalib_prep',$
         'photred_updatelists','photmatch','phot_overlap','printstr','push','touchzero','writeline','wcsfit',$
         'wcsfit_imacs','airmass','bh_rdfort','djs_int2bin','djs_angpos','djs_ceil','hdr2wcstnx',$
         'odd','rndint','sexig2ten','signs','stress','strep','strmult','strtrim0','wcs_getval',$
         'combine_structs','imfwhm','iraf_imalign','iraf_imcombine','loadals','loadcoo','loadinput',$
         'loadraw','mkopt','mk_imacsmask','randomize','rotsph','rotsphcen','wcstnx2hdr','wcstnx_rd2xy',$
         'wcstnx_xy2rd','wcs_coord2pix','writecol','arr2str','badpar','parsetnx','wcstnxcor','sign',$
         'printline','rd2xieta','xieta2rd','writeals']
test = PROG_TEST(progs)
if min(test) eq 0 then begin
  bd = where(test eq 0,nbd)
  print,'SOME NECESSARY PROGRAMS MISSING'
  print,progs[bd]
  return
endif

; Make sure we have the right printlog.pro, not Markwardt's version
tempprogs = strsplit(!path,':',/extract)+'/printlog.pro'
test = file_test(tempprogs)
ind = where(test eq 1,nind)
bd = where(stregex(tempprogs[ind],'markwardt',/boolean) eq 1,nbd)
if nbd gt 0 then begin
  baddir = file_dirname(tempprogs[ind[bd]])
  print,"There is a version of Markwardt's PRINTLOG.PRO in "+baddir
  print,'Please rename this program (i.e. printlog.pro.orig)'
  return
endif


; LOAD THE SETUP FILE
;--------------------
; This is a 2xN array.  First colume are the keywords
; and the second column are the values.
; Use READPAR.PRO to read it
PHOTRED_LOADSETUP,setup,count=count
if (count lt 1) then return

; Are we redoing?
doredo = READPAR(setup,'REDO')
if keyword_set(redo) or (doredo ne '-1' and doredo ne '0') then redo=1


; Get the IRAF directory from the setup file
;-------------------------------------------
irafdir = READPAR(setup,'IRAFDIR')
if irafdir eq '0' or irafdir eq '-1' then irafdir = '~/iraf/'
irafdir = FILE_SEARCH(irafdir,/fully_qualify,count=nirafdir)
if nirafdir eq 0 then begin
  print,'NO IRAF DIRECTORY'
  return
endif


; Run CHECK_IRAF.PRO to make sure that you can run IRAF from IDL
;---------------------------------------------------------------
CHECK_IRAF,iraftest,irafdir=irafdir
if iraftest eq 0 then begin
  print,'IRAF TEST FAILED.  EXITING'
  return
endif


;#########################################
;#  STARTING THE PROCESSING
;#########################################


;-------
; RENAME
;-------
dorename = READPAR(setup,'RENAME')
if dorename ne '0' then $
PHOTRED_RENAME,redo=redo


;---------
; SPLIT
;---------
dosplit = READPAR(setup,'SPLIT')
if dosplit ne '0' then $
PHOTRED_SPLIT,redo=redo


;----
; WCS
;----
dowcs = READPAR(setup,'WCS')
if dowcs ne '0' then $
PHOTRED_WCS,redo=redo


;--------
; DAOPHOT
;--------
;   this includes acting as a pleione "daemon" and checking
;   up on the processes
dodaophot = READPAR(setup,'DAOPHOT')
if dodaophot ne '0' then $
PHOTRED_DAOPHOT,redo=redo


;------
; MATCH
;------
domatch = READPAR(setup,'MATCH')
if domatch ne '0' then $
PHOTRED_MATCH,redo=redo


;---------
; ALLFRAME
;---------
doallframe = READPAR(setup,'ALLFRAME')
if doallframe ne '0' then $
PHOTRED_ALLFRAME,redo=redo


;--------------------
; APERTURE CORRECTION
;--------------------
; Do an entire night together
doapcor = READPAR(setup,'APCOR')
if doapcor ne '0' then $
PHOTRED_APCOR,redo=redo


;-----------
; ASTROMETRY
;-----------
doastrom = READPAR(setup,'ASTROM')
if doastrom ne '0' then $
PHOTRED_ASTROM,redo=redo


;------------
; CALIBRATION
;------------
docalib = READPAR(setup,'CALIB')
if docalib ne '0' then $
PHOTRED_CALIB,redo=redo


;--------------
; COMBINE CHIPS
;--------------
docombine = READPAR(setup,'COMBINE')
if docombine ne '0' then $
PHOTRED_COMBINE,redo=redo


;---------
; DEREDDEN
;---------
doderedden = READPAR(setup,'DEREDDEN')
if doderedden ne '0' then $
PHOTRED_DEREDDEN,redo=redo


;-----------
; SAVE FILES
;-----------
dosave = READPAR(setup,'SAVE')
if dosave ne '0' then $
PHOTRED_SAVE,redo=redo


;-----------
; HTML
;-----------
dohtml = READPAR(setup,'HTML')
if dohtml ne '0' then $
PHOTRED_HTML,redo=redo


print,'PHOTRED FINISHED'


; Run PHOTRED_SUMMARY
;--------------------
PHOTRED_SUMMARY


; End logfile
;------------
JOURNAL

if keyword_set(stp) then stop

end
