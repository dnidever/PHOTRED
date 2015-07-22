pro stdred,dirs,redo=redo,stp=stp

;+
;
; STDRED
;
; This is the Standard Star portion of the PHOTRED photometry
; pipeline.  It takes final reduced images of standard star
; fields and output the transformation equations.
;
; STDRED should be run in ONE directory with all the standard
; star data from an ENTIRE run.
;
; A "stdred.setup" parameter should be in the directory with
; the necessary parameters.  here is an example:
;
;   ##### REQUIRED #####
;   scriptsdir  /net/home/dln5q/daophot/
;   irafdir     /net/home/dln5q/iraf/
;   telescope   Blanco
;   instrument  MOSAIC
;   observatory CTIO
;   ##### OPTIONAL #####
;   keepmef     0
;   redo        0
;   #wcsup       N
;   #wcsleft     E
;   #pixscale    0.5
;   wcsrefname  2MASS-PSC
;   #searchdist  60
;   wcsrmslim   0.5
;   matchdist   0.8
;   ##### STAGES #####
;   rename
;   split
;   wcs
;   aperphot
;   daogrow
;   astrom
;   matchcat
;   combinecat
;   fitdata
;
;
; There are 9 stages to STDRED.  Each one has it's own
; program (STDRED_stagename.PRO).  They can be run in order by
; putting them in the "stdred.setup" file, or they can be run
; separately from the command line.
;
; RENAME
; SPLIT
; WCS
; APERPHOT
; DAOGROW
; ASTROM
; MATCHCAT
; COMBINECAT
; FITDATA
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
;           Dn1-ccd1012.fits       Dn1-ccd1012.fits
;           Mn2-ccd1024.fits       Mn2-ccd1024.fits
;           Tn3-ccd1053.fits       Tn3-ccd1053c1.fits  IMACS
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
;           Dn1-ccd1012.fits       Dn1-ccd1012.fits
;           Mn2-ccd1024.fits       Mn2-ccd1024.fits
;           Tn3-ccd1053.fits       Tn3-ccd1053c1.fits  IMACS
;   OUTLIST: FITS  - All split files that are split okay, or
;                    single-chip files
;           Dn1-ccd1012.fits       Dn1-ccd1012_1.fits, Dn1-ccd1012_2.fits, ...
;           Mn2-ccd1024.fits       Mn2-ccd1024_1.fits, Mn2-ccd1024_2.fits, ...
;           Tn3-ccd1053.fits       Tn3-ccd1053c1.fits
;
; WCS
;   This program gets the correct WCS for images.
;
;   It takes all files from split.outlist and puts them into wcs.inlist.
;   All files that succeeded get put in wcs.outlist.
;   INLIST:  FITS  - Moved from rename.outlist
;           Single-Chip            Split Multi-chip
;           Dn1-ccd1012.fits       Dn1-ccd1012_1.fits
;           Mn2-ccd1024.fits       Mn2-ccd1024_1.fits
;           Tn3-ccd1053.fits       Tn3-ccd1053c1.fits
;   OUTLIST: FITS  - All object files that are given a proper wcs.
;           Dn1-ccd1012.fits       Dn1-ccd1012_1.fits
;           Mn2-ccd1024.fits       Mn2-ccd1024_1.fits
;           Tn3-ccd1053.fits       Tn3-ccd1053c1.fits
;
; APERPHOT
;   This program runs DAOPHOT PHOTOMETRY which gets aperture photometry for
;   all the stars in the image.
;
;   It takes all files from wcs.outlist and puts them into aperphot.inlist
;   All fits files that successfully run through daophot get put into
;   aperphot.outlist
;   INLIST:  FITS  - Moved from wcs.outlist
;           Single-Chip            Split Multi-chip
;           Dn1-ccd1012.fits       Dn1-ccd1012_1.fits
;           Mn2-ccd1024.fits       Mn2-ccd1024_1.fits
;           Tn3-ccd1053.fits       Tn3-ccd1053c1.fits
;   OUTLIST: AP  - All files in inlist that are succcessfully
;                  processed
;           Single-Chip            Split Multi-chip
;           Dn1-ccd1012.ap         Dn1-ccd1012_1.ap
;           Mn2-ccd1024.ap         Mn2-ccd1024_1.ap
;           Tn3-ccd1053.ap         Tn3-ccd1053c1.ap
;
; DAOGROW
;   This program corrects the aperture photometry for the aperture correction
;   (star by star) using DAOGROW.
;
;   It takes all files from aperphot.outlist and puts them into daogrow.inlist
;   All ap files that successfully run through daogrow get put into daogrow.outlist
;   INLIST:  AP  - Moved from aperphot.outlist 
;           Single-Chip            Split Multi-chip
;           Dn1-ccd1012.ap         Dn1-ccd1012_1.ap
;           Mn2-ccd1024.ap         Mn2-ccd1024_1.ap
;           Tn3-ccd1053.ap         Tn3-ccd1053c1.ap
;   OUTLIST: TOT  - Every ap files in inlist that is successfully processed and
;                      has a tot file from daogrow
;           Dn1-ccd1012.tot        Dn1-ccd1012_1.tot
;           Mn2-ccd1024.tot        Mn2-ccd1024_1.tot
;           Tn3-ccd1053.tot        Tn3-ccd1053c1.tot
;
; ASTROM
;   This program gets coordinates for all stars from the WCS in the image.
;
;   Takes all of the .tot files from daogrow.outlist and puts them into
;   astrom.inlist.  All files that are successfully given coordinates
;   are put into the astrom.outlist (with .ast endings).  There will
;   be a separate file for each field chip/amp.
;   INLIST:  TOT  - Moved from daogrow.outlist
;           Single-Chip            Split Multi-chip
;           Dn1-ccd1012.tot        Dn1-ccd1012_1.tot
;           Mn2-ccd1024.tot        Mn2-ccd1024_1.tot
;           Tn3-ccd1053.tot        Tn3-ccd1053c1.tot
;   OUTLIST: AST  - Every TOT file that was processed successfully
;           Dn1-ccd1012.ast        Dn1-ccd1012_1.ast
;           Mn2-ccd1024.ast        Mn2-ccd1024_1.ast
;           Tn3-ccd1053.ast        Tn3-ccd1053c1.ast
;
; MATCHCAT
;   This takes all of the .tot files and matches the stars to the standard
;   stars in the standard star catalogs.  It uses the coordinates to figure
;   out which field it is (i.e. SA98, SA110, etc.).
;
;   INLIST:  AST  - Moved from astrom.outlist
;           Single-Chip            Split Multi-chip
;           Dn1-ccd1012.ast        Dn1-ccd1012_1.ast
;           Mn2-ccd1024.ast        Mn2-ccd1024_1.ast
;           Tn3-ccd1053.ast        Tn3-ccd1053c1.ast
;   OUTLIST: CAT  - Each ast file that has some standard star matches.
;           Dn1-ccd1012.cat        Dn1-ccd1012_1.cat
;           Mn2-ccd1024.cat        Mn2-ccd1024_1.cat
;           Tn3-ccd1053.cat        Tn3-ccd1053c1.cat
;
; COMBINECAT
;   This takes all of the .cat files for each filter and combines them
;   into a combined .cat file.
;
;   INLIST:  CAT  - Moved from matchcat.outlist
;           Single-Chip            Split Multi-chip
;           Dn1-ccd1012.cat        Dn1-ccd1012_1.cat
;           Mn2-ccd1024.cat        Mn2-ccd1024_1.cat
;           Tn3-ccd1053.cat        Tn3-ccd1053c1.cat
;   OUTLIST: CAT  - One combined cat file for each filter observed
;           D.cat                  D.cat
;           M.cat                  M.cat
;           T.cat                  T.cat
;
; FITDATA
;   This takes the matched standard star catalogs for each filter and
;   fits the photometric transformation equations.
;
;   This takes all the files in combinecat.outlist and puts them into 
;   fitdata.inlist.
;   INLIST:  CAT  - Moved from combinecat.outlist
;           D.cat
;           M.cat
;           T.cat
;   OUTLIST: TRANS  - Transformation equations for each filter and
;                     each night.
;           D.trans
;           M.trans
;           T.trans
;           n1.trans
;           n2.trans
;           n3.trans
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
;
; Each stage can be run separately from the command line if desired.
;
;
; By D. Nidever  May 2008
;-



COMMON photred,setup


; Start the logfile
;----------------
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
logfile = 'stdred.'+smonth+sday+syear+shour+sminute+ssecond+'.log'
JOURNAL,logfile


; Check that all of the required programs are available
;------------------------------------------------------
;  Each sub-program will do its own test.
progs = ['stdred_rename','photred_split','photred_wcs','stdred_aperphot','stdred_daogrow',$
         'stdred_astrom','stdred_matchcat','stdred_combinecat','stdred_fitdata','readline',$
         'strsplitter','readpar','check_iraf','photred_getfilter','photred_getexptime','photred_getuttime',$
         'photred_getairmass','photred_getdate','photred_getgain','photred_getrdnoise',$
         'readlist','printlog','photred_loadsetup','photred_mkopt','mkopt','stdred_transphot',$
         'apcorrect','roi_cut','srcmatch','undefine','array_indices2','first_el','importascii','mad',$
         'maketemp','maxloc','minloc','mktemp','mpfit','photred_getinput','range','scale','scale_vector',$
         'slope','stringize','add_tag','getpixscale','head_xyad','head_adxy','loadaper',$
         'photred_updatelists','printstr','push','touchzero','wcsfit','wcsfit_imacs','matchstars',$
         'writecol','writeline','airmass','hdr2wcstnx','rndint','sexig2ten','signs','stress','strep',$
         'strmult','strtrim0','combine_structs','imfwhm','loadcoo','loadinput','mk_imacsmask','printline',$
         'randomize','rotsph','rotsphcen','wcstnx2hdr','wcstnx_rd2xy','wcstnx_xy2rd','parsetnx','badpar',$
         'sign','wcstnxcor','rd2xieta','xieta2rd']
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
PHOTRED_LOADSETUP,setup,count=count,/std
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
STDRED_RENAME,redo=redo


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


;---------
; APERPHOT
;---------
doaperphot = READPAR(setup,'APERPHOT')
if doaperphot ne '0' then $
STDRED_APERPHOT,redo=redo


;--------
; DAOGROW
;--------
dodaogrow = READPAR(setup,'DAOGROW')
if dodaogrow ne '0' then $
STDRED_DAOGROW,redo=redo


;-----------
; ASTROMETRY
;-----------
doastrom = READPAR(setup,'ASTROM')
if doastrom ne '0' then $
STDRED_ASTROM,redo=redo


;---------
; MATCHCAT
;---------
domatchcat = READPAR(setup,'MATCHCAT')
if domatchcat ne '0' then $
STDRED_MATCHCAT,redo=redo


;-----------
; COMBINECAT
;-----------
docombinecat = READPAR(setup,'COMBINECAT')
if docombinecat ne '0' then $
STDRED_COMBINECAT,redo=redo


;--------
; FITDATA
;--------
dofitdata = READPAR(setup,'FITDATA')
if dofitdata ne '0' then $
STDRED_FITDATA,redo=redo


print,'STDRED FINISHED'

; Run STDRED_SUMMARY
;--------------------
STDRED_SUMMARY


; End logfile
;------------
JOURNAL

if keyword_set(stp) then stop

end
