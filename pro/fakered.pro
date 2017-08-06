pro fakered,redo=redo


;FAKERED:
;-make another pipeline called FAKERED that does artificial star tests for
; fields.
;-You can run it in a fakered/ directory in a night's directory (e.g. n1/fakered/)
; or a separate directory for all nights.
;-The first stage should be COPY and it copies all of the files from that it needs
; (FITS images, mch, als, psf, alf, apcor, etc).  It also renames the files appropriately.
; You need to make a fakered.copyfields file that has the names of the fields you want
; and the field name.  For example,
; ../F5   84S251
; ../n2/F10   130L169a
; It also makes it's own "fields" file with renamed fields, e.g.
; F1  84S251
; F2  130L169a
; It might be good to keep the "original" filename somewhere so we can there where
; it came from.
;-There are many things that we don't need to redo: split, wcs, finding psf, 
; daomatch/daomaster transformations, etc.  We do however need to find the sources and
; extract the photometry etc. with ALLSTAR and do much of the ALLFRAME steps (including
; shifting/combining frames, Sextractor, ALLFRAME, etc.).
;-We can probably use modified versions of the PHOTRED programs for FAKERED.
;-For CALIB we always want to have the instrumental photometry, so that should be hardwired.
;
;STAGES:
;COPY
;ADDSTAR
;DAOPHOT  basically just find and allstar
;MATCH
;ALLFRAME
;CALIB
;COMBINE
;SAVE
;COMPLETE  calculates completeness
;
;
;There are two difficult parts about FAKERED: (1) what mag/color distribution to use for the
;artificial stars, (2) how to reliably tell if we recovered/found a star or not.
;It's important to cover the magnitude and color space as much as possible and go to brighter/fainter
;magnitudes in EACH filter.
;
;ADDSTAR (or fakestar)  this adds the artificial stars to multiple copies of each file.  Because each
;   time we do this we can only do ~100 stars and so we need to do this ~10x or so.  This could be parameter
;   you can specify in "fakered.setup".
;   -this should figure out the color/magnitude distribution to use for the entire field and make the
;    appropriate fake images for all amps/filters  (even if there are multiple exposures per filter).
;   -maybe the fake images should start with "T" or "A" or "AF" instead of "F" so we know it's a "fake" field.
;    also can make a separate "fakefields" (or something) file that tells you which "real" field
;    it "belongs" to:
;    T1 F1
;    T2 F1
;    T3 F1
;    T4 F2
;    T5 F2
;   -the aritificial stars need to be put into a separate file that can then be read by the COMPLETE stage
;    to see if the stars were recovered.
;   -this stage also need to make the psf, mch, apcor, alf combined psf and other files for each "fake" field.
;DAOPHOT
;  same as photred_daophot.pro but it runs a different script that just runs FIND and then uses the psf
;  (that we already have) and runs allstar.
;  I'M NOT SURE I NEED THIS.
;MATCH
;  just use the "mch" file and combine all of the photometry files into a "raw" file.
;ALLFRAME
;  run just like photred_allframe.pro  It's important that this is run with the same settings as the original
;  PHOTRED was run.  Don't need to do IMALIGN.  Just combine the frames.  Use the previously found combined psf.
;CALIB
;  run like in photred_calib.pro but always keep the "instrumental" magnitudes since the "input" magnitudes
;  are "instrumental".
;COMBINE
;  run like normal, but we won't have RA/DEC astronometry.  should we save it as an IDL save file?
;  photred_combine renames the IDs to "REALFIELD_AMP.NUMBER".  do we want that?
;COMPLETE
;  this figures out the completeness for each field.  For all the input files to this stage it uses the
;  "fakefields" file and figures out which "fake" fields belong to the "real" fields in the "fields" file.
;  -it has to use the original file created by ADDSTAR/FAKESTAR to see if a star was recovered or not.
;  -make a final completeness curve.  one separate one for each filter and then maybe one for color/filter
;   pairs.  this could be specified in "fakered.setup".  It could also make some plots.
;
;Really need to do things EXACTLY the way they were run in PHOTRED or the completeness won't be right.
;It's especially important if we run allframe or just allstar and HOW we run allframe (which program for
;detection, number of interations, etc.).
;
;Look at "startest.pro" that I wrote to do this.  Other ones to look at:
; checkfakestars.pro
; checkfakestars2.pro
; copyfakestars.pro   copies/makes images for fake files
; fakestar.pro        recovers the fake stars
; makefake.pro
; makefake2.pro       add fake stars to iamges
; makefakestars.pro
;Not sure which ones are the "final" or "last" versions.


;WCS
;-don't need to run
;DAOPHOT
;-keep opt, als.opt
;-keep psf
;-redo find, allstar
;-keep a.ap, a.als
;MATCH
;-keep mch
;-remake raw
;ALLFRAME
;-mch
;-psf
;-comb.psf
;-scale/weights/zero/shift
;-redo:
;  -shift,mask,combine
;  -keep the same zero, weights, scale, etc.
;  -keep the combined psf
;  -redo detection/find with daophot/sextractor
;  -run allframe
;  -recreate .mag file
;APCOR
;-keep apcor.lst
;ASTROM
;-run as normal
;CALIB
;-do we need this???
;COMBINE
;-run as normal
;DEREDDEN
;-don't need this
;SAVE
;-don't need this



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
logfile = 'fakered.'+smonth+sday+syear+shour+sminute+ssecond+'.log'
JOURNAL,logfile


; Check that all of the required programs are available
;------------------------------------------------------
;  Each sub-program will do its own test.
progs = ['photred_daophot','photred_match','photred_allframe','photred_astrom','photred_calib','photred_combine','readline',$
         'strsplitter','readpar','check_python','check_iraf','photred_getfilter','photred_getexptime','photred_getuttime',$
         'photred_getairmass','photred_getdate','photred_getgain','photred_getrdnoise',$
         'readlist','printlog','photred_loadsetup','photred_mkopt','mkopt','stdred_transphot',$
         'apcorrect','roi_cut','srcmatch','undefine','array_indices2','first_el','importascii','mad',$
         'maketemp','maxloc','minloc','mktemp','mpfit','photred_getinput','range','scale','scale_vector',$
         'slope','stringize','add_tag','getpixscale','head_xyad','head_adxy','loadaper',$
         'photred_updatelists','printstr','push','touchzero','wcsfit','wcsfit_imacs','matchstars',$
         'writecol','writeline','airmass','hdr2wcstnx','rndint','sexig2ten','signs','stress','strep',$
         'strmult','strtrim0','combine_structs','imfwhm','loadcoo','loadinput','mk_imacsmask','printline',$
         'randomize','rotsph','rotsphcen','wcstnx2hdr','wcstnx_rd2xy','wcstnx_xy2rd','parsetnx','badpar',$
         'sign','wcstnxcor','rd2xieta','xieta2rd', 'runshellcmd']
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
PHOTRED_LOADSETUP,setup,count=count,/fake
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


; do we need an "INIT" stage, that does the copying???

;-------
; COPY
;-------
; Copy the original files to the mock directories
docopy = READPAR(setup,'COPY')
if docopy ne '0' then $
FAKERED_COPY,redo=redo

;---------
; ADDSTAR
;---------
; Add the artificial sources
doaddstar = READPAR(setup,'ADDSTAR')
if doaddstar ne '0' then $
FAKERED_ADDSTAR,redo=redo

;---------
; DAOPHOT
;---------
; Basically just find and allstar
dodaophot = READPAR(setup,'DAOPHOT')
if dodaophot ne '0' then $
PHOTRED_DAOPHOT,/fake,redo=redo

;-------
; MATCH
;-------
; Recreate raw file
domatch = READPAR(setup,'MATCH')
if domatch ne '0' then $
PHOTRED_MATCH,/fake,redo=redo

;----------
; ALLFRAME
;----------
; Remake stack, detection and allframe
; Check if input list for ALLFRAME is coming from ADDSTAR 
; (when DAOPHOT and DAOMATH are NOT run)
; If so, copy ADDSTAR.outlist to ALLFRAME.inlist
if doaddstar ne '0' and dodaophot eq '0' and domatch eq '0' then $
  file_copy,'logs/ADDSTAR.outlist','logs/ALLFRAME.inlist',/OVERWRITE
doallframe = READPAR(setup,'ALLFRAME')
if doallframe ne '0' then $
PHOTRED_ALLFRAME,/fake,redo=redo

;--------
; ASTROM
;--------
; Add coordinates
doastrom = READPAR(setup,'ASTROM')
if doastrom ne '0' then $
PHOTRED_ASTROM,redo=redo

;-------------
; CALIBRATION
;-------------
; Calibrate
docalib = READPAR(setup,'CALIB')
if docalib ne '0' then $
PHOTRED_CALIB,redo=redo

;---------
; COMBINE
;---------
; Combine chips
docombine = READPAR(setup,'COMBINE')
if docombine ne '0' then $
PHOTRED_COMBINE,redo=redo

; Save, DO WE NEED THIS??
;PHOTRED_SAVE

;----------
; COMPLETE
;----------
; Calculate completeness
docomplete = READPAR(setup,'COMPLETE')
if docomplete ne '0' then $
FAKERED_COMPLETE,redo=redo


print,'FAKERED FINISHED'


; Run FAKERED_SUMMARY
;--------------------
FAKERED_SUMMARY


; End logfile
;------------
JOURNAL

if keyword_set(stp) then stop

end

