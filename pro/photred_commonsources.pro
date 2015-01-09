;+
;
; PHOTRED_COMMONSOURCES
;
; This gets confirmed celestial sources for a frame by making sure
; they appear in other frames of the same field.  These will be used
; to pick PSF stars by daophot.sh.  It uses the field information
; appended to the filename to find other frames of the same field.
;
; INPUTS:
;  input    The FITS file(s) for which to find confirmed celestial
;             sources.  It should have the field information
;             appended, i.e. F6-ccd1001.fits. 
;             Three formats can be used (1) Name of file
;             with a list of input filenames.  Must start with an '@';
;             (2) A name with wildcard characters, such as '*';
;             (3) An array of filenames.
;  =minsources  The number of minimum sources required.  The default
;                 is 6 sources because that's what DAOPHOT needs to
;                 fit a PSF.  If there are less than this number then
;                 no .cmn.lst file will be created.
;  =maxframes   Only use a maximum number of common frames. Default is 5.
;  /gaussfit    Fit Gaussians to the common sources to weed out junk.
;  /redo        Redo if it already exists, otherwise exit.
;  /stp         Stop at the end of the program.
;
; OUTPUTS:
;  A file will be created with the list of confirmed celestial
;  sources, i.e. sources that were detected in other images as well.
;  The filename will be BASE.cmn.lst and be in a "standard" DAOPHOT
;  format with three header lines and columns of:
;  ID, X, Y, MAG, ERR, SKY, SKYSIG, SHARP, ROUND, ROUND2
;  =error   The error message if one occured
;
; USAGE:
;  IDL>photred_commonsources,file,minsources=minsources,gaussfit=gaussfit,error=error,stp=stp'
;
; By D. Nidever   June 2008
;-

pro photred_commonsources,input,minsources=minsources0,redo=redo,gaussfit=gaussfit,$
                          error=error,stp=stp,maxframes=maxframes

COMMON photred,setup

undefine,error

; Not enough inputs
ninput = n_elements(input)
if ninput eq 0 then begin
  print,'Syntax - photred_commonsources,input,minsources=minsources,gaussfit=gaussfit,'
  print,'                                maxframes=maxframes,error=error,stp=stp'
  error = 'Not enough inputs'
  return
endif

; Loading input
LOADINPUT,input,file,count=nfiles

; Not enough inputs
if nfiles eq 0 then begin
  print,'No files'
  return
endif

; More than one name input
if nfiles gt 1 then begin
  for i=0,nfiles-1 do begin
    PHOTRED_COMMONSOURCES,file[i],minsources=minsources0,redo=redo,gaussfit=gaussfit,$
                          error=error,stp=stp,maxframes=maxframes
    if keyword_set(verbose) then print,''
  endfor
  return
endif


; Is this a FITS file
len = strlen(file)
if (strmid(strtrim(file,2),len-5,5) ne '.fits') then begin
  print,file,' is NOT a FITS file'
  error = file+' is NOT a FITS file'
  return
endif

; Does the FITS file exist
if FILE_TEST(file) eq 0 then begin
  print,file,' NOT FOUND'
  error = file+' NOT FOUND'
  return
endif

; Defaults
if n_elements(minsources0) eq 0 then minsources=6 else minsources=minsources0
minsources = long(minsources) > 1
if n_elements(maxframes) eq 0 then maxframes=5

; LOAD THE SETUP FILE if not passed
;-----------------------------------
; This is a 2xN array.  First colume are the keywords
; and the second column are the values.
; Use READPAR.PRO to read it
if n_elements(setup) eq 0 then begin
  PHOTRED_LOADSETUP,setup,count=count
  if count lt 1 then return
endif

; Log files
;----------
;  write to DAOPHOT logfile
logfile = 'logs/DAOPHOT.log'
logfile = FILE_EXPAND_PATH(logfile)  ; want absolute filename
if file_test(logfile) eq 0 then SPAWN,'touch '+logfile,out

; Telescope, Instrument
telescope = READPAR(setup,'TELESCOPE')
telescope = strupcase(telescope)
instrument = READPAR(setup,'INSTRUMENT')
instrument = strupcase(instrument)
; PSFCOMGAUSS
dogaussfit = READPAR(setup,'PSFCOMGAUSS')
if keyword_set(gaussfit) or (dogaussfit ne '-1' and dogaussfit ne '0') then gaussfit=1

; Get the scripts directory from setup
scriptsdir = READPAR(setup,'SCRIPTSDIR')
if scriptsdir eq '' then begin
  printlog,logfile,'NO SCRIPTS DIRECTORY'
  return
endif

; LOAD THE "imagers" FILE
;----------------------------
;printlog,logfile,'Loading imager information'
imagerstest = FILE_TEST(scriptsdir+'/imagers')
if (imagerstest eq 0) then begin
  printlog,logfile,'NO >>imagers<< file in '+scriptsdir+'  PLEASE CREATE ONE!'
  return
endif
; The columns need to be: Telescope, Instrument, Namps, separator
imagers_fieldnames = ['telescope','instrument','observatory','namps','separator']
imagers_fieldtpes = [7,7,7,3,7]
imagers = IMPORTASCII(scriptsdir+'/imagers',fieldnames=imagers_fieldnames,$
                      fieldtypes=imagers_fieldtypes,comment='#',/silent)
imagers.telescope = strupcase(strtrim(imagers.telescope,2))
imagers.instrument = strupcase(strtrim(imagers.instrument,2))
imagers.observatory = strupcase(strtrim(imagers.observatory,2))
singleind = where(imagers.namps eq 1,nsingle)
if nsingle gt 0 then imagers[singleind].separator = ''
if (n_tags(imagers) eq 0) then begin
  printlog,logfile,'NO imagers in '+scriptsdir+'/imagers'
  return
endif

; What IMAGER are we using??
;---------------------------
ind_imager = where(imagers.telescope eq telescope and imagers.instrument eq instrument,nind_imager)
if nind_imager eq 0 then begin
  printlog,logfile,'TELESCOPE='+telescope+' INSTRUMENT='+instrument+' NOT FOUND in >>imagers<< file'
  return
endif
thisimager = imagers[ind_imager[0]]



; Get the field information
base = FILE_BASENAME(file,'.fits')
basarr = strsplit(base,'-',/extract)
dash = strpos(base,'-')
if (strmid(base,0,1) ne 'F' or dash[0] eq -1) then begin
  print,'NO field information in ',base
  error = 'NO field information in '+base
  return
endif
field = basarr[0]


; Does the file exist already?
test = FILE_TEST(base+'.cmn.lst')
if test eq 1 and not keyword_set(redo) then begin
  print,'FILE '+base+'.cmn.lst EXISTS ALREADY.'
  return
endif


print,'Making list of CONFIRMED CELESTIAL SOURCES for ',file


;#########################################
;#   Find other files of this field
;#########################################


; Check the logs/RENAME.outlist
nfieldfiles = 0
if FILE_TEST('logs/RENAME.outlist') eq 1 then begin
  READLIST,'logs/RENAME.outlist',fieldfiles,/exist,/unique,/fully,/silent,count=nfieldfiles
  if nfieldfiles gt 0 then begin
    gf = where(stregex(FILE_BASENAME(fieldfiles),'^'+field+'-',/boolean) eq 1,ngf)
    if ngf gt 0 then begin
      fieldfiles = fieldfiles[gf]
      nfieldfiles = ngf
    endif else begin
      undefine,fieldfiles
      nfieldfiles = 0
    endelse
  endif else begin
    undefine,fieldfiles
    nfieldfiles = 0
  endelse
endif

; Search for files in the directory.
; If less than 2 files found
if (nfieldfiles lt 2) then begin
  fieldfiles = FILE_SEARCH(field+'-*.fits',count=nfieldfiles)

  ; Remove a.fits, s.fits, _comb.fits and other "temporary" files.
  if nfieldfiles gt 0 then begin
    fbases = FILE_BASENAME(fieldfiles,'.fits')
    bad = where(stregex(fbases,'a$',/boolean) eq 1 or $         ; psf stars image
                stregex(fbases,'s$',/boolean) eq 1 or $         ; allstar subtracted file
                stregex(fbases,'_comb$',/boolean) eq 1 or $     ; stacked field image
                stregex(fbases,'_comb_sub$',/boolean) eq 1 or $ ; allstar subtracted stacked image
                stregex(fbases,'j$',/boolean) eq 1 or $         ; allframe temp file
                stregex(fbases,'k$',/boolean) eq 1 or $         ; allframe temp file
                stregex(fbases,'jnk$',/boolean) eq 1,nbad)      ; daophot? temp file
    if nbad gt 0 then begin
      if nbad eq nfieldfiles then begin
        undefine,fieldfiles
        nfieldfiles = 0
      endif else begin
        REMOVE,bad,fieldfiles
        nfieldfiles = n_elements(fieldfiles)
      endelse
    endif ; some ones to remove
  endif ; some fieldfiles

endif

; Make sure the current file is in the FIELDFILES list
allfiles = [fieldfiles, file]
fieldfiles = FILE_SEARCH(allfiles,/fully_qualify,count=nfieldfiles)
ui = uniq(fieldfiles,sort(fieldfiles))
fieldfiles = fieldfiles[ui]
nfieldfiles = n_elements(fieldfiles)


; If multiple chips, get the correct chip files
if thisimager.namps gt 1 then begin

  amp = first_el(strsplit(base,thisimager.separator,/extract),/last)
  fbases = FIlE_BASENAME(fieldfiles,'.fits')
  gdbase = where(stregex(fbases,thisimager.separator+amp+'$',/boolean) eq 1,ngdbase)

  if (ngdbase gt 0) then begin
    fieldfiles = fieldfiles[gdbase]
    nfieldfiles = ngdbase
  endif else begin
    print,'NO frames for this Field/amp combination  ',field,'/',amp
    error = 'NO frames for this Field/amp combination  '+field+'/'+amp
    return
  endelse

endif


;; If MOSAIC or LBC, get the correct CHIP files
;underscore = strpos(base,'_')
;if (instrument eq 'MOSAIC' or instrument eq 'LBC' and underscore[0] ne -1) then begin
;  amp = first_el(strsplit(base,'_',/extract),/last)
;  fbases = FIlE_BASENAME(fieldfiles,'.fits')
;  gdbase = where(stregex(fbases,'_'+amp+'$',/boolean) eq 1,ngdbase)
;
;  if (ngdbase gt 0) then begin
;    fieldfiles = fieldfiles[gdbase]
;    nfieldfiles = ngdbase
;  endif else begin
;    print,'NO frames for this Field/amp combination  ',field,'/',amp
;    error = 'NO frames for this Field/amp combination  '+field+'/'+amp
;    return
;  endelse
;
;endif
;
;; If IMACS, get the correct CHIP files
;imacs = stregex(base,'c[1-8]$',/boolean)
;if (instrument eq 'IMACS' and imacs eq 1) then begin
;
;  len = strlen(base)
;  amp = strmid(base,len-1,1)
;  fbases = FIlE_BASENAME(fieldfiles,'.fits')
;  gdbase = where(stregex(fbases,'c'+amp+'$',/boolean) eq 1,ngdbase)
;
;  if (ngdbase gt 0) then begin
;    fieldfiles = fieldfiles[gdbase]
;    nfieldfiles = ngdbase
;  endif else begin
;    print,'NO frames for this Field/chip combination  ',field,'/',amp
;    error = 'NO frames for this Field/chip combination  '+field+'/'+amp
;    return
;  endelse
;
;endif


print,'Found ',strtrim(nfieldfiles,2),' frames of FIELD=',field


;###############################################
;#   Get Sources with DAOPHOT FIND/PHOTOMETRY
;###############################################

; Only take a certain number of frames
if keyword_set(maxframes) then begin
  ; Want maxframes plust the current file
  if maxframes+1 lt nfieldfiles then begin
    ; Make sure the current file is in the FIELDFILES list
    fieldfiles = [fieldfiles[0:maxframes],file]
    fieldfiles = FILE_SEARCH(fieldfiles,/fully_qualify,count=nfieldfiles)
    ui = uniq(fieldfiles,sort(fieldfiles))
    fieldfiles = fieldfiles[ui]
    nfieldfiles = n_elements(fieldfiles)
  endif
endif


; Run DAOPHOT FIND/PHOT on ALL the field files
For i=0,nfieldfiles-1 do begin

  ifile = fieldfiles[i]
  ibase = FILE_BASENAME(ifile,'.fits')

  ; Do the CMN.COO and CMN.AP files already exist?
  coofile = ibase+'.cmn.coo'
  cootest = FILE_TEST(coofile)
  if cootest eq 1 then coolines=FILE_LINES(coofile) else coolines=0
  apfile = ibase+'.cmn.ap'
  aptest = FILE_TEST(apfile)
  if aptest eq 1 then aplines=FILE_LINES(apfile) else aplines=0

  ; Run DAOPHOT FIND and PHOTOMETRY
  if (coolines lt 4 or aplines lt 4) then begin

    ;; Use the .COO and .AP files
    ;coofile1 = ibase+'.coo'
    ;cootest1 = FILE_TEST(coofile1)
    ;if cootest1 eq 1 then coolines1=FILE_LINES(coofile1) else coolines1=0
    ;apfile1 = ibase+'.ap'
    ;aptest1 = FILE_TEST(apfile1)
    ;if aptest1 eq 1 then aplines1=FILE_LINES(apfile1) else aplines1=0
    ;if coolines1 gt 3 and aplines1 gt 3 then begin
    ;  print,'Using .COO and .AP files for ',ibase
    ;  FILE_COPY,[coofile1,apfile1],[coofile,apfile],/overwrite,/allow
    ;  goto,BOMB_DAOPHOT
    ;endif


    ; RUN DAOPHOT
    print,'Getting sources for ',ibase,' using DAOPHOT'

    ; Make a .opt file
    optfile = ibase+'.opt'
    opttest = FILE_TEST(optfile)
    if opttest eq 0 then begin
      PHOTRED_MKOPT,ifile,fwhm=fwhm,error=mkopterror
    endif else begin
      READLINE,optfile,optlines
      fwhmind = where(stregex(optlines,'FW =',/boolean) eq 1,nfwhmind)
      fwhmarr = strsplit(optlines[fwhmind[0]],'=',/extract)
      fwhm = float(fwhmarr[1])
    endelse

    ; Problems with .opt file
    if FILE_TEST(optfile) eq 0 then goto,BOMB_DAOPHOT
 
    ; Copy the .opt file daophot.opt 
    if FILE_TEST('daophot.opt') eq 0 then $
      FILE_COPY,ibase+'.opt','daophot.opt',/over,/allow
 
    ; Make temporary photo.opt file
    tphotofile = MKTEMP('photo')
    undefine,photlines
    push,photlines,'A1 = '+STRING(fwhm*3.0,format='(F7.4)')
    push,photlines,'IS = 45.0000'
    push,photlines,'OS = 50.0000'
    WRITELINE,tphotofile,photlines
 
    ; Make a temporary script to run FIND
    undefine,lines
    push,lines,'#!/bin/sh'
    ;push,lines,'daophot="/net/astro/bin/daophot"'
    push,lines,'export image=${1}'
    push,lines,'rm ${image}.cmn.log      >& /dev/null'
    push,lines,'rm ${image}.cmn.coo      >& /dev/null'
    push,lines,'rm ${image}.cmn.ap       >& /dev/null'
    push,lines,'daophot << END_DAOPHOT >> ${image}.cmn.log'
    push,lines,'OPTIONS'
    push,lines,'${image}.opt'
    push,lines,''
    push,lines,'ATTACH ${image}.fits'
    push,lines,'FIND'
    push,lines,'1,1'
    push,lines,'${image}.cmn.coo'
    push,lines,'y'
    push,lines,'PHOTOMETRY'
    push,lines,FILE_BASENAME(tphotofile)
    push,lines,''
    ; Erase any PSF file for this frame
    psffile = ibase+'.psf'
    if FILE_TEST(psffile) eq 1 then begin
      FILE_DELETE,psffile,/allow,/quiet
      ;if FILE_LINES(psffile) gt 0 then $
      ;  push,lines,'END-OF-FILE'  ; don't use profile-fitting
    endif
    push,lines,'${image}.cmn.coo'
    push,lines,'${image}.cmn.ap'
    push,lines,'EXIT'
    push,lines,'END_DAOPHOT'
    ;tempfile = maketemp('dao','.sh')
    tempfile = MKTEMP('dao')    ; absolute path
    WRITELINE,tempfile,lines
    FILE_CHMOD,tempfile,'755'o

    ; Run the program
    SPAWN,tempfile+' '+ibase,out,errout
    ;FILE_DELETE,tempfile    ; delete the temporary script
 
    ; Test the coo and ap file
    cootest = FILE_TEST(ibase+'.cmn.coo')
    if cootest eq 1 then coolines=FILE_LINES(ibase+'.cmn.coo') else coolines=0
    aptest = FILE_TEST(ibase+'.cmn.ap')
    if aptest eq 1 then aplines=FILE_LINES(ibase+'.cmn.ap') else aplines=0
 
    ; Remove the temporary files
    FILE_DELETE,tempfile    ; delete the temporary script
    FILE_DELETE,'daophot.opt',/allow
    FILE_DELETE,tphotofile,/allow
    junk = FILE_TEST(ibase+'jnk.fits')
    if junk eq 1 then FILE_DELETE,ibase+'jnk.fits'
    ;FILE_DELETE,ibase+['.cmn.log','.cmn.coo','.cmn.ap'],/allow
    ;FILE_DELETE,ibase+['.opt','.als.opt'],/allow

  endif

  BOMB_DAOPHOT:

Endfor



;###############################################
;#   Match the sources
;###############################################


; Load the CURRENT file
coofile = base+'.cmn.coo'
cootest = FILE_TEST(coofile)
if cootest eq 1 then coolines=FILE_LINES(coofile) else coolines=0
apfile = base+'.cmn.ap'
aptest = FILE_TEST(apfile)
if aptest eq 1 then aplines=FILE_LINES(apfile) else aplines=0

; DAOPHOT ran properly
if (coolines ge 4 and aplines ge 4) then begin

  ; Get the header
  fitsfile = base+'.fits'
  head = headfits(fitsfile)

  ; Load the coordinates file
  LOADCOO,coofile,coo,coohead1
  ncoo = n_elements(coo)

  ; Load the aperture photometry file
  LOADAPER,apfile,aper,aperhead1

  ; Get the coordinates
  EXTAST,head,astr
  ra=coo.x*0.-999999. & dec=ra    ; BAD until proven okay
    if n_elements(astr) gt 0 then $
  HEAD_XYAD,head,coo.x,coo.y,ra,dec,/deg

  ; Create the CAT structure
  dum = {id:0L,x:0.0d0,y:0.0d0,mag:0.0,err:0.0,sky:0.0,skysig:0.0,sharp:0.0,round:0.0,round2:0.0,$
         ra:0.0d0,dec:0.0d0}
  cat = replicate(dum,ncoo)
  cat.id = coo.id
  cat.x = coo.x
  cat.y = coo.y
  cat.sharp = coo.sharp
  cat.round = coo.round
  cat.round2 = coo.round2
  cat.mag = aper.mag[0]
  cat.err = aper.err[0]
  cat.sky = aper.sky
  cat.skysig = aper.skysig
  cat.ra = ra
  cat.dec = dec

  ; Only keep sources with "decent" photometry
  cat_orig = cat
  gdcat = where(cat.mag lt 50. and cat.err lt 5.0,ngdcat)
  ;gdcat = where(cat.mag lt 50. and cat.err lt 5.0 and cat.sharp lt 2.0,ngdcat)
  if (ngdcat eq 0) then begin
    print,base,' has NO sources with good photometry'
    error = base+' has NO sources with good photometry'
    return
  endif
  cat = cat[gdcat]
  ncat = n_elements(cat)

; DAOPHOT did NOT run properly
endif else begin
  print,'DAOPHOT did NOT run properly on ',base
  error = 'DAOPHOT did NOT run properly on '+base
  return
endelse


; Match the sources
ndetected = lonarr(ncat)  ; The number of detections in other frames
For i=0,nfieldfiles-1 do begin

  ifile = fieldfiles[i]
  ibase = FILE_BASENAME(ifile,'.fits')

  if ibase eq base then goto,BOMB  ; Skip the current file

  ; Do the CMN.COO and CMN.AP files exist?
  coofile = ibase+'.cmn.coo'
  cootest = FILE_TEST(coofile)
  if cootest eq 1 then coolines=FILE_LINES(coofile) else coolines=0
  apfile = ibase+'.cmn.ap'
  aptest = FILE_TEST(apfile)
  if aptest eq 1 then aplines=FILE_LINES(apfile) else aplines=0

  ; DAOPHOT ran properly
  if (coolines ge 4 and aplines ge 4) then begin
    undefine,coo,aper

    ; Get the header
    fitsfile = ibase+'.fits'
    head1 = headfits(fitsfile)

    ; Load the coordinates file
    LOADCOO,coofile,coo,coohead
    ncoo = n_elements(coo)

    ; Load the aperture photometry file
    LOADAPER,apfile,aper,aperhead

    ; Get the coordinates
    EXTAST,head1,astr1
    ra=coo.x*0.-999999. & dec=ra    ; BAD until proven okay
    if n_elements(astr1) gt 0 then $
      HEAD_XYAD,head1,coo.x,coo.y,ra,dec,/deg

    ; Create the CAT structure
    dum = {id:0L,x:0.0d0,y:0.0d0,mag:0.0,err:0.0,sky:0.0,skysig:0.0,sharp:0.0,round:0.0,round2:0.0,$
           ra:0.0d0,dec:0.0d0}
    refcat = replicate(dum,ncoo)
    refcat.id = coo.id
    refcat.x = coo.x
    refcat.y = coo.y
    refcat.sharp = coo.sharp
    refcat.round = coo.round
    refcat.round2 = coo.round2
    refcat.mag = aper.mag[0]
    refcat.err = aper.err[0]
    refcat.sky = aper.sky
    refcat.skysig = aper.skysig
    refcat.ra = ra
    refcat.dec = dec

    ; Only keep sources with "decent" photometry
    refcat_orig = refcat
    gdref = where(refcat.mag lt 50. and refcat.err lt 5.0,ngdref)
    ;gdref = where(refcat.mag lt 50. and refcat.err lt 5.0 and refcat.sharp lt 2.0,ngdref)
    if (ngdref eq 0) then begin
      print,ibase,' has NO sources with good photometry'
      goto,BOMB
    endif
    refcat = refcat[gdref]


    ; MATCH sources with X/Y
    MATCHSTARS,cat.x,cat.y,refcat.x,refcat.y,ind1_xy,ind2_xy,trans,count=nmatch_xy,rms=rms_xy,/silent

    ; MATCH sources with RA/DEC
    nmatch_ad = -1
    rms_ad = -1.0
    if max(cat.ra) gt -1000. and max(refcat.ra) gt -1000. then begin
      SRCMATCH,cat.ra,cat.dec,refcat.ra,refcat.dec,1.0,ind1_ad,ind2_ad,count=nmatch_ad,/sph
      if (nmatch_ad gt 10) then begin
        raresid = (cat[ind1_ad].ra-refcat[ind2_ad].ra)*cos(cat[ind1_ad].dec/!radeg)*3600.
        decresid = (cat[ind1_ad].dec-refcat[ind2_ad].dec)*3600.
        resid = sqrt( raresid^2.0 + decresid^2.0 )
        rms_ad =  sqrt( mean( resid^2. ) )
      endif
    end

    ; No matches 
    if nmatch_xy eq 0 and nmatch_ad eq 0 then begin
      print,'NO matches'
      error = 'NO matches'
      return
    endif

    ; Use X/Y or RA/DEC Matches?
    mcase = 0   ; X/Y by default
    if nmatch_ad lt 0 then mcase=0    ; NO RA/DEC, use X/Y
    if nmatch_ad gt 10 then mcase=1   ; Use RA/DEC
    if nmatch_xy eq 0 and nmatch_ad gt 0 then mcase=1
    CASE mcase of
      0: begin   ; use X/Y
          ind1 = ind1_xy
          ind2 = ind2_xy
          ndetected[ind1_xy]++      ; increment the star Ndetected array
          print,ibase,'  Nmatch = ',strtrim(nmatch_xy,2),'  RMS=',strtrim(rms_xy,2),' pixels'
      end
      1:  begin  ; use RA/DEC
          ind1 = ind1_ad
          ind2 = ind2_ad
          ndetected[ind1_ad]++      ; increment the star Ndetected array
          print,ibase,'  Nmatch = ',strtrim(nmatch_ad,2),'  RMS=',strtrim(rms_ad,2),' arcsec'
      end
      else:
    ENDCASE

    ;stop

    ;; Now match the sources with CAT
    ;MATCHSTARS,cat.x,cat.y,refcat.x,refcat.y,ind1,ind2,trans,count=nmatch,rms=rms,/silent
    ;
    ;; We have a X/Y MATCHING frame
    ;if (nmatch gt 10 and rms lt 1.0) then begin
    ;  print,ibase,'  Nmatch = ',strtrim(nmatch,2),'  RMS=',strtrim(rms,2),' pixels'
    ;  ndetected[ind1]++      ; increment the detected array
    ;
    ;; No MATCH with X/Y
    ;endif else begin
    ;
    ;  print,'No MATCHES with X/Y'
    ;
    ;  ; Do we have a WCS
    ;  maxra1 = max(cat.ra)
    ;  maxra2 = max(refcat.ra)
    ;  if maxra1 lt -1000.  and maxra2 lt -1000. then begin
    ;    print,'ONE of the images does NOT have a WCS.  NO MATCHES for ',ibase
    ;    goto,BOMB
    ;  endif
    ;
    ;  ; Try matching with the coordinates
    ;  print,'Checking RA/DEC matches'
    ;  rms2 = 999999.
    ;  SRCMATCH,cat.ra,cat.dec,refcat.ra,refcat.dec,1.0,ind1,ind2,count=nmatch2,/sph
    ;
    ;  ; We have an RA/DEC MATCH
    ;  if (nmatch2 gt 10) then begin
    ;
    ;    ; Calculate RMS
    ;    raresid = (cat[ind1].ra-refcat[ind2].ra)*cos(cat[ind1].dec/!radeg)*3600.
    ;    decresid = (cat[ind1].dec-refcat[ind2].dec)*3600.
    ;    resid = sqrt( raresid^2.0 + decresid^2.0 )
    ;    rms2 =  sqrt( mean( resid^2. ) )
    ;
    ;    print,ibase,'  Nmatch = ',strtrim(nmatch2,2),'  RMS=',strtrim(rms2,2),' arcsec'
    ;    ndetected[ind1]++      ; increment the detected array
    ;
    ;  ; NO MATCH WHATSOEVER
    ;  endif else begin
    ;    print,' NO MATCHES for ',ibase
    ;  endelse
    ;
    ;endelse ; no X/Y match

  endif ; we have coo/ap file

  BOMB:

end  ; fieldfiles loop


; Only keep sources that are seen in other frames
gd2 = where(ndetected ge 2,ngd2)  ; detected in TWO other frames
gd1 = where(ndetected ge 1,ngd1)  ; detected in ONE other frame


; Keeping detections in TWO other frames
if (ngd2 ge 100) then begin
  gd = gd2
  ngd = ngd2
  print,strtrim(ngd,2),' sources detected in at least TWO other frames'
endif
; Keeping detections in ONE other frame
if (ngd2 lt 100 and ngd1 gt 0) then begin
  gd = gd1
  ngd = ngd1
  print,strtrim(ngd,2),' sources detected in at least ONE other frame'
endif
; No detections in other frames
if (ngd1 eq 0) then begin
  print,'NO sources detected in other frames'
  error = 'NO sources detected in other frames'
  return
endif


; Fit Gaussians to the sources
if keyword_set(gaussfit) then begin
  print,'Fitting Gaussians to common sources'
  FITS_READ,file,im,head
  sz = size(im)
  x = lindgen(sz[1])
  y = lindgen(sz[2])
  gstr = REPLICATE({x:0.0,y:0.0,pars:fltarr(7),perror:fltarr(7),chisq:0.0,dof:0L,status:0L},ngd)
  For i=0,ngd-1 do begin

    ix = cat[gd[i]].x-1
    iy = cat[gd[i]].y-1
    xlo = (round(ix)-10)>0
    xhi = (round(ix)+10)<(sz[1]-1)
    ylo = (round(iy)-10)>0
    yhi = (round(iy)+10)<(sz[2]-1)
    subim = im[xlo:xhi,ylo:yhi]
    xarr = x[xlo:xhi]
    yarr = y[ylo:yhi]

    parinfo = replicate({limited:[0,0],limits:[0,0],fixed:0},7)
    parinfo[4].limited=[1,1]    ; constrain X and Y
    ;parinfo[4].limits=[-1,1]+ix
    parinfo[4].limits = [min(xarr),max(xarr)]
    parinfo[5].limited=[1,1]
    ;parinfo[5].limits=[-1,1]+iy
    parinfo[5].limits = [min(yarr),max(yarr)]
    ;estimates = []
    fit = MPFIT2DPEAK(subim,pars,xarr,yarr,chisq=chisq,dof=dof,perror=perror,/gaussian,$
                      parinfo=parinfo,status=status)
    gstr[i].x = ix
    gstr[i].y = iy
    gstr[i].pars = pars
    gstr[i].perror = perror
    gstr[i].chisq = chisq
    gstr[i].dof = dof
    gstr[i].status = status

    ; The 2D Gaussian parameters are:
    ;   A(0)   Constant baseline level
    ;   A(1)   Peak value
    ;   A(2)   Peak half-width (x) -- gaussian sigma or half-width at half-max
    ;   A(3)   Peak half-width (y) -- gaussian sigma or half-width at half-max
    ;   A(4)   Peak centroid (x)
    ;   A(5)   Peak centroid (y)
    ;   A(6)   Rotation angle (radians) if TILT keyword set

    ;display,subim,position=[0,0,0.5,1.0]
    ;display,fit,position=[0.5,0,1.0,1.0],/noerase
    ;wait,0.5
    ;stop

  End

  ; Now pick out the "good" ones
  medpar2 = MEDIAN(gstr.pars[2])
  sigpar2 = MAD(gstr.pars[2])
  medpar3 = MEDIAN(gstr.pars[3])
  sigpar3 = MAD(gstr.pars[3])
  medchisq = MEDIAN(gstr.chisq)
  sigchisq = MAD(gstr.chisq)
  ;  Positive amplitude,  Low Chisq, and Half-widths within "normal" values
  okay = where(gstr.pars[1] gt 0.0 AND gstr.chisq lt medchisq+3*sigchisq AND $
               abs(gstr.pars[2]-medpar2) lt 3*sigpar2 AND $
               abs(gstr.pars[3]-medpar3) lt 3*sigpar3,nokay)

  ; No stars passed, increase threshold, 4 sigma
  if nokay eq 0 then $
    okay = where(gstr.pars[1] gt 0.0 AND gstr.chisq lt medchisq+4*sigchisq AND $
                 abs(gstr.pars[2]-medpar2) lt 4*sigpar2 AND $
                 abs(gstr.pars[3]-medpar3) lt 4*sigpar3,nokay)
  ; No stars passed, increase threshold, 5 sigma
  if nokay eq 0 then $
    okay = where(gstr.pars[1] gt 0.0 AND gstr.chisq lt medchisq+5*sigchisq AND $
                 abs(gstr.pars[2]-medpar2) lt 5*sigpar2 AND $
                 abs(gstr.pars[3]-medpar3) lt 5*sigpar3,nokay)

  if nokay gt 0 then begin
    gd_orig = gd
    gd = gd[okay]
    print,strtrim(nokay,2),'/',strtrim(ngd,2),' sources passed the Gaussian fitting tests'
  endif

endif ; Gaussian fitting


; Create the COMMON file
; Only create it if there are enough sources
if (ngd ge minsources) then begin
  com = cat[gd]
  print,'Creating the CONFIRMED CELESTIAL SOURCES file = ',base+'.cmn.lst'
  WRITECOL,base+'.cmn.lst',com.id,com.x,com.y,com.mag,com.err,com.sky,com.skysig,$
           com.sharp,com.round,com.round2,fmt='(I7,2F9.2,3F9.3,F9.2,3F9.3)'

  ; Prepend the COO header
  WRITELINE,base+'.cmn.lst',[coohead1,''],/prepend

; Not enough sources
endif else begin
  print,'NO CMN.LST file created.  Need at least ',strtrim(minsources,2),' sources.'
endelse

if keyword_set(stp) then stop

end
