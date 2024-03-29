;+
;
; PHOTRED_FIELDSUMMARY
;
; Create summary files for each field
;
; INPUTS:
;  field      The short field name, i.e. F1, F2
;  =setupdir  The main directory
;  /redo      Redo
;
; OUTPUTS:
;  A summary field file
;
; USAGE:
;  IDL>photred_fieldsummary,'F1',setupdir=setupdir
;
; By D.Nidever  May 2015
;-

pro photred_fieldsummary,field,setupdir=setupdir,redo=redo,quick=quick
  
COMMON photred,setup

resolve_routine,'mrdfits',/compile_full_file,/either

undefine,error

; Not enough inputs
ninput = n_elements(field)
if ninput eq 0 then begin
  print,'Syntax - photred_fieldsummary,field,setupdir=setupdir,redo=redo'
  error = 'Not enough inputs'
  return
endif


; LOAD THE SETUP FILE if not passed
;-----------------------------------
; This is a 2xN array.  First colume are the keywords
; and the second column are the values.
; Use READPAR.PRO to read it
if n_elements(setup) eq 0 then begin
  PHOTRED_LOADSETUP,setup,setupdir=setupdir,count=count
  if count lt 1 then return
endif

; Log files
;----------
;  write to SAVE logfile
if n_elements(setupdir) gt 0 then logfile=setupdir+'/logs/SAVE.log' else $
  logfile = 'logs/SAVE.log'
logfile = FILE_EXPAND_PATH(logfile)  ; want absolute filename
if file_test(logfile) eq 0 then SPAWN,'touch '+logfile,out

; Telescope, Instrument
telescope = READPAR(setup,'TELESCOPE')
telescope = strupcase(telescope)
instrument = READPAR(setup,'INSTRUMENT')
instrument = strupcase(instrument)

; Get the scripts directory from setup
scriptsdir = READPAR(setup,'SCRIPTSDIR')
if scriptsdir eq '' then begin
  printlog,logfile,'NO SCRIPTS DIRECTORY'
  return
endif

; MCHUSETILES
mchusetiles = READPAR(setup,'MCHUSETILES')
if mchusetiles eq '0' or mchusetiles eq '' or mchusetiles eq '-1' then undefine,mchusetiles
tilesep = '+'

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

; Getting the transformation filename
transfile = READPAR(setup,'TRANS')
if transfile ne '0' and transfile ne '-1' then begin
  transfile = FILE_SEARCH(setupdir+'/'+transfile,/fully,count=ntransfile)
  if ntransfile gt 0 then begin
    ; Load the transformation equations
    READ_TRANS,transfile,trans,/silent

    ;; Do we have NIGHT information
    transnightinfo = 0
    if tag_exist(trans,'NIGHT') then begin
      gdtransnight = where(trans.night ge 0,ngdtransnight)
      if ngdtransnight gt 0 then transnightinfo=1
    endif
    ;; Do we have CHIP information
    transchipinfo = 0
    if tag_exist(trans,'CHIP') then begin
      gdtranschip = where(trans.chip ge 0,ngdtranschip)
      if ngdtranschip gt 0 then transchipinfo=1
    endif
    ;; Do we have FILE information
    transfileinfo = 0
    if tag_exist(trans,'FILE') then begin
      gdtransfile = where(trans.file ne '',ngdtransfile)
      if ngdtransfile gt 0 then transfileinfo=1
    endif
  endif
endif


printlog,logfile
printlog,logfile,'--- Making summary file for Field = '+field+' ---'
if keyword_set(quick) then printlog,logfile,'QUICK SUMMARY'
printlog,logfile


;#########################################
;#   Find all files for this field
;#########################################

;; Get files for this field from the success list
READLIST,setupdir+'/logs/DAOPHOT.success',fieldfiles,/unique,/fully,setupdir=setupdir,count=nfieldfiles,logfile=logfile,/silent
;; Get the files for this field
if nfieldfiles gt 0 then begin
  ; Get base names
  fbases = strarr(nfieldfiles)
  fext = strarr(nfieldfiles)
  for i=0,nfieldfiles-1 do begin
    if strmid(fieldfiles[i],6,7,/reverse_offset) eq 'fits.fz' then begin
      fbases[i] = FILE_BASENAME(fieldfiles[i],'.fits.fz')
      fext[i] = 'fits.fz'
    endif else begin
      fbases[i] = FILE_BASENAME(fieldfiles[i],'.fits')
      fext[i] = 'fits'
    endelse
  endfor
  ;; Get files only for this field
  flen = strlen(field)
  gd = where(strmid(fbases,0,flen+1) eq field+'-',nfieldfiles)
  if nfieldfiles gt 0 then begin
    fieldfiles = fieldfiles[gd]
    fbases = fbases[gd]
    fext = fext[gd]
  endif

  ; Make sure they are unique
  ;  sometimes you can accidentally get a fits and fits.fz files
  ;  for the same image
  if nfieldfiles gt 0 then begin
    dbl = doubles(fbases,count=ndbl)
    dblall = doubles(fbases,/all)
    badmask = bytarr(nfieldfiles)
    for i=0,ndbl-1 do begin
      MATCH,fbases[dblall],fbases[dbl[i]],ind1,ind2,/sort,count=nmatch
      ind = dblall[ind1]
      badmask[ind] = 1
      ; keep the fits.fz ending over the fits one
      ;  else take the first fits file
      gd = where(fext[ind] eq 'fits.fz',ngd)
      if ngd gt 0 then badmask[ind[gd[0]]]=0 else badmask[ind[0]]=0     
    endfor
    bad = where(badmask eq 1,nbad)
    if nbad gt 0 then REMOVE,bad,fieldfiles,fbases,fext
    nfieldfiles = n_elements(fieldfiles)
  endif
endif

if nfieldfiles eq 0 then begin
  printlog,logfile,'No ',field,' files found in current directory'
  return
endif

printlog,logfile,'Found ',strtrim(nfieldfiles,2),' frames of FIELD=',field

; Getting names of final/dat/dered files
finalfile = file_search('*.final',count=nfinalfile)  ; e.g. Field57sh.final
;; if there are multiple, then take the most recent one
if nfinalfile gt 0 then begin
  info = file_info(finalfile)
  si = reverse(sort(info.mtime))  ;; most recent one first
  finalfile = finalfile[si[0]]
endif
finalname = file_basename(finalfile,'.final')
datfile = finalname+'.dat'
deredfile = file_search(field+'-*.dered',count=ndered)
refbase = file_basename(deredfile,'.dered')    ; e.g. F1-00277654.dered

; Check if the output file already exists
outfile = finalname+'_summary.fits'
if file_test(outfile) eq 1 and not keyword_set(redo) then begin
  printlog,logfile,outfile,' already EXISTS and /redo not set'
  return
endif

; Load the FINAL structure for this field
printlog,logfile,'Loading Final structure'
if not keyword_set(quick) then begin
  if file_test(datfile) eq 1 then restore,datfile else $
    final=importascii(finalfile,/header)
endif

; Load the aperture correction "apcor.lst" file
undefine,apcor
apcorfile = setupdir+'/apcor.lst'
if file_test(apcorfile) eq 1 then begin
  ; Getting the aperture correction structure
  apcor = IMPORTASCII(apcorfile,fieldnames=['name','value'],/noprint)
  ; Remove the 'a.del' endings for the names
  apcor_orig = apcor
  apcor.name = repstr(apcor.name,'a.del','')  ; base names
endif

;; Load all the ALS files in the MCH files
;;  this makes a difference if we are using tiles
if keyword_set(mchusetiles) then begin
  READLIST,setupdir+'/logs/ASTROM.success',mchfiles,/unique,/fully,setupdir=setupdir,count=nmchfiles,logfile=logfile,/silent
  ;; Replace any .mag with .mch
  mchfiles = repstr(mchfiles,'.mag','.mch')
  ;; Get files for this FIELD
  flen = strlen(field)
  gd = where(strmid(file_basename(mchfiles),0,flen+1) eq field+'-',ngd)
  if ngd gt 0 then begin
    mchfiles = mchfiles[gd]
    nmchfiles = n_elements(mchfiles)
    ;; Load all of the ALS files from the MCH files
    undefine,alsinfo
    for i=0,nmchfiles-1 do begin
      LOADMCH,mchfiles[i],alsfiles,count=nalsfiles
      alsinfo1 = replicate({dir:'',mchfile:'',alsfile:''},nalsfiles)
      alsinfo1.dir = file_dirname(mchfiles[i])
      alsinfo1.mchfile = file_basename(mchfiles[i])
      alsinfo1.alsfile = alsfiles
      push,alsinfo,alsinfo1
    endfor
 endif else printlog,logfile,'No ALS files in the MCH files for this field'
endif else printlog,logfile,'No MCH tiles for this field'

; Loop through all files for this field
;  and gather all of the necessary information
nan = !values.f_nan
dnan = !values.d_nan
chipstr = replicate({field:'NAN',file:'NAN',expnum:'NAN',chip:-1L,base:'NAN',filter:'NAN',exptime:nan,utdate:'NAN',uttime:'NAN',$
                     airmass:nan,gain:nan,rdnoise:nan,nx:-1L,ny:-1L,wcstype:'NAN',pixscale:nan,ra:dnan,dec:dnan,wcsrms:nan,fwhm:nan,$
                     skymode:nan,skysig:nan,dao_nsources:-1L,dao_depth:nan,dao_npsfstars:-1L,dao_psftype:'NAN',dao_psfboxsize:-1L,$
                     dao_psfvarorder:-1L,dao_psfchi:nan,alf_nsources:-1L,alf_depth:nan,calib_depth:nan,calib_color:'NAN',calib_zpterm:nan,$
                     calib_zptermsig:nan,calib_amterm:nan,calib_amtermsig:nan,calib_colorterm:nan,calib_colortermsig:nan,$
                     calib_magname:'NAN',apcor:nan,ebv:nan},nfieldfiles)
printlog,logfile,''
printlog,logfile,'Chip-level information'
printlog,logfile,''
printlog,logfile,'  NUM       Filename    Filt  Exptime WCSRMS   FWHM   Skymode  Skysig  DAO_Nsrc DAO_Depth PSFtype NPSFstars PSFchi ALF_Nsrc APcor'
For i=0,nfieldfiles-1 do begin
  undefine,fitsfile,base,shfield,expnum,chip,filter
  undefine,exptime,utdate,uttime,airmass,gain,rdnoise
  undefine,head,nx,ny,ctype1,astr,scale,ra,dec,wcsrms
  undefine,optlines,optarr,loglines,skymode,skysig
  undefine,als,alshead,hist,xhist,minarr,maxarr,alsdepth
  undefine,lstlines,psflines,psfarr,psfva,psfloglines
  undefine,chilines,chilinesarr,chiarr,psfan,minpsftype,maxpsftype,psfchi

  fitsfile = fieldfiles[i]
  fitsdir = file_dirname(fitsfile)+'/'
  if strmid(fitsfile,6,7,/reverse_offset) eq 'fits.fz' then fpack=1 else fpack=0
  if fpack eq 1 then base=file_basename(fitsfile,'.fits.fz') else $
    base = file_basename(fitsfile,'.fits')   
  shfield = first_el(strsplit(base,'-',/extract))
  if thisimager.namps gt 1 then begin
    tmp = first_el(strsplit(base,'-',/extract),/last)
    expnum = first_el(strsplit(tmp,thisimager.separator,/extract))
    chip = fix( first_el(strsplit(tmp,thisimager.separator,/extract),/last) )
  endif else begin
    chip = 1
    expnum = first_el(strsplit(base,'-',/extract),/last)
  endelse
  chipstr[i].field = field
  chipstr[i].file = fitsfile
  chipstr[i].expnum = expnum
  chipstr[i].chip = chip
  chipstr[i].base = base

  ;; Read in the header
  if fpack eq 1 then begin
    head = PHOTRED_READFILE(fitsfile,exten=1,/header)
    ; Fix the NAXIS1/NAXIS2 in the header
    sxaddpar,head,'NAXIS1',sxpar(head,'ZNAXIS1')
    sxaddpar,head,'NAXIS2',sxpar(head,'ZNAXIS2')
  endif else begin
    head = PHOTRED_READFILE(fitsfile,/header)
  endelse

  ; Filter, Exptime, Date/Time, airmass, gain, rdnoise
  filter = PHOTRED_GETFILTER(fitsfile,head=head)
  chipstr[i].filter = filter
  exptime = PHOTRED_GETEXPTIME(fitsfile,head=head)
  chipstr[i].exptime = exptime
  utdate = PHOTRED_GETDATE(fitsfile,head=head)
  chipstr[i].utdate = utdate
  uttime = PHOTRED_GETUTTIME(fitsfile,head=head)
  chipstr[i].uttime = uttime  
  airmass = PHOTRED_GETAIRMASS(fitsfile,head=head)
  chipstr[i].airmass = airmass
  gain = PHOTRED_GETGAIN(fitsfile,head=head)
  chipstr[i].gain = gain
  rdnoise = PHOTRED_GETRDNOISE(fitsfile,head=head)
  chipstr[i].rdnoise = rdnoise

  ; From FITS header
  ;  nx, ny, ctype, scale, ra, dec
  nx = sxpar(head,'NAXIS1',count=n_nx)
  if n_nx gt 0 then chipstr[i].nx = nx
  ny = sxpar(head,'NAXIS2',count=n_ny)
  if n_ny gt 0 then chipstr[i].ny = ny
  ctype1 = sxpar(head,'CTYPE1',count=n_ctype1)
  if n_ctype1 gt 0 then begin
    WCS_CHECK_CTYPE,ctype1,wcstype
    chipstr[i].wcstype = wcstype
  endif
  EXTAST,head,astr,noparams
  if noparams ge 1 then begin
    GETPIXSCALE,fitsfile,scale,head=head
    chipstr[i].pixscale = scale
    head_xyad,head,nx/2,ny/2,ra,dec,/deg
    chipstr[i].ra = ra
    chipstr[i].dec = dec
  endif

  ; WCS RMS
  wcsind = where(stregex(head,'WCSFIT: RMS',/boolean) eq 1,nwcsfit)
  if nwcsfit gt 0 then begin
    ; HISTORY WCSFIT: RMS=0.216 arcsec on Fri Apr 24 10:20:29 2015
    wcsline = head[wcsind[0]]
    lo = strpos(wcsline,'RMS=')
    tmp = strmid(wcsline,lo+4)
    wcsrms = float( first_el(strsplit(tmp,' ',/extract)) )
    chipstr[i].wcsrms = wcsrms
  endif
    
  ; Load DAOPHOT option file
  optfile = fitsdir+base+'.opt'
  if file_test(optfile) eq 1 then begin
    READLINE,optfile,optlines
    optarr = strsplitter(optlines,'=',/extract)
    ind = where(strtrim(optarr[0,*],2) eq 'FW',nind)
    if nind gt 0 then chipstr[i].fwhm=float(optarr[1,ind[0]])
  endif
     
  ; Load DAOPHOT log file
  daologfile = fitsdir+base+'.log'
  if file_test(daologfile) eq 1 then begin
    READLINE,daologfile,loglines
    ; Sky mode and standard deviation =   48.608    4.110
    ind = where(stregex(loglines,'Sky mode',/boolean,/fold_case) eq 1,nind)
    if nind gt 0 then begin
      line = loglines[ind[0]]
      lo = strpos(line,'=')
      skymode = float( first_el(strsplit(strmid(line,lo+1),' ',/extract)) )
      skysig = float( first_el(strsplit(strmid(line,lo+1),' ',/extract),/last) )
      chipstr[i].skymode = skymode
      chipstr[i].skysig = skysig
    endif
  endif
  ;; Get the sky level from the ALS below instead,
  ;;    much faster if I need to recreate the FITS files
  ;if n_elements(skymode) eq 0 then begin 
  ;  im = PHOTRED_READFILE(fitsfile,head)
  ;  SKY,im,skymode,skysig,/silent
  ;  chipstr[i].skymode = skymode
  ;  chipstr[i].skysig = skysig
  ;endif

  ; Load ALS file
  alsfile = fitsdir+base+'.als'
  if file_test(alsfile) eq 1 then begin
    if not keyword_set(quick) or n_elements(skymode) eq 0 then begin
      LOADALS,alsfile,als,alshead,count=nals
      chipstr[i].dao_nsources = nals
      ;; Calculate sky values
      ;;   much faster than reading from FITS if I need to recreate the FITS files 
      if n_elements(skymode) eq 0 then begin
        chipstr[i].skymode = median([als.sky])
        chipstr[i].skysig = mad([als.sky])
      endif
      ;; Calculate "depth"
      if nals gt 0 then begin
        hist = histogram(als.mag,bin=0.2,locations=xhist,min=0,max=30)
        xhist += 0.5*0.2
        DLN_MAXMIN,hist,minarr,maxarr
        alsdepth = xhist[first_el(maxarr,/last)]  ; use last maximum
        chipstr[i].dao_depth = alsdepth   ; instrumental "depth"    
      endif
    endif else chipstr[i].dao_nsources = file_lines(alsfile)-3
  endif

  ; Check list of DAOPHOT PSF stars
  lstfile = fitsdir+base+'.plst'
  if file_test(lstfile) eq 1 then begin
    READLINE,lstfile,lstlines
    chipstr[i].dao_npsfstars = n_elements(lstlines)-3
  endif

  ; Load DAOPHOT PSF file
  psffile = fitsdir+base+'.psf'
  if file_test(psffile) eq 1 then begin
    READLINE,psffile,psflines
    if n_elements(psflines) gt 0 then begin
      ; PENNY1     69    4    6    0   14.048       1091.621   1022.5   2046.5
      psfarr = strsplit(psflines[0],' ',/extract)
      chipstr[i].dao_psftype = psfarr[0]
      chipstr[i].dao_psfboxsize = long(psfarr[1])
      ; Get PSF spatial variation from OPT file (VA)
      if n_elements(optarr) gt 0 then begin
        indva = where(strtrim(optarr[0,*],2) eq 'VA',nindva)
        if nindva gt 0 then psfva=long(optarr[1,indva[0]])
      endif
      ; Get PSF spatial variation from PSF file
      if n_elements(psfva) eq 0 then begin
        ; NEXP=1 (VA=0), NEXP=3 (VA=1), NEXP=6 (VA=2)
        nexp = long(psfarr[3])
        psfva = nexp/3 
      endif
      if n_elements(psfva) gt 0 then chipstr[i].dao_psfvarorder=psfva
    endif else begin
      ;; no PSF lines
      printlog,logfile,psffile,' IS EMPTY'
      chipstr[i].dao_psftype = -1
      chipstr[i].dao_psfboxsize = -1
      chipstr[i].dao_psfvarorder = -99
    endelse
  endif

  ; PSF chi-value
  psflogfile = fitsdir+base+'.psf.log'
  if file_test(psflogfile) eq 1 then begin
    READLINE,psflogfile,psfloglines
    ;
    ;Chi    Parameters...
    ;>>   0.0261   2.04305   2.13328
    ;>>   0.0666   1.36261   1.30098  -0.33620
    ;>>   0.0227   1.78936   1.82379  -0.05575
    ;>>   0.0161   1.91927   2.01636  -0.03596
    ;>>   0.0587   1.71325   1.82209  -0.05311
    ;>>   0.0151   1.99451   2.09012   0.81983  -0.03768
    ;
    ;
    ;Profile errors:
    ;
    ;46  0.020      1130  0.018      2294  0.017      3055  0.017      4279  0.012
    ;49  0.019      1224  0.012      2375  0.015      3063  0.014      4280  0.011
    ;75  0.013      1246  0.015      2397  0.011      3169  0.013      4354  0.012
    ind1 = where(stregex(psfloglines,'Chi',/boolean) eq 1 and stregex(psfloglines,'Parameters',/boolean) eq 1,nind1)
    ind2 = where(stregex(psfloglines,'Profile errors',/boolean) eq 1,nind2)
    if nind1 gt 0 and nind2 gt 0 then begin
      chilines = psfloglines[ind1[0]+1:ind2[0]-1]
      ; Replace "Failed to converge" with high chi
      bdchilines = where(stregex(chilines,'Failed to converge',/fold_case,/boolean) eq 1,nbdchilines)
      if nbdchilines gt 0 then chilines[bdchilines]='>>  99.99'
      gdchilines = where(strtrim(chilines,2) ne '' and stregex(chilines,'>>',/boolean) eq 1,ngdchilines)
      if ngdchilines gt 0 then chilines=chilines[gdchilines] else undefine,chilines
      chilinesarr = strsplitter(strmid(chilines,3),' ',/extract)
      chiarr = float(reform(chilinesarr[0,*]))
      ; If AN (Analytic model PSF) is negative, then try all PSF types
      ; from 1 to |AN|, inclusive.  PHOTRED default is "-6"
      ; from NPARAM in mathsubs.f
      ; LABEL     IPSTYPE  NPARAM
      ; GAUSSIAN    1         2
      ; MOFFAT15    2         3
      ; MOFFAT25    3         3
      ; MOFFAT35    4         3
      ; LORENTZ     5         3
      ; PENNY1      6         4
      ; PENNY2      7         5
      ; subroutine GETPSF in psf.f loops from MAX(1,AN) to |AN|, so for
      ;  PHOTRED it should be 1 to 6
      psflabels = ['GAUSSIAN','MOFFAT15','MOFFAT25','MOFFAT35','LORENTZ','PENNY1','PENNY2']
      psftype = [1,2,3,4,5,6,7]
      if n_elements(optarr) gt 0 then begin
        indpsfan = where(strtrim(optarr[0,*],2) eq 'AN',nindpsfan)
        psfan = long(optarr[1,indpsfan[0]])
        minpsftype = max([1,psfan])
        maxpsftype = abs(psfan)
      endif else begin  ; no opt file assuming AN=-6
        minpsftype = 1
        maxpsftype = 6
      endelse
      ; Get the chi for the right LABEL that was used
      psflabelused = chipstr[i].dao_psftype
      indpsf = first_el(where(psflabels eq psflabelused,nindpsf))
      psftypeused = psftype[indpsf]
      ;if (psftypeused-minpsftype+1) gt n_elements(chiarr) then psfchi=min(chiarr) else $
      ;  psfchi = chiarr[psftypeused-minpsftype]
      ;if abs(psfchi-min(chiarr)) gt 0.01 then stop,'Problem with PSF chi'
      psfchi = min(chiarr) ; just use the MINIMUM chi value!
      chipstr[i].dao_psfchi = psfchi
    endif
  endif

  ; Load ALF file
  nalf = -1
  if keyword_set(mchusetiles) eq 0 then begin
    alffile = fitsdir+base+'.alf'
    if file_test(alffile) eq 1 then begin
      if not keyword_set(quick) then LOADALS,alffile,alf,alfhead,count=nalf else nalf=file_lines(alffile)-3
    endif
  ;; Using TILES
  endif else begin
    gd = where(stregex(alsinfo.alsfile,'^'+base,/boolean) eq 1,ngd)
    undefine,alf
    if keyword_set(quick) then nalf=0
    for j=0,ngd-1 do begin
      alffile1 = alsinfo[gd[j]].dir+'/'+repstr(alsinfo[gd[j]].alsfile,'.als','.alf')
      if file_test(alffile1) eq 1 then begin
        if not keyword_set(quick) then begin
          LOADALS,alffile1,alf1,alfhead1,count=nalf1 
          if nalf1 gt 0 then PUSH,alf,alf1
        endif else nalf+=file_lines(alffile1)-3
      endif
    endfor
    if not keyword_set(quick) then nalf = n_elements(alf)
  endelse
  if nalf ge 0 then chipstr[i].alf_nsources = nalf

  ; Calculate ALF "depth"
  if nalf gt 1 and not keyword_set(quick) then begin
    hist = histogram(alf.mag,bin=0.2,locations=xhist,min=0,max=30)
    xhist += 0.5*0.2
    DLN_MAXMIN,hist,minarr,maxarr
    alfdepth = xhist[first_el(maxarr,/last)]  ; use last maximum
    chipstr[i].alf_depth = alfdepth   ; instrumental "depth"    
  endif

  ; Get aperture correction
  if n_elements(apcor) gt 0 then begin
    indapcor = where(apcor.name eq base,nindapcor)    
    if nindapcor gt 0 then chipstr[i].apcor=apcor[indapcor[0]].value
  endif

  printlog,logfile,i+1,base,chipstr[i].filter,chipstr[i].exptime,chipstr[i].wcsrms,chipstr[i].fwhm,chipstr[i].skymode,chipstr[i].skysig,chipstr[i].dao_nsources,chipstr[i].dao_depth,$
           chipstr[i].dao_psftype,chipstr[i].dao_npsfstars,chipstr[i].dao_psfchi,chipstr[i].alf_nsources,chipstr[i].apcor,format='(I5,A18,A5,F8.2,2F8.4,F8.2,F8.3,I9,F9.2,A11,I8,F8.3,I8,F8.3)'
Endfor

;; QUICK, get median ALS depth
g = where(chipstr.dao_nsources gt 0,ng)
if ng gt 0 then begin
  ;; Get the median chip file in DAO_NSOURCES
  ind = where(chipstr[g].dao_nsources eq median([chipstr[g].dao_nsources]))
  bestind = g[ind[0]]
  LOADALS,file_dirname(chipstr[bestind].file)+'/'+chipstr[bestind].base+'.als',als,alshead,count=nals
  ; Calculate "depth"
  if nals gt 0 then begin
    hist = histogram(als.mag,bin=0.2,locations=xhist,min=0,max=30)
    xhist += 0.5*0.2
    DLN_MAXMIN,hist,minarr,maxarr
    alsdepth = xhist[first_el(maxarr,/last)]  ; use last maximum
    chipstr[g].dao_depth = alsdepth   ; instrumental "depth"    
    printlog,logfile,'Median DAO depth = '+stringize(alsdepth,ndec=3)+' mag'
  endif  
endif
;; QUICK, get median ALF depth
g = where(chipstr.alf_nsources gt 0,ng)
if ng gt 0 then begin
  ;; Get the median chip file in DAO_NSOURCES
  ind = where(chipstr[g].alf_nsources eq median(chipstr[g].alf_nsources))
  bestind = g[ind[0]]
  fitsdir = file_dirname(chipstr[bestind].file)
  base = chipstr[bestind].base
  nalf = -1
  if keyword_set(mchusetiles) eq 0 then begin
    alffile = fitsdir+base+'.alf'
    if file_test(alffile) eq 1 then LOADALS,alffile,alf,alfhead,count=nalf
  ;; Using TILES
  endif else begin
    gd = where(stregex(alsinfo.alsfile,'^'+base,/boolean) eq 1,ngd)
    undefine,alf
    for j=0,ngd-1 do begin
      alffile1 = alsinfo[gd[j]].dir+'/'+repstr(alsinfo[gd[j]].alsfile,'.als','.alf')
      if file_test(alffile1) eq 1 then begin
        LOADALS,alffile1,alf1,alfhead1,count=nalf1 
        if nalf1 gt 0 then PUSH,alf,alf1
      endif
    endfor
    nalf = n_elements(alf)
  endelse
  ; Calculate "depth"
  if nalf gt 0 then begin
    hist = histogram(alf.mag,bin=0.2,locations=xhist,min=0,max=30)
    xhist += 0.5*0.2
    DLN_MAXMIN,hist,minarr,maxarr
    alfdepth = xhist[first_el(maxarr,/last)]  ; use last maximum
    chipstr[g].alf_depth = alfdepth   ; instrumental "depth"    
    printlog,logfile,'Median ALF depth = '+stringize(alfdepth,ndec=3)+' mag'
  endif  
endif


; Get calibrated depth and transformation equation on a chip-basis
; ------------------------------------------------------------------
;  Use the .phot calibrated photometry file, photcalib_prep .input files,
;  and .mch files to get the calibrated photometry for each chip.

printlog,logfile,' '
printlog,logfile,'Getting calibrated photometry from the .phot files'

;; Loop over the groups/phot files
READLIST,setupdir+'/logs/COMBINE.success',photfiles,/unique,/fully,setupdir=setupdir,count=nphotfiles,logfile=logfile,/silent
flen = strlen(field)
gd = where(strmid(file_basename(photfiles),0,flen+1) eq field+'-',ngd)
if ngd gt 0 then begin
  photfiles = photfiles[gd]
  nphotfiles = ngd
endif else begin
  printlog,logfile,'NO phot files for field '+field
  undefine,photfiles
  nphotfiles = 0
endelse
For i=0,nphotfiles-1 do begin
  ;; Load the phot and mch files
  if not keyword_set(quick) then begin
    phot = PHOTRED_READFILE(photfiles[i],meta,count=nphot)
  endif else begin
    if file_isfits(photfiles[i]) eq 1 then begin
      meta = PHOTRED_READFILE(photfiles[i],/header)
      ;; use mrd_table.pro to load a blank table
      mrd_table,meta,structyp,0,[0,9],rsize,phot,nrows,nfld,typarr,fnames,$
                fvalues,vcls, vtpes, scales, offsets, 0, status,/silent
      nphot = sxpar(meta,'naxis2')
    endif else begin
      phot = importascii(photfiles[i],/header,nlineread=10,count=nphot,/silent)
    endelse
  endelse
  phtags = tag_names(phot)
  mchfile = repstr(photfiles[i],'.phot','.mch')
  LOADMCH,mchfile,alsfiles,transmch,count=nalsfiles
  inputfile = repstr(photfiles[i],'.phot','.input')
  READLINE,inputfile,inputlines
  inputarr = strsplit(inputlines[0],' ',/extract)
  ;   Old or New format, check if the second value is
  ;     an integer (new format, NIGHT) or character (old format,
  ;     band/filter name)
  ;    All lines need to use the same format (old or new)
  newformat = valid_num(inputarr[1],/integer)
  ; -- NEW Format --
  ; photometry filename, NIGHT, CHIP, filter, airmass, exptime, aperture
  ;   correction, Band2 ... 
  if newformat then begin
    alsfilter = inputarr[3:*:6]  ; filters
    night = inputarr[1:*:6]      ; MJD 58766
  ; -- OLD Format --
  ; photometry filename, Band1 name, Band1 airmass, Band1 exptime, Band1
  ;   aperture correction, Band2 ... 
  endif else begin
    alsfilter = inputarr[1:*:4]  ; filters
  endelse

  printlog,logfile,i+1,file_basename(photfiles[i]),nphot,format='(I5,A22,I8)'

  ; Creating the column names (code from photcalib.pro)
  ui = uniq(alsfilter,sort(alsfilter))
  ubands = alsfilter[ui]
  nubands = n_elements(ubands)
  magname = strarr(n_elements(alsfiles))
  for j=0,nubands-1 do begin
    gdbands = where(alsfilter eq ubands[j],ngdbands)
    ; More than one observation in this band
    if (ngdbands gt 1) then begin
      magname[gdbands] = strupcase(ubands[j])+'MAG'+strtrim(indgen(ngdbands)+1,2)
    ; Only ONE obs in this band
    endif else begin
      magname[gdbands[0]] = strupcase(ubands[j])+'MAG'
    endelse
  endfor
  ; But these column/mag names are only correct if it
  ;  saved the individual exposure magnitudes

  ;; Get the CHIPSTR indices for the als files for this group
  if keyword_set(mchusetiles) then begin
    alsbase = reform((strsplitter(alsfiles,tilesep+'T',/extract))[0,*])
  endif else alsbase=file_basename(alsfiles,'.als')
  MATCH,chipstr.base,alsbase,indgrp,alsind,/sort,count=nindgrp

  ; Loop through the matched images for this group
  for j=0,nindgrp-1 do begin
    indgrp1 = indgrp[j]
    alsind1 = alsind[j]
    imagname = magname[alsind1]
    phind = where(phtags eq imagname,nphind)
    if nphind gt 0 and not keyword_set(quick) then begin
      ; Calculate "depth"
      hist = histogram(phot.(phind[0]),bin=0.2,locations=xhist,min=0,max=50)
      xhist += 0.5*0.2
      DLN_MAXMIN,hist,minarr,maxarr
      depth = xhist[first_el(maxarr,/last)] ; use last maximum
      chipstr[indgrp1].calib_depth = depth   ; calibrated "depth"
      chipstr[indgrp1].calib_magname = imagname  ; keep track of magname
      ;; CALIB_MAGNAME will not be unique if using TILES!!!

      ; Get EBV for these stars
      gdmag = where(phot.(phind[0]) lt 50,ngdmag)
      if ngdmag gt 0 then begin
        SRCMATCH,final.ra,final.dec,phot[gdmag].ra,phot[gdmag].dec,0.2,ind1,ind2,/sph,count=nmatch
        if nmatch gt 0 then begin
          med_ebv = median([final[ind1].ebv],/even)
          chipstr[indgrp1].ebv = med_ebv
        endif
      endif
    endif ; we have the proper column/tag
    ;; Get rough E(B-V) quickly
    if keyword_set(quick) then begin
      glactc,chipstr[indgrp1].ra,chipstr[indgrp1].dec,2000.0,glon,glat,1,/deg
      ebv = dust_getval(glon,glat,/noloop,/interp)
      chipstr[indgrp1].ebv = median([ebv])
    endif
  Endfor  ; images in group loop

  ; Fill in transformation equation info
  if n_elements(trans) gt 0 then begin
    ; Loop through the images in this group
    For j=0,nindgrp-1 do begin
      indgrp1 = indgrp[j]
      alsind1 = alsind[j]

      ;; There are four options for matching the file:
      ;; 1) No chip/night/filename information
      ;; 2) Only chip given
      ;; 3) Night+chip given
      ;; 4) Filename given
      ;; This information can be different for each file, i.e. multiple
      ;; formats can be used in a given trans file.  Try MOST specific
      ;; (#4) to LEAST specific (#1).
      nmatch = 0
      ;; Try filename + band
      if transfileinfo eq 1 then begin
         ;; First check for an exact match
         MATCH,trans.file+':'+trans.band,chipstr[indgrp1].base+':'+chipstr[indgrp1].filter,ind1,ind2,/sort,count=nmatch
         ;; Check if FILE is in FIELD-EXPOSURE format and match to all chips of this exposure
         if nmatch eq 0 then begin
           match = bytarr(n_elements(trans))
           for k=0,n_elements(trans)-1 do $
             match[k] = (stregex(chipstr[indgrp1].base,trans[k].file,/boolean) eq 1 and trans[k].file ne '' and trans[k].band eq chipstr[indgrp1].filter)
           ind1 = where(match eq 1,nmatch)
         endif
      endif
      mjd = PHOTRED_GETMJD('',thisimager.observatory,dateobs=chipstr[indgrp1].utdate+'T'+chipstr[indgrp1].uttime,error=errormjd)
      ;; Try night+chip + band
      ;;  only want to match lines with file=''
      if nmatch eq 0 and transchipinfo eq 1 and transnightinfo eq 1 then $
         MATCH,trans.file+':'+strtrim(trans.night,2)+':'+strtrim(trans.chip,2)+':'+trans.band,$
               ':'+strtrim(mjd,2)+':'+strtrim(chipstr[indgrp1].chip,2)+':'+chipstr[indgrp1].filter,ind1,ind2,/sort,count=nmatch
      ;; Try chip + band
      ;;  only want to match lines with file='' and night=-1
      if nmatch eq 0 and transchipinfo eq 1 then $
         MATCH,trans.file+':'+strtrim(trans.night,2)+':'+strtrim(trans.chip,2)+':'+trans.band,$
               ':-1:'+strtrim(chipstr[indgrp1].chip,2)+':'+chipstr[indgrp1].filter,ind1,ind2,/sort,count=nmatch
      ;; Try just the band
      ;;  only want to match lines with file='', night=-1 and chip=-1
      if nmatch eq 0 then $
         MATCH,trans.file+':'+strtrim(trans.night,2)+':'+strtrim(trans.chip,2)+':'+trans.band,$
               ':-1:-1:'+chipstr[indgrp1].filter,ind1,ind2,/sort,count=nmatch

      ;; Found the transformation for this file+band
      if (nmatch gt 0) then begin
        trans1 = trans[ind1[0]]
      endif else begin
        printlog,logf,'NO TRANSFORMATION INPUT FOR  ',chipstr[indgrp1].file
        return
      endelse
      chipstr[indgrp1].calib_color = trans1.color
      chipstr[indgrp1].calib_zpterm = trans1.zpterm
      chipstr[indgrp1].calib_zptermsig = trans1.zptermsig
      chipstr[indgrp1].calib_amterm = trans1.amterm
      chipstr[indgrp1].calib_amtermsig = trans1.amtermsig
      chipstr[indgrp1].calib_colorterm = trans1.colterm
      chipstr[indgrp1].calib_colortermsig = trans1.coltermsig
    Endfor
  endif 

  ; Get EBV another way if necessary
  bdebv = where(finite(chipstr[indgrp].ebv) eq 0,nbdebv)
  for j=0,nbdebv-1 do begin
    fitsfile = chipstr[indgrp[bdebv[j]]].file
    fitstest = file_test(fitsfile)
    alsfile = chipstr[indgrp[bdebv[j]]].base+'.als'
    alstest = file_test(alsfile)
    coofile = chipstr[indgrp[bdebv[j]]].base+'.coo'
    cootest = file_test(coofile)
    if fitstest eq 1 and (alstest eq 1 or cootest eq 1) then begin
      if alstest eq 1 then LOADALS,alsfile,cat else $
        LOADCOO,coofile,cat
      if strmid(fitsfile,6,7,/reverse_offset) eq 'fits.fz' then begin
        head = PHOTRED_READFILE(fitsfile,exten=1,/header)
        ; Fix the NAXIS1/NAXIS2 in the header
        sxaddpar,head,'NAXIS1',sxpar(head,'ZNAXIS1')
        sxaddpar,head,'NAXIS2',sxpar(head,'ZNAXIS2')
      endif else begin
        head = PHOTRED_READFILE(fitsfile,/header)
      endelse
      HEAD_XYAD,head,cat.x-1,cat.y-1,ra,dec,/deg
      SRCMATCH,final.ra,final.dec,ra,dec,0.2,ind1,ind2,/sph,count=nmatch
      if nmatch gt 0 then begin
        med_ebv = median([final[ind1].ebv],/even)
        chipstr[indgrp[bdebv[j]]].ebv = med_ebv
      endif
    endif
  endfor  ; bad EBV loop
Endfor  ; group loop


; Create the exposure level information structure
;------------------------------------------------
printlog,logfile,''
ui = uniq(chipstr.expnum,sort(chipstr.expnum))
expnum = chipstr[ui].expnum
nexpnum = n_elements(expnum)
printlog,logfile,'There are '+strtrim(nexpnum,2),' unique exposures'
printlog,logfile,'Exposure-level information'

; Loop through the exposures
expstr = replicate({expnum:'',nchips:-1L,filter:'NAN',exptime:nan,utdate:'',uttime:'',airmass:nan,wcstype:'NAN',ra:dnan,dec:dnan,wcsrms:nan,fwhm:nan,skymode:nan,skysig:nan,$
                    dao_nsources:-1L,dao_depth:nan,dao_psfchi:nan,alf_nsources:-1L,alf_depth:nan,apcor:nan,ebv:nan,magname:'NAN'},nexpnum)
printlog,logfile,''
printlog,logfile,' NUM    EXPNUM   Filter Exptime        DATE/TIME          WCSRMS   FWHM   Skymode Skysig  DAO_Nsrc DAO_Depth PSFchi ALF_Nsrc ALF_Depth APcor    EBV'
For i=0,nexpnum-1 do begin
  iexpnum = expnum[i]
  ind = where(chipstr.expnum eq iexpnum,nind)
  chipstr1 = chipstr[ind]
  expstr[i].expnum = iexpnum
  expstr[i].nchips = nind
  expstr[i].filter = chipstr1[0].filter
  expstr[i].exptime = chipstr1[0].exptime
  expstr[i].utdate = chipstr1[0].utdate
  expstr[i].uttime = chipstr1[0].uttime
  gd = where(finite(chipstr1.airmass) eq 1,ngd)
  if ngd gt 0 then expstr[i].airmass = median([chipstr1[gd].airmass],/even)
  expstr[i].wcstype = chipstr1[0].wcstype
  gd = where(finite(chipstr1.ra) eq 1,ngd)
  if ngd gt 0 then begin
    ra = chipstr1[gd].ra
    if range(ra) gt 100 then begin
      bd = where(ra gt 180,nbd)
      if nbd gt 0 then ra[bd] -= 360
    endif
    mnra = mean(minmax(ra))
    if mnra lt 0 then mnra+=360
    expstr[i].ra = mnra
  endif
  gd = where(finite(chipstr1.dec) eq 1,ngd)
  if ngd gt 0 then expstr[i].dec = median([chipstr1[gd].dec],/even)
  gd = where(finite(chipstr1.wcsrms) eq 1,ngd)
  if ngd gt 0 then expstr[i].wcsrms = median([chipstr1[gd].wcsrms],/even)
  gd = where(finite(chipstr1.fwhm) eq 1,ngd)
  if ngd gt 0 then expstr[i].fwhm = median([chipstr1[gd].fwhm],/even)
  gd = where(finite(chipstr1.skymode) eq 1,ngd)
  if ngd gt 0 then expstr[i].skymode = median([chipstr1[gd].skymode],/even)
  gd = where(finite(chipstr1.skysig) eq 1,ngd)
  if ngd gt 0 then expstr[i].skysig = median([chipstr1[gd].skysig],/even)
  gd = where(chipstr1.dao_nsources ge 0,ngd)
  if ngd gt 0 then expstr[i].dao_nsources = total(chipstr1[gd].dao_nsources)
  gd = where(finite(chipstr1.dao_depth) eq 1,ngd)
  if ngd gt 0 then expstr[i].dao_depth = median([chipstr1[gd].dao_depth],/even)
  gd = where(finite(chipstr1.dao_psfchi) eq 1,ngd)
  if ngd gt 0 then expstr[i].dao_psfchi = median([chipstr1[gd].dao_psfchi],/even)
  gd = where(chipstr1.alf_nsources ge 0,ngd)
  if ngd gt 0 then expstr[i].alf_nsources = total(chipstr[gd].alf_nsources)
  gd = where(finite(chipstr1.alf_depth) eq 1,ngd)
  if ngd gt 0 then expstr[i].alf_depth = median([chipstr1[gd].airmass],/even)
  gd = where(finite(chipstr1.apcor) eq 1,ngd)
  if ngd gt 0 then expstr[i].apcor = median([chipstr1[gd].apcor],/even)
  gd = where(finite(chipstr1.ebv) eq 1,ngd)
  if ngd gt 0 then expstr[i].ebv = median([chipstr1[gd].ebv],/even)
  gd = where(chipstr1.calib_magname ne 'NAN',ngd)
  if ngd gt 0 then expstr[i].magname=chipstr1[gd[0]].calib_magname
  
  printlog,logfile,i+1,iexpnum,expstr[i].filter,expstr[i].exptime,'  ',expstr[i].utdate+'  '+expstr[i].uttime,expstr[i].wcsrms,expstr[i].fwhm,expstr[i].skymode,expstr[i].skysig,$
           expstr[i].dao_nsources,expstr[i].dao_depth,expstr[i].dao_psfchi,expstr[i].alf_nsources,expstr[i].alf_depth,expstr[i].apcor,expstr[i].ebv,$
           format='(I5,A12,A5,F8.2,A2,A24,2F8.4,F8.2,F8.3,I9,F9.2,F8.3,I9,F9.2,F8.3,F8.3)'
Endfor


; Write the summary file
;-----------------------
undefine,head
MKHDR,head,0
sxaddhist,systime(0),head
info = GET_LOGIN_INFO()
sxaddhist,info.user_name+' on '+info.machine_name,head
sxaddhist,'IDL '+!version.release+' '+!version.os+' '+!version.arch,head
sxaddhist,' ',head
sxaddhist,'PHOTRED field summary file for '+field+' - '+finalname,head
sxaddhist,'HDU1 constains a FITS binary table with exposure-level information',head
sxaddhist,'HDU2 constains a FITS binary table with chip-level information',head
sxaddhist,'',head
sxaddhist,'HDU1 columns:',head
sxaddhist,'-----------------------',head
sxaddhist,'Expnum: Exposure number, e.g. 00277653',head
sxaddhist,'Nchips: Number of chips with data, e.g. 59',head
sxaddhist,'Filter: The short filter name, e.g. g',head
sxaddhist,'Exptime: The exposure time in seconds, e.g. 30.0',head
sxaddhist,'UTDate: The UT date, e.g. 2014-01-25',head
sxaddhist,'UTTime: The UT time, e.g. 04:59:07.689347',head
sxaddhist,'Airmass: The median airmass, e.g. 1.12959',head
sxaddhist,'WCSType: The World Coordinate System projection, e.g. TPV',head
sxaddhist,'RA: The Right Ascension at the field center in degrees, e.g. 108.10813',head
sxaddhist,'DEC: The Declination at the field center in degrees, e.g. -53.717605 ',head
sxaddhist,'WCSRMS: Median RMS scatter of the WCS solution in arcsec., e.g. 0.216',head
sxaddhist,'FWHM: Median seeing or average width of sources in pixels, e.g. 4.17',head
sxaddhist,'Skymode: Median sky background in ADU, e.g. 48.6',head
sxaddhist,'Skysig: Median sky scatter in ADU, e.g. 4.11',head
sxaddhist,'DAO_Nsources: Total sources in DAOPHOT ALLSTAR file, e.g. 149265',head
sxaddhist,'DAO_Depth: Median instrumental depth estimate of ALLSTAR file, e.g. 20.3',head
sxaddhist,'DAO_PSFchi: Median chi-squared value of PSF fit by DAOPHOT, e.g. 0.0151',head
sxaddhist,'ALF_Nsources: Total sources in DAOPHOT ALLFRAME file, e.g. 14500',head
sxaddhist,'ALF_Depth: Median instrumental depth of ALLFRAME file, e.g. 20.8',head
sxaddhist,'Apcor: Median additive aperture correction, e.g. 0.0141987',head
sxaddhist,'EBV: Median extinction in B-V, E(B-V), e.g. 0.124',head
sxaddhist,'Magname: Column name in final phot file, e.g. GMAG1',head
sxaddhist,'',head
sxaddhist,'HDU2 columns:',head
sxaddhist,'-----------------------',head
sxaddhist,'Field: The short field name, e.g. F1',head
sxaddhist,'File: The FITS filename, e.g. F1-00277653_01.fits',head
sxaddhist,'Expnum: Exposure number, e.g. 00277653',head
sxaddhist,'Chip: The chip number, e.g. 1',head
sxaddhist,'Base: The base filename, e.g. F1-00277653_01',head
sxaddhist,'Filter: The short filter name, e.g. g',head
sxaddhist,'Exptime: The exposure time in seconds, e.g. 30.0',head
sxaddhist,'UTDate: The UT date, e.g. 2014-01-25',head
sxaddhist,'UTTime: The UT time, e.g. 04:59:07.689347',head
sxaddhist,'Airmass: The airmass, e.g. 1.12959',head
sxaddhist,'Gain: The gain in electrons/ADU, e.g. 4.33970',head
sxaddhist,'Rdnoise: The readout noise in electrons, e.g. 6.60126',head
sxaddhist,'Nx: Number of pixels in the X-dimension, e.g. 2046',head
sxaddhist,'Ny: Number of pixels in the Y-dimension, e.g. 4094',head
sxaddhist,'WCSType: The World Coordinate System projection, e.g. TPV',head
sxaddhist,'Pixscale: The size of one pixel in arcseconds, e.g. 0.262595',head
sxaddhist,'RA: The Right Ascension at the chip center in degrees, e.g. 108.10813',head
sxaddhist,'DEC: The Declination at the chip center in degrees, e.g. -53.717605 ',head
sxaddhist,'WCSRMS: The RMS scatter of the WCS solution in arcseconds, e.g. 0.216',head
sxaddhist,'FWHM: The seeing or average width of sources in pixels, e.g. 4.17',head
sxaddhist,'Skymode: The mode of the sky background in ADU, e.g. 48.6',head
sxaddhist,'Skysig: The scatter of the sky background in ADU, e.g. 4.11',head
sxaddhist,'DAO_Nsources: Number of sources in DAOPHOT ALLSTAR file, e.g. 2549',head
sxaddhist,'DAO_Depth: Instrumental depth estimate of ALLSTAR file, e.g. 20.3',head
sxaddhist,'DAO_NPSFstars: Number of sources that defined the DAOPHOT PSF, e.g. 97',head
sxaddhist,'DAO_PSFType: The type of analytical PSF used by DAOPHOT, e.g. PENNY1',head
sxaddhist,'DAO_PSFboxsize: Size of square PSF box in pixels, e.g. 69',head
sxaddhist,'DAO_PSFvarorder: Order of spatial PSF variations, e.g. 2',head
sxaddhist,'DAO_PSFchi: Chi-squared value of PSF fit by DAOPHOT, e.g. 0.0151',head
sxaddhist,'ALF_Nsources: Number of sources in DAOPHOT ALLFRAME file, e.g. 2500',head
sxaddhist,'ALF_Depth: Instrumental depth estimate of ALLFRAME file, e.g. 20.8',head
sxaddhist,'Calib_Depth: Calibrated depth estimate, e.g. 20.8',head
sxaddhist,'Calib_Color: Name of the color used in the calibration, e.g. u-g',head
sxaddhist,'Calib_ZPterm: Calibration zero-point term, e.g. -0.1087',head
sxaddhist,'Calib_Amterm: Calibration airmass/extinction term, e.g. 0.1012',head
sxaddhist,'Calib_Colorterm: Calibration color term, e.g. -0.0849',head
sxaddhist,'Calib_Magname: Column name in phot file, e.g. GMAG1',head
sxaddhist,'Apcor: The additive aperture correction, e.g. 0.0141987',head
sxaddhist,'EBV: Extinction in B-V, E(B-V), e.g. 0.124',head

printlog,logfile,''
printlog,logfile,'Writing summary file to ',outfile
MWRFITS,0,outfile,head,/create
; Write exposure level FITS binary table in HDU1
MWRFITS,expstr,outfile,/silent
; Write chip-level FITS binary table in HDU2
MWRFITS,chipstr,outfile,/silent

; Copy file to main directory
FILE_COPY,outfile,setupdir,/over,/verbose,/allow

;stop

end
