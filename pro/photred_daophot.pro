;+
;
; PHOTRED_DAOPHOT
;
; This runs daophot for photred
; Can be on pleione or just a single machine depending
; on where this is launched from.
; See daophot_setup.pro for how to make the pleione scripts
;
; INPUTS:
;  /redo Redo files that were already done.
;  /fake Run with artificial star tests.
;  /stp  Stop at the end of the program.
;
; OUTPUTS:
;  The DAOPHOT/ALLSTAR output files.
;
; By D.Nidever  Feb 2008
;-

pro photred_daophot,redo=redo,fake=fake,nmulti=nmulti,stp=stp

COMMON photred,setup

; Don't use images that end in 'a.fits', 's.fits', or 'j.fits'

print,''
print,'########################'
print,'RUNNING PHOTRED_DAOPHOT'
print,'########################'
print,''

CD,current=curdir

; Does the logs/ directory exist?
testlogs = file_test('logs',/directory)
if testlogs eq 0 then FILE_MKDIR,'logs'

; Log files
;----------
thisprog = 'DAOPHOT'
logfile = 'logs/'+thisprog+'.log'
logfile = FILE_EXPAND_PATH(logfile)  ; want absolute filename
if file_test(logfile) eq 0 then SPAWN,'touch '+logfile,out


; Print date and time to the logfile
printlog,logfile,''
printlog,logfile,'Starting PHOTRED_'+thisprog+'  ',systime(0)


; Check that all of the required programs are available
progs = ['readline','readlist','readpar','photred_getinput','photred_updatelists','photred_loadsetup',$
         'photred_mkopt','imfwhm','resistant_mean','undefine','meanclip','pbs_daemon','printlog',$
         'pbs_checkstat','pbs_makescript','mktemp','photred_getgain','photred_getrdnoise','push',$
         'strsplitter','loadinput','touchzero','writeline','first_el','mad','maketemp','rndint',$
         'mpfit','mpfitfun','mpfit2dfun','mpfit2dpeak']
test = PROG_TEST(progs)
if min(test) eq 0 then begin
  bd = where(test eq 0,nbd)
  printlog,logfile,'SOME NECESSARY PROGRAMS ARE MISSING:'
  printlog,logfile,progs[bd]
  return
endif


; Check that the DAOPHOT/ALLSTAR programs exist
SPAWN,'which daophot',out,errout
daophotfile = FILE_SEARCH(out,count=ndaophotfile)
if (ndaophotfile eq 0) then begin
  print,'DAOPHOT PROGRAM NOT AVAILABLE'
  return
endif
SPAWN,'which allstar',out,errout
allstarfile = FILE_SEARCH(out,count=nallstarfile)
if (nallstarfile eq 0) then begin
  print,'ALLSTAR PROGRAM NOT AVAILABLE'
  return
endif


; LOAD THE SETUP FILE if not passed
;-----------------------------------
; This is a 2xN array.  First colume are the keywords
; and the second column are the values.
; Use READPAR.PRO to read it
if n_elements(setup) eq 0 then begin
  PHOTRED_LOADSETUP,setup,count=count
  if count lt 1 then return
endif


; Are we redoing?
if n_elements(redo) eq 0 then begin
  doredo = READPAR(setup,'REDO',count=ndoredo)
  if ndoredo gt 0 and doredo ne '0' and doredo ne '' and doredo ne '-1' then redo=1
  if ndoredo gt 0 and strtrim(doredo,2) eq '0' then redo=0
  if ndoredo eq 0 then redo=0   ; don't redo by default
endif

; Telescope, Instrument
telescope = READPAR(setup,'TELESCOPE',count=ntelescope)
instrument = READPAR(setup,'INSTRUMENT',count=ninstrument)

; Are we running PSFSTARS?
psfcomsrc = READPAR(setup,'PSFCOMSRC',count=npsfcomsrc)
if npsfcomsrc gt 0 and psfcomsrc ne '0' and psfcomsrc ne '' and psfcomsrc ne '-1' then psfcomsrc=1
if npsfcomsrc gt 0 and strtrim(psfcomsrc,2) eq '0' then psfcomsrc=0
if npsfcomsrc eq 0 then psfcomsrc=1  ; use by default

if psfcomsrc ne '0' and psfcomsrc ne '' and psfcomsrc ne '-1' then psfcomsrc=1
if strtrim(psfcomsrc,2) eq '0' then psfcomsrc=0
if n_elements(psfcomsrc) eq 0 then psfcomsrc=1

; Global Common PSF stars?
psfcomglobal = READPAR(setup,'PSFCOMGLOBAL',count=npsfcomglobal)
if npsfcomglobal gt 0 and psfcomglobal ne '0' and psfcomglobal ne '' and psfcomglobal ne '-1' then psfcomglobal=1
if npsfcomglobal gt 0 and strtrim(psfcomglobal,2) eq '0' then psfcomglobal=0
if npsfcomglobal eq 0 then psfcomglobal=1  ; use it by default

; Hyperthread?
hyperthread = READPAR(setup,'hyperthread',count=nhyperthread)
if nhyperthread gt 0 and hyperthread ne '0' and hyperthread ne '' and hyperthread ne '-1' then hyperthread=1
if nhyperthread gt 0 and strtrim(hyperthread,2) eq '0' then hyperthread=0
if nhyperthread eq 0 then hyperthread=0  ; don't use by default

; DAOPHOT .opt values
daopsfva = READPAR(setup,'DAOPSFVA',count=ndaopsfva)
if ndaopsfva eq 0 then undefine,daopsfva
daofitradfwhm = READPAR(setup,'DAOFITRADFWHM',count=ndaofitradfwhm)
if ndaofitradfwhm eq 0 then undefine,daofitradfwhm

; Temporary working directory
workdir = READPAR(setup,'WORKDIR',count=nworkdir)
if nworkdir eq 0 then undefine,workdir

; Clean intermediate files at the end
clean = READPAR(setup,'CLEAN',count=nclean)
if nclean eq 0 then undefine,clean


;###################
; GETTING INPUTLIST
;###################
; INLIST         FITS or FITS.FZ files (must be split)
; OUTLIST        ALS files
; SUCCESSLIST    FITS files

; Get input
;-----------
precursor = 'WCS'
lists = PHOTRED_GETINPUT(thisprog,precursor,redo=redo,ext=['fits','fits.fz'])
ninputlines = lists.ninputlines


; No files to process
;---------------------
if ninputlines eq 0 then begin
  printlog,logfile,'NO FILES TO PROCESS'
  return
endif

inputlines = lists.inputlines


fitsdirlist = FILE_DIRNAME(inputlines)
fitsbaselist = FILE_BASENAME(inputlines)
nfitsbaselist = n_elements(fitsbaselist)

; Unique directories
uidirs = uniq(fitsdirlist,sort(fitsdirlist))
uidirs = uidirs[sort(uidirs)]
dirs = fitsdirlist[uidirs]
ndirs = n_elements(dirs)

; Get the scripts directory from setup
scriptsdir = READPAR(setup,'SCRIPTSDIR')
if scriptsdir eq '' then begin
  printlog,logfile,'NO SCRIPTS DIRECTORY'
  return
endif


; Copy the scripts to the directories
;-------------------------------------
;; Checking that the scripts exist
scripts = ['photo.opt','apcor.opt','daophot.sh','daophot_fake.sh','lstfilter.f','goodpsf.pro','srcfilter.pro']
nscripts = n_elements(scripts)
for i=0,nscripts-1 do begin
  scriptfile = FILE_SEARCH(scriptsdir+'/'+scripts[i],count=nscriptfile)
  if (nscriptfile eq 0) then begin
    printlog,logfile,scriptsdir+'/'+scripts[i],' NOT FOUND'
    return
  endif
endfor
;; Copy the scripts to the directories
for i=0,ndirs-1 do FILE_COPY,scriptsdir+'/'+scripts,dirs[i],/overwrite

;; Compile lstfilter.f it it wasn't already compiled
;;   compiling it locally allows for different architectures/machines
;;   using the same repository
if file_test(curdir+'/lstfilter') eq 0 and file_test(scriptsdir+'/lstfilter.f') then begin
  printlog,logfile,'Compiling lstfilter.f'
  file_copy,scriptsdir+'/lstfilter.f',curdir
  ;; Check which fortran compiler we have
  compiler = ''
  spawn,['which','gfortran'],out,errout,/noshell
  if file_test(strtrim(out[0],2)) eq 1 and errout[0] eq '' then compiler='gfortran'
  spawn,['which','g77'],out,errout,/noshell
  if file_test(strtrim(out[0],2)) eq 1 and errout[0] eq '' then compiler='g77'
  if compiler eq '' then begin
    printlog,logfile,'NO fortran compiler found'
    return
  endif
  ;; Compile
  spawn,[compiler,'lstfilter.f','-o','lstfilter'],out,errout,/noshell
  if file_test('lstfilter') eq 0 or errout[0] ne '' then begin
    printlog,logfile,'ERROR in compiling lstfilter.f'
    return
  endif
  ;; Copy the binary to all the directories
  for i=0,ndirs-1 do file_copy,'lstfilter',dirs[i],/overwrite
endif




; Getting NMULTI from setup file if not given on command line
if n_elements(nmulti) eq 0 then begin
  nmulti = READPAR(setup,'NMULTI')
  nmulti = long(nmulti)

  ; Use NMULTI_DAOPHOT if set
  nmultidaophot = READPAR(setup,'NMULTI_DAOPHOT')
  if nmultidaophot ne '0' and nmultidaophot ne '' and nmultidaophot ne '-1' then nmulti=long(nmultidaophot)
  nmulti = nmulti > 1  ; must be >=1
endif


; LOAD THE "imagers" FILE
;----------------------------
printlog,logfile,'Loading imager information'
imagerstest = FILE_TEST(scriptsdir+'/imagers')
if (imagerstest eq 0) then begin
  printlog,logfile,'NO >>imagers<< file in '+scriptsdir+'  PLEASE CREATE ONE!'
  return
endif
; The columns need to be: Telescope, Instrument, Naps, separator
imagers_fieldnames = ['telescope','instrument','observatory','namps','separator']
imagers_fieldtpes = [7,7,7,3,7]
imagers = IMPORTASCII(scriptsdir+'/imagers',fieldnames=imagers_fieldnames,$
                      fieldtypes=imagers_fieldtypes,comment='#')
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
ind_imager = where(imagers.telescope eq strupcase(telescope) and imagers.instrument eq strupcase(instrument),nind_imager)
if nind_imager eq 0 then begin
  printlog,logfile,'TELESCOPE='+telescope+' INSTRUMENT='+instrument+' NOT FOUND in >>imagers<< file'
  return
endif
thisimager = imagers[ind_imager[0]]
; print out imager info
printlog,logfile,''
printlog,logfile,'USING IMAGER:'
printlog,logfile,'Telescope = '+thisimager.telescope
printlog,logfile,'Instrument = '+thisimager.instrument
printlog,logfile,'Namps = '+strtrim(thisimager.namps,2)
printlog,logfile,"Separator = '"+thisimager.separator+"'"
printlog,logfile,''


;##################################################
;#  PROCESSING THE FILES
;##################################################
; Make the DAPHOT/ALLSTAR option files (.opt and als.opt)
printlog,logfile,''
printlog,logfile,'-----------------------'
printlog,logfile,'PROCESSING THE FILES'
printlog,logfile,'-----------------------'
printlog,logfile,''
printlog,logfile,systime(0)

successarr = intarr(nfitsbaselist)-1                ; 0-bad, 1-good



; Check that the FITS header information can be properly interpreted
;-------------------------------------------------------------------
headerproblem = 0
for i=0,ninputlines-1 do begin

  file = inputlines[i]
  if strmid(file,6,7,/reverse_offset) eq 'fits.fz' then begin
    base = FILE_BASENAME(file,'.fits.fz')
    head = PHOTRED_READFILE(file,exten=1,/header)
  endif else begin
    base = FILE_BASENAME(file,'.fits')
    head = PHOTRED_READFILE(file,/header)
  endelse
  com=''

  ; Checking GAIN
  gain = PHOTRED_GETGAIN(file)
  if gain le 0.0 then com=com+' GAIN ERROR'

  ; Checking READNOISE
  rdnoise = PHOTRED_GETRDNOISE(file)
  if rdnoise le 0.0 then com=com+' READNOISE ERROR'

  ; There were header problems
  if com ne '' then begin
    printlog,logfile,base,com
    PUSH,failurelist,file
    headerproblem = 1
    testing = 1
  endif

endfor

; UPDATE the Lists
PHOTRED_UPDATELISTS,lists,outlist=outlist,successlist=successlist,$
                    failurelist=failurelist,setupdir=curdir,/silent

; There were HEADER problems
if (headerproblem eq 1) then begin
  printlog,logfile,''
  printlog,logfile,'HEADER problems.  RETURNING'
  printlog,logfile,''
  RETALL
endif

; FAKE, artificial star tests
if keyword_set(fake) then begin
  procbaselist = strarr(nfitsbaselist)
  for i=0,nfitsbaselist-1 do begin
    if strmid(fitsbaselist[i],6,7,/reverse_offset) eq 'fits.fz' then $
      procbaselist[i]=FILE_BASENAME(fitsbaselist[i],'.fits.fz') else $
      procbaselist[i]=FILE_BASENAME(fitsbaselist[i],'.fits')
  endfor
  procdirlist = fitsdirlist
  successarr[*] = 1
  goto,rundaophot
endif


printlog,logfile,'Checking which files to make DAOPHOT/ALLSTAR option files for'

; Figure out which files to make OPT files for
undefine,tomakeoptlist
; Loop through directories
for i=0,ndirs-1 do begin

  ; CD to the directory
  CD,dirs[i]

  ; Getting files that are in this directory
  gd = where(fitsdirlist eq dirs[i],ngd)

  ; Some files in this directory
  If (ngd gt 0) then begin

    printlog,logfile,''
    printlog,logfile,'Checking ',dirs[i]
    ;printlog,logfile,'Making OPTION files for ',dirs[i]
    printlog,logfile,''

    ; Looping through files in this directory
    for j=0,ngd-1 do begin

      ind = gd[j]
      fil = fitsbaselist[ind]
      if strmid(fil,6,7,/reverse_offset) eq 'fits.fz' then begin
        fpack = 1
        base = FILE_BASENAME(fil,'.fits.fz')
      endif else begin
        fpack = 0
        base = FILE_BASENAME(fil,'.fits')
      endelse
      printlog,logfile,fil

      ; Make sure that the FITS files are FLOAT
      ;----------------------------------------
      ; Make sure that |BITPIX| > 16
      head = PHOTRED_READFILE(fil,/header)
      bitpix = long(SXPAR(head,'BITPIX',/silent))
      if (bitpix eq 8 or bitpix eq 16) and (fpack eq 0) then begin
        printlog,logfile,'BIXPIX = ',strtrim(bitpix,2),'.  Making image FLOAT'

        ; Read in the image
        im = PHOTRED_READFILE(fil,head,error=error)

        ; Make sure BZERO=0
        bzero = sxpar(head,'BZERO',count=nbzero,/silent)
        if nbzero gt 0 then sxaddpar,head,'BZERO',0.0

        ; Write the FLOAT image
        if n_elements(error) eq 0 and size(im,/type) lt 4 then $
          FITS_WRITE_RESOURCE,fil,float(im),head

        ; There was a problem reading the image
        if n_elements(error) gt 0 then begin
          printlog,logfile,'PROBLEM READING IN ',fil
          successarr[ind] = 0
        endif
      endif

      ; Do the OPT files exist already
      testopt = FILE_TEST(base+'.opt')
      testals = FILE_TEST(base+'.als.opt')


      ; SHOULD we specify the SATURATION limits for the various
      ; CHIPS????


      ; Add to list of files to make OPT file for
      ;--------------------------
      if testopt eq 0 or testals eq 0 or keyword_set(redo) then PUSH,tomakeoptlist,dirs[i]+'/'+fil

    Endfor ; files loop

  endif ; some good files

  ; cd back to original directory
  CD,curdir

endfor ; directories loop


printlog,logfile,'Making DAOPHOT/ALLSTAR option files'
printlog,logfile,systime(0)

ntomakeoptlist = n_elements(tomakeoptlist)
printlog,logfile,strtrim(ntomakeoptlist,2),' files need OPT files'

; Make the OPT files
;-------------------
if ntomakeoptlist gt 0 then begin
  tomakeoptlist_base = file_basename(tomakeoptlist)
  tomakeoptlist_dir = file_dirname(tomakeoptlist)

  ; We could run them in groups of 5.  Lose less "runtime" between checks

  ; Make commands for daophot
  cmd = "PHOTRED_MKOPT,'"+tomakeoptlist_base+"'"
  if n_elements(daopsfva) gt 0 then cmd+=',va='+strtrim(daopsfva,2)
  if n_elements(daofitradfwhm) gt 0 then cmd+=',fitradius_fwhm='+strtrim(daofitradfwhm,2)
  ;cmd = "cd,'"+tomakeoptlist_dir+"' & PHOTRED_MKOPT,'"+tomakeoptlist_base+"'"
  ; Submit the jobs to the daemon
  PBS_DAEMON,cmd,tomakeoptlist_dir,nmulti=nmulti,prefix='dopt',hyperthread=hyperthread,$
             /idle,waittime=1,/cdtodir,scriptsdir=scriptsdir
endif

; Check all OPT files
;---------------------
for i=0,nfitsbaselist-1 do begin

  CD,fitsdirlist[i]

  if strmid(fitsbaselist[i],6,7,/reverse_offset) eq 'fits.fz' then $
    base = file_basename(fitsbaselist[i],'.fits.fz') else $
    base = file_basename(fitsbaselist[i],'.fits')

  ; Check OPT file
  ;----------------
  optfile = base+'.opt'
  testopt = FILE_TEST(optfile)
  if testopt eq 1 then begin
    READCOL,optfile,opttags,dum,optvals,format='A,A,F',/silent

    ; Checking rdnoise, gain and fwhm
    rdnoise = optvals[0]           ; RDNOISE, 1st line
    gain = optvals[1]              ; GAIN, 2nd line
    fwhm = optvals[4]              ; FWHM, 5th line
    ps = optvals[12]               ; PSF radius, 13th line

    ; Checking RDNOISE, GAIN, FWHM
    ; FWHM<=20 for daophot
    if (rdnoise gt 50.) then successarr[i]=0
    if (gain gt 50.) then successarr[i]=0
    if (fwhm gt 20.) then successarr[i]=0
    if (ps gt 51.) then successarr[i]=0

    if (rdnoise gt 50.) then printlog,logfile,optfile,' RDNOISE BAD.  TOO LARGE.'
    if (gain gt 50.) then printlog,logfile,optfile,' GAIN BAD.  TOO LARGE.'
    if (fwhm gt 20.) then printlog,logfile,optfile,' FWHM BAD.  TOO LARGE.'
    if (ps gt 51.) then printlog,logfile,optfile,' PS BAD.  TOO LARGE.'

  ; No OPT file
  endif else begin
    successarr[i]=0
    printlog,logfile,optfile,' NOT FOUND'
  endelse

  ; Check ALS.OPT file
  ;--------------------
  alsfile = base+'.als.opt'
  testals = FILE_TEST(alsfile)
  if testals eq 1 then begin
    READCOL,alsfile,alstags,dum,alsvals,format='A,A,F',/silent

    ; Checking FI, IS, OS
    fi = alsvals[0]               ; fitting radius, 1st line
    is = alsvals[1]               ; inner sky radius, 2nd line
    os = alsvals[2]               ; outer sky radius, 3rd line
    ma = alsvals[9]               ; maximum group size, 10th line

    ; Checking FI, IS, OS
    ; FI<=51 for allstar
    ; IS<=35 for allstar
    ; OS<=100 for allstar
    ; MA<=100 for allstar
    if (fi gt 51.) then successarr[i]=0
    if (is gt 35.) then successarr[i]=0
    if (os gt 100.) then successarr[i]=0
    if (ma gt 100.) then successarr[i]=0

    if (fi gt 51.) then printlog,logfile,optfile,' FI BAD.  TOO LARGE.'
    if (is gt 35.) then printlog,logfile,optfile,' IS BAD.  TOO LARGE.'
    if (os gt 100.) then printlog,logfile,optfile,' OS BAD.  TOO LARGE.'
    if (ma gt 100.) then printlog,logfile,optfile,' MA BAD.  TOO LARGE.'

  ; No ALS.OPT file
  endif else begin
    successarr[i]=0
    printlog,logfile,alsfile,' NOT FOUND'
  endelse

  ; Is everything okay?
  ; If not "bad" then it is "good"
  if successarr[i] ne 0 then successarr[i]=1

  ; cd back to original directory
  CD,curdir

Endfor ; files loop

; Make the Confirmed celestial sources lists, GLOBAL METHOD
;-----------------------------------------------------------
if (psfcomsrc eq 1) and keyword_set(psfcomglobal) then begin

  printlog,logfile,''
  printlog,logfile,'Making CONFIRMED CELESTIAL SOURCES lists - GLOBAL METHOD'
  printlog,logfile,''
  printlog,logfile,systime(0)
  
  ; Get field and chip information for each file
  fitsfieldlist = strarr(nfitsbaselist)
  fitschiplist = strarr(nfitsbaselist)
  for i=0,nfitsbaselist-1 do begin
    if strmid(fitsbaselist[i],6,7,/reverse_offset) eq 'fits.fz' then $
      base = file_basename(fitsbaselist[i],'.fits.fz') else $
      base = file_basename(fitsbaselist[i],'.fits')
    fitsfieldlist[i] = first_el(strsplit(base,'-',/extract))
    fitschiplist[i] = first_el(strsplit(base,thisimager.separator,/extract),/last)
  endfor

  ; Group them into dir/field
  alldirfield = fitsdirlist+' '+fitsfieldlist
  ui = uniq(alldirfield,sort(alldirfield))
  dirfield = alldirfield[ui]
  ndirfield = n_elements(dirfield)
  
  ; Field loop
  for i=0,ndirfield-1 do begin
    idir = first_el(strsplit(dirfield[i],' ',/extract))
    ifield = first_el(strsplit(dirfield[i],' ',/extract),/last)

    CD,idir  ; cd to the appropriate directory
    PHOTRED_COMMONSOURCES_GLOBAL,ifield,setupdir=curdir,redo=redo
    CD,curdir
  endfor
  
endif

; Make the Confirmed celestial sources lists, SINGLE-CHIP, OLD METHOD
;--------------------------------------------------------------------
if (psfcomsrc eq 1) and not keyword_set(psfcomglobal) then begin

  ; Make the confirmed CELESTIAL SOURCES list to be used to make PSF stars
  ;---------------------------------------------------------------------
  printlog,logfile,''
  printlog,logfile,'Making CONFIRMED CELESTIAL SOURCES lists - SINGLE-CHIP METHOD'
  printlog,logfile,''
  printlog,logfile,systime(0)
  
  ; Get field and chip information for each file
  fitsfieldlist = strarr(nfitsbaselist)
  fitschiplist = strarr(nfitsbaselist)
  for i=0,nfitsbaselist-1 do begin
    if strmid(fitsbaselist[i],6,7,/reverse_offset) eq 'fits.fz' then $
      base = file_basename(fitsbaselist[i],'.fits.fz') else $
      base = file_basename(fitsbaselist[i],'.fits')
    fitsfieldlist[i] = first_el(strsplit(base,'-',/extract))
    fitschiplist[i] = first_el(strsplit(base,thisimager.separator,/extract),/last)
  endfor

  ; Group them into dir/field/chip groups
  ;  the different groups shouldn't conflict with each other
  alldirfieldchip = fitsdirlist+' '+fitsfieldlist+' '+fitschiplist
  ui = uniq(alldirfieldchip,sort(alldirfieldchip))
  dirfieldchip = alldirfieldchip[ui]
  ndirfieldchip = n_elements(dirfieldchip)

  ; Make the command files
  cmd = strarr(ndirfieldchip)
  cmnprocdirs = strarr(ndirfieldchip)
  for i=0,ndirfieldchip-1 do begin

    idirfieldchip = dirfieldchip[i]
    ind = where(alldirfieldchip eq idirfieldchip,nind)

    ; Make the CONFIRMED CELESTIAL SOURCES list to be used to make PSF stars
    icmd = "PHOTRED_COMMONSOURCES,['"+strjoin(fitsbaselist[ind],"','")+"'],setupdir='"+curdir+"'"
    if keyword_set(redo) then icmd=icmd+',/redo'
    ;  PHOTRED_COMMONSOURCES,fil,error=psferror,redo=redo
    cmd[i] = icmd
    cmnprocdirs[i] = fitsdirlist[ind[0]]
  endfor

  ; Submit the jobs to the daemon
  PBS_DAEMON,cmd,cmnprocdirs,nmulti=nmulti,prefix='dcmn',hyperthread=hyperthread,$
             /idle,waittime=1,/cdtodir,scriptsdir=scriptsdir
endif


;##############################
;#  RUNNING DAOPHOT
;##############################

; Files with okay option files
gd = where(successarr eq 1,ngd)
if ngd eq 0 then begin
  printlog,logfile,'NO GOOD .OPT and .ALS.OPT FILES FOUND'

  PUSH,failurelist,inputlines

  ; UPDATE the Lists
  PHOTRED_UPDATELISTS,lists,outlist=outlist,successlist=successlist,$
                      failurelist=failurelist,setupdir=curdir,/silent

  return
endif

; These are the files to process
procbaselist = strarr(ngd)
procextlist = strarr(ngd)
for i=0,ngd-1 do begin
  if strmid(fitsbaselist[gd[i]],6,7,/reverse_offset) eq 'fits.fz' then $
    procextlist[i]='.fits.fz' else procextlist[i]='.fits'
  procbaselist[i]=FILE_BASENAME(fitsbaselist[gd[i]],procextlist[i])
endfor
procdirlist = fitsdirlist[gd]
nprocbaselist = n_elements(procbaselist)


; Have some been done already??
;------------------------------
If not keyword_set(redo) then begin

  ; Checking possible output files 
  donearr = lonarr(nprocbaselist)
  for i=0,nprocbaselist-1 do begin
    base = procdirlist[i]+'/'+procbaselist[i]
    ; Check if this file has an ALS file
    alstest = FILE_TEST(base+'.als')
    if alstest eq 1 then alslines=FILE_LINES(base+'.als') else alslines=0
    ; Check if this file has an A.ALS file
    aalstest = FILE_TEST(base+'a.als')
    if aalstest eq 1 then aalslines=FILE_LINES(base+'a.als') else aalslines=0
    ; Done properly?
    if (alslines gt 3 and aalslines gt 3) then donearr[i]=1
  endfor

  ; Some done already. DO NOT REDO
  bd = where(donearr eq 1,nbd)
  if nbd gt 0 then begin
    ; Print out the names
    printlog,logfile,''
    for i=0,nbd-1 do $
      printlog,logfile,procdirlist[bd[i]]+'/'+procbaselist[bd[i]]+' DAOPHOT ALREADY DONE'

    ; Add these to the "success" and "outlist" list
    PUSH,successlist,procdirlist[bd]+'/'+procbaselist[bd]+procextlist[bd]
    PUSH,outlist,procdirlist[bd]+'/'+procbaselist[bd]+'.als'
    PHOTRED_UPDATELISTS,lists,outlist=outlist,successlist=successlist,$
                        failurelist=failurelist,setupdir=curdir,/silent

    ; Remove them from the arrays
    if nbd lt nprocbaselist then REMOVE,bd,procbaselist,procdirlist
    if nbd eq nprocbaselist then UNDEFINE,procbaselist,procdirlist
    nprocbaselist = n_elements(procbaselist)

    printlog,logfile,''
    printlog,logfile,'REMOVING '+strtrim(nbd,2)+' files from INLIST.  '+$
                     strtrim(nprocbaselist,2)+' files left to PROCESS'
    printlog,logfile,''
  endif

  ; No files to run
  if nprocbaselist eq 0 then begin
    printlog,logfile,'NO FILES TO PROCESS'
    return
  endif

Endif ; some done already?


;----------------------
; RUNNING THE COMMANDS
;----------------------
rundaophot:

printlog,logfile,''
printlog,logfile,'Running DAOPHOT on '+strtrim(n_elements(procbaselist),2)+' files'
printlog,logfile,''
printlog,logfile,systime(0)

; Make commands for daophot
if not keyword_set(fake) then begin
  cmd = './daophot.sh '+procbaselist
endif else begin
  cmd = './daophot_fake.sh '+procbaselist
endelse
if n_elements(workdir) gt 0 then cmd+=' '+workdir  ; temporary working directory
  
; Submit the jobs to the daemon
PBS_DAEMON,cmd,procdirlist,nmulti=nmulti,prefix='dao',hyperthread=hyperthread,$
           waittime=1,/cdtodir,scriptsdir=scriptsdir

; IT WOULD BE BETTER TO UPDATE THE LISTS
; AFTER EACH FILE IS PROCESSED!!!!!
; CAN PBS_DAEMON DO THAT?????

;-------------------
; Checking OUTPUTS
;-------------------

; Loop through all files in fitsbaselist
for i=0,nfitsbaselist-1 do begin

  ; Only check "good" files
  if (successarr[i] eq 1) then begin

    CD,fitsdirlist[i]

    fil = fitsbaselist[i]
    if strmid(fil,6,7,/reverse_offset) eq 'fits.fz' then base=FILE_BASENAME(fil,'.fits.fz') else $
      base = FILE_BASENAME(fil,'.fits')
    ; Check that this file has an ALS file
    alstest = FILE_TEST(base+'.als')
    if alstest eq 1 then alslines=FILE_LINES(base+'.als') else alslines=0
    if (alstest eq 0 or alslines lt 3) then successarr[i]=0
    if (alstest eq 0) then printlog,logfile,base+'.als NOT FOUND'
    ; Check the A.ALS file
    if not keyword_set(fake) then begin
      aalstest = FILE_TEST(base+'a.als')
      if aalstest eq 1 then aalslines=FILE_LINES(base+'a.als') else aalslines=0
      if (aalstest eq 0 or aalslines lt 3) then successarr[i]=0
      if (aalstest eq 0) then printlog,logfile,base+'a.als NOT FOUND'
    endif

  endif

  CD,curdir

endfor


;##########################################
;#  UPDATING LIST FILES
;##########################################
undefine,outlist,successlist,failurelist

; Success List
ind = where(successarr eq 1,nind)
if nind gt 0 then successlist = inputlines[ind] else UNDEFINE,successlist

; Output List
; Creating the new output array, ALS files
if (nind gt 0) then begin
  bases = strarr(nind)
  for i=0,nind-1 do begin
    if strmid(fitsbaselist[ind[i]],6,7,/reverse_offset) eq 'fits.fz' then $
      bases[i]=FILE_BASENAME(fitsbaselist[ind[i]],'.fits.fz') else $
      bases[i]=FILE_BASENAME(fitsbaselist[ind[i]],'.fits')
  endfor
  outlist = fitsdirlist[ind]+'/'+bases+'.als'
endif else UNDEFINE,outlist

; Failure List
bd = where(successarr eq 0,nbd)
if (nbd gt 0) then begin
  failurelist = inputlines[bd]
endif else UNDEFINE,failurelist

PHOTRED_UPDATELISTS,lists,outlist=outlist,successlist=successlist,$
                    failurelist=failurelist,setupdir=curdir

;;######################
;;  CLEANING UP
;;######################
if keyword_set(clean) then begin
  printlog,logfile,'CLEANING UP.  CLEAN='+strtrim(clean,2)

  ;; Only clean up for successful files
  READLIST,curdir+'/logs/DAOPHOT.success',fitsfiles,/unique,/fully,setupdir=curdir,count=nfitsfiles,logfile=logfile,/silent
  for i=0,nfitsfiles-1 do begin
    dir1 = file_dirname(fitsfiles[i])
    base1 = file_basename(fitsfiles[i])
    if strmid(base1,6,7,/reverse_offset) eq 'fits.fz' then base1 = FILE_BASENAME(base1,'.fits.fz') else $
      base1 = FILE_BASENAME(base1,'.fits')
    ;; lst, lst1, lst2, lst1.chi, grp, nst, lst2.chi, plst.chi,
    ;; nei, als.inp, a.fits, cmn.log, cmn.coo, cmn.ap, cmn.lst
    FILE_DELETE,dir1+'/'+base1+['.lst','.lst1','.lst2','.lst1.chi','.lst2.chi','.grp','.nst','.plst.chi',$
                                '.nei','.als.inp','.cmn.log','.cmn.coo','.cmn.ap','.cmn.lst','a.fits','a.fits.fz'],/allow
    ;; If CLEAN=2, then also remove coo, ap, s.fits as well
    if clean ge 2 then FILE_DELETE,dir1+'/'+base1+['.coo','.ap','s.fits','s.fits.fz'],/allow
  endfor
endif

printlog,logfile,'PHOTRED_DAOPHOT Finished  ',systime(0)

if keyword_set(stp) then stop

end
