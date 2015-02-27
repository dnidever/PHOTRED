pro stdred_aperphot,redo=redo,stp=stp

;+
;
; STDRED_APERPHOT
;
; This runs daophot aperture photometry for stdred
;
; INPUTS:
;  /redo Redo files that were already done.
;  /stp  Stop at the end of the program.
;
; OUTPUTS:
;  The DAOPHOT/PHOTOMETRY .tot output files.
;
; By D.Nidever  May 2008
;-

COMMON photred,setup


print,''
print,'########################'
print,'RUNNING STDRED_APERPHOT'
print,'########################'
print,''

CD,current=curdir

; Does the logs/ directory exist?
testlogs = file_test('logs',/directory)
if testlogs eq 0 then FILE_MKDIR,'logs'

; Log files
;----------
thisprog = 'APERPHOT'
logfile = 'logs/'+thisprog+'.log'
logfile = FILE_EXPAND_PATH(logfile)  ; want absolute filename
if file_test(logfile) eq 0 then SPAWN,'touch '+logfile,out


; Print date and time to the logfile
printlog,logfile,''
printlog,logfile,'Starting STDRED_'+thisprog+'  ',systime(0)


; Check that all of the required programs are available
progs = ['readline','readlist','readpar','photred_getinput','photred_updatelists','photred_loadsetup',$
         'photred_mkopt','imfwhm','resistant_mean','undefine','meanclip','mktemp','photred_getgain',$
         'photred_getrdnoise','printlog','push','roi_cut','srcmatch','writeline','array_indices2',$
         'first_el','mad','maxloc','minloc','mpfit','range','scale','scale_vector','slope','stringize',$
         'strsplitter','loadinput','touchzero','signs','strmult','strtrim0','sign']
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
;SPAWN,'which allstar',out,errout
;allstarfile = FILE_SEARCH(out,count=nallstarfile)
;if (nallstarfile eq 0) then begin
;  print,'ALLSTAR PROGRAM NOT AVAILABLE'
;  return
;endif


; LOAD THE SETUP FILE if not passed
;-----------------------------------
; This is a 2xN array.  First colume are the keywords
; and the second column are the values.
; Use READPAR.PRO to read it
if n_elements(setup) eq 0 then begin
  PHOTRED_LOADSETUP,setup,count=count,/std
  if count lt 1 then return
endif


; Are we redoing?
doredo = READPAR(setup,'REDO')
if keyword_set(redo) or (doredo ne '-1' and doredo ne '0') then redo=1

; Telescope, Instrument
telescope = READPAR(setup,'TELESCOPE')
instrument = READPAR(setup,'INSTRUMENT')
observatory = READPAR(setup,'OBSERVATORY')
if strlowcase(telescope) eq 'blanco' then observatory='ctio'
if strlowcase(telescope) eq 'swope' then observatory='lco'
if strlowcase(telescope) eq 'magellan' then observatory='lco'
if strlowcase(telescope) eq 'lbt' then observatory='mgio'
if observatory eq '0' or observatory eq '-1' or observatory eq '' then begin
  printlog,logfile,'OBSERVATORY not input'
  return
endif
; Get observatory structure
OBSERVATORY,observatory,obs_struct
if obs_struct.name eq '' then begin
  printlog,logfile,'Observatory >>',observatory,'<< NOT FOUND'
  return
endif

; Hyperthread?
hyperthread = READPAR(setup,'hyperthread')
if hyperthread ne '0' and hyperthread ne '' and hyperthread ne '-1' then hyperthread=1
if strtrim(hyperthread,2) eq '0' then hyperthread=0
if n_elements(hyperthread) eq 0 then hyperthread=0


;###################
; GETTING INPUTLIST
;###################
; INLIST         FITS files (must be split)
; OUTLIST        ALS files
; SUCCESSLIST    FITS files

; Get input
;-----------
precursor = 'WCS'
lists = PHOTRED_GETINPUT(thisprog,precursor,redo=redo,ext='fits')
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
; Checking that the scripts exist
scripts = ['photo.opt','aperphot.sh']
;scripts = ['photo.opt','apcor.opt','daophot.sh','lstfilter','goodpsf.pro']
nscripts = n_elements(scripts)
for i=0,nscripts-1 do begin
  scriptfile = FILE_SEARCH(scriptsdir+'/'+scripts[i],count=nscriptfile)
  if (nscriptfile eq 0) then begin
    print,scriptdir+'/'+scripts[i],' NOT FOUND'
    return
  endif
end
; Copy the scripts to the directories
for i=0,ndirs-1 do begin
  FILE_COPY,scriptsdir+'/'+scripts,dirs[i],/overwrite
end

; Getting NMULTI
nmulti = READPAR(setup,'NMULTI')
nmulti = long(nmulti)


;##################################################
;#  PROCESSING THE FILES
;##################################################
; Make the DAPHOT/ALLSTAR option files (.opt and als.opt)
printlog,logfile,''
printlog,logfile,'-----------------------'
printlog,logfile,'PROCESSING THE FILES'
printlog,logfile,'-----------------------'
printlog,logfile,''

; Initialzing some arrays
UNDEFINE,outlist,successlist,failurelist

successarr = intarr(nfitsbaselist)-1                ; 0-bad, 1-good



; Check that the FITS header information can be properly interpreted
;-------------------------------------------------------------------
headerproblem = 0
for i=0,ninputlines-1 do begin

  file = inputlines[i]
  base = FILE_BASENAME(file,'.fits')
  head = HEADFITS(file)
  com=''

  ; Checking GAIN
  gain = PHOTRED_GETGAIN(file)
  if gain lt 0.0 then com=com+' GAIN ERROR'

  ; Checking READNOISE
  rdnoise = PHOTRED_GETRDNOISE(file)
  if rdnoise lt 0.0 then com=com+' READNOISE ERROR'

  ; There were header problems
  if com ne '' then begin
    printlog,logfile,base,com
    PUSH,failurelist,file
    headerproblem = 1
    testing = 1
  endif

end

; UPDATE the Lists
PHOTRED_UPDATELISTS,lists,outlist=outlist,successlist=successlist,$
                    failurelist=failurelist,/silent

; There were HEADER problems
if (headerproblem eq 1) then begin
  printlog,logfile,''
  printlog,logfile,'HEADER problems.  RETURNING'
  printlog,logfile,''
  RETALL
endif



printlog,logfile,''
printlog,logfile,'-----------------------------------'
printlog,logfile,'Making DAOPHOT/ALLSTAR option files'
printlog,logfile,'-----------------------------------'
printlog,logfile,''


; Loop through directories
for i=0,ndirs-1 do begin

  ; CD to the directory
  CD,dirs[i]

  ; Getting files that are in this directory
  gd = where(fitsdirlist eq dirs[i],ngd)

  ; Some files in this directory
  If (ngd gt 0) then begin

    printlog,logfile,''
    printlog,logfile,'Making OPTION files for ',dirs[i]
    printlog,logfile,''

    ; Looping through files in this directory
    for j=0,ngd-1 do begin

      ind = gd[j]
      longfile = fitsdirlist[ind]+'/'+fitsbaselist[ind]
      fil = fitsbaselist[ind]
      base = FILE_BASENAME(fil,'.fits')
      printlog,logfile,fil

      ; Make sure that the FITS files are FLOAT
      ;----------------------------------------
      ; Make sure that |BITPIX| > 16
      head = HEADFITS(fil)
      bitpix = long(SXPAR(head,'BITPIX'))
      if (bitpix eq 8 or bitpix eq 16) then begin
        printlog,logfile,'BIXPIX = ',strtrim(bitpix,2),'.  Making image FLOAT'

        ; Read in the image
        FITS_READ,fil,im,head,/no_abort,message=message

        ; Write the FLOAT image
        if (message[0] eq '') then $
        FITS_WRITE,fil,float(im),head

        ; There was a problem reading the image
        if (message[0] ne '') then begin
          printlog,logfile,'PROBLEM READING IN ',fil
          successarr[ind] = 0
        endif
      endif

      ; Do the OPT files exist already
      testopt = FILE_TEST(base+'.opt')
      testals = FILE_TEST(base+'.als.opt')


      ; SHOULD we specify the SATURATION limits for the various
      ; CHIPS????


      ; Make daophot option files
      ;--------------------------
      ;MKOPT,base+'.fits'
      if testopt eq 0 or testals eq 0 then begin
        PHOTRED_MKOPT,base+'.fits',error=error
        if n_elements(error) gt 0 then begin
          printlog,logfile,'PHOTRED_MKOPT Error'
          successarr[ind]=0
        endif
      endif

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
        if (rdnoise gt 50.) then successarr[ind]=0
        if (gain gt 50.) then successarr[ind]=0
        if (fwhm gt 20.) then successarr[ind]=0
        if (ps gt 51.) then successarr[ind]=0

        if (rdnoise gt 50.) then printlog,logfile,optfile,' RDNOISE BAD.  TOO LARGE.'
        if (gain gt 50.) then printlog,logfile,optfile,' GAIN BAD.  TOO LARGE.'
        if (fwhm gt 20.) then printlog,logfile,optfile,' FWHM BAD.  TOO LARGE.'
        if (ps gt 51.) then printlog,logfile,optfile,' PS BAD.  TOO LARGE.'

      ; No OPT file
      endif else begin
        successarr[ind]=0
        printlog,logfile,optfile,' NOT FOUND'
      endelse


      ; Is everything okay?
      ; If not "bad" then it is "good"
      if successarr[ind] ne 0 then successarr[ind]=1

      ; Failure
      if successarr[ind] eq 0 then PUSH,failurelist,longfile


    End ; files loop


  ; No files in this directory
  Endif else begin
    ; Do nothing
  Endelse

  ; cd back to original directory
  CD,curdir



end

; Updating the lists
PHOTRED_UPDATELISTS,lists,outlist=outlist,successlist=successlist,$
                    failurelist=failurelist,/silent




;##############################
;#  RUNNING DAOPHOT
;##############################

; Files with okay option files
gd = where(successarr eq 1,ngd)
if ngd eq 0 then begin
  print,'NO GOOD .OPT and .ALS.OPT FILES FOUND'
  return
endif

printlog,logfile,''
printlog,logfile,'--------------------------'
printlog,logfile,'RUNNING DAOPHOT/PHOTOMETRY'
printlog,logfile,'--------------------------'
printlog,logfile,''

; These are the files to process
procbaselist = FILE_BASENAME(fitsbaselist[gd],'.fits')
procdirlist = fitsdirlist[gd]
nprocbaselist = n_elements(procbaselist)




;----------------------
; RUNNING THE COMMANDS
;----------------------

; Make commands for daophot
cmd = './aperphot.sh '+procbaselist

; Submit the jobs to the daemon
PBS_DAEMON,cmd,procdirlist,nmulti=nmulti,prefix='dao',hyperthread=hyperthread


;-------------------
; Checking OUTPUTS
;-------------------

; Loop through all files in procbaselist
For i=0,nprocbaselist-1 do begin

  ; CD to the appropriate directory
  CD,procdirlist[i]

  longfile = procdirlist[i]+'/'+procbaselist[i]+'.fits'
  base = procbaselist[i]

  ; Check the outputs
  cootest = FILE_TEST(base+'.coo')
  if cootest eq 1 then coolines=FILE_LINES(base+'.coo') else coolines=0
  aptest = FILE_TEST(base+'.ap')
  if aptest eq 1 then aplines=FILE_LINES(base+'.ap') else aplines=0
  acootest = FILE_TEST(base+'a.coo')
  if acootest eq 1 then acoolines=FILE_LINES(base+'a.coo') else acoolines=0
  aaptest = FILE_TEST(base+'a.ap')
  if aaptest eq 1 then aaplines=FILE_LINES(base+'a.ap') else aaplines=0


  ; Successful
  if (coolines ge 4 and aplines ge 4 and acoolines ge 4 and aaplines ge 4) then begin
    printlog,logfile,'Nstars = ',strtrim(coolines-3,2),' Found'
    PUSH,successlist,longfile
    PUSH,outlist,procdirlist[i]+'/'+base+'.ap'

  ; Failure
  endif else begin
    printlog,logfile,'DAOPHOT PROBLEMS FOR ',base
    PUSH,failurelist,longfile
  endelse

  CD,curdir

Endfor

;#####################
; UPDATE the Lists
;#####################
PHOTRED_UPDATELISTS,lists,outlist=outlist,successlist=successlist,$
                    failurelist=failurelist,/silent


printlog,logfile,'STDRED_APERPHOT Finished  ',systime(0)

if keyword_set(stp) then stop

end
