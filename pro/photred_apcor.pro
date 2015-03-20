pro photred_apcor,redo=redo,stp=stp

;+
;
; PHOTRED_APCOR
;
; This does the aperture correction for photred
; See cp_daogrow.pro, mk_daogrow_all.pro, mk_daogrow.pro, and
; daogrow.sh
;
; INPUTS:
;  /redo Redo files that were already done.
;  /stp  Stop at the end of the program.
;
; OUTPUTS:
;  The "apcor.lst" file that contains the aperture corrections
;  for all the ALS files
;
; By D.Nidever  Mar 2008
;-

COMMON photred,setup


print,''
print,'########################'
print,'RUNNING PHOTRED_APCOR'
print,'########################'
print,''

CD,current=curdir

; Does the logs/ directory exist?
testlogs = file_test('logs',/directory)
if testlogs eq 0 then FILE_MKDIR,'logs'

; Log files
;----------
thisprog = 'APCOR'
logfile = 'logs/'+thisprog+'.log'
logfile = FILE_EXPAND_PATH(logfile)  ; want absolute filename
if FILE_TEST(logfile) eq 0 then SPAWN,'touch '+logfile,out


; Print date and time to the logfile
printlog,logfile,''
printlog,logfile,'Starting PHOTRED_'+thisprog+'  ',systime(0)


; Check that all of the required programs are available
progs = ['readline','readlist','readpar','photred_getfilter','photred_getuttime',$
         'photred_getexptime','photred_getinput','photred_updatelists','photred_loadsetup',$
         'undefine','push','printlog','mkdel','apcor','loadinput','resistant_mean',$
         'photred_getairmass','photred_getdate','badpar','airmass','first_el','strsplitter',$
         'touchzero','writeline','mktemp','sexig2ten']
test = PROG_TEST(progs)
if min(test) eq 0 then begin
  bd = where(test eq 0,nbd)
  printlog,logfile,'SOME NECESSARY PROGRAMS ARE MISSING:'
  printlog,logfile,progs[bd]
  return
endif


; Check that the DAOGROW programs exist
SPAWN,'which daogrow',out,errout
daogrowfile = FILE_SEARCH(out,count=ndaogrowfile)
if (ndaogrowfile eq 0) then begin
  print,'DAOGROW PROGRAM NOT AVAILABLE'
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
doredo = READPAR(setup,'REDO')
if keyword_set(redo) or (doredo ne '-1' and doredo ne '0') then redo=1

; Telescope, Instrument
telescope = READPAR(setup,'TELESCOPE')
instrument = READPAR(setup,'INSTRUMENT')
observatory = READPAR(setup,'OBSERVATORY')
if observatory eq '0' or observatory eq '-1' or observatory eq '' then undefine,observatory
if strlowcase(telescope) eq 'blanco' then observatory='ctio'
if strlowcase(telescope) eq 'swope' then observatory='lco'
if strlowcase(telescope) eq 'magellan' then observatory='lco'
if strlowcase(telescope) eq 'lbt' then observatory='mgio'


;###################
; GETTING INPUTLIST
;###################
; INLIST         FITS files
; OUTLIST        FITS files
; SUCCESSLIST    FITS files

; Takes all of the FITS files from DAOPHOT.success(!!) (because the files in DAOPHOT.outlist will already
; have been moved by MATCH)and COPIES them into APCOR.inlist.
; All FITS files that have an aperture correction in the final apcor.lst get put into the apcor.outlist.
; successlist and outlist are the same, although files will accumulate in the successlist over multiple
; runs while the outlist will get remade each time.
;
; Remove stars in success list from inlist unless /redo set

; Get input
;-----------
precursor = 'DAOPHOT.success'
lists = PHOTRED_GETINPUT(thisprog,precursor,redo=redo,ext='fits')
ninputlines = lists.ninputlines


; No files to process
;---------------------
if ninputlines eq 0 then begin
  printlog,logfile,'NO FILES TO PROCESS'
  return
endif

inputlines = lists.inputlines




;########################################
;#  PROCESSING THE FILES
;########################################
printlog,logfile,''
printlog,logfile,'-----------------------'
printlog,logfile,'PROCESSING THE FILES'
printlog,logfile,'-----------------------'
printlog,logfile,''


; Make a daogrow directory
test = file_test('daogrow',/directory)
; directory exists already, erase and start again
if test eq 1 then FILE_DELETE,'daogrow',/recursive
FILE_MKDIR,'daogrow'
daogrowdir = curdir+'/daogrow/'

; Get the scripts directory from setup
scriptsdir = READPAR(setup,'SCRIPTSDIR')
if scriptsdir eq '' then begin
  printlog,logfile,'NO SCRIPTS DIRECTORY'
  return
endif


successarr = intarr(ninputlines)-1          ; 0-bad, 1-good

;-----------------------------------------------------
; Copying the A.ALS and A.AP file to daogrow/
; and Checking that the FITS file exists
;-----------------------------------------------------
printlog,logfile,'Copying ALS and AP files to daogrow/'
for i=0,ninputlines-1 do begin

  longfile = inputlines[i]
  file = FILE_BASENAME(longfile)
  filedir = FILE_DIRNAME(longfile)

  CD,filedir


  ; Check that the FITS file exists
  ;--------------------------------
  fitstest = FILE_TEST(file)
  if (fitstest eq 0) then begin
    successarr[i]=0
    printlog,logfile,file,' NOT FOUND'
    goto,BOMB
  endif

  ; Check that this file has the FITS ending
  ;-----------------------------------------
  len = strlen(file)
  ending = strmid(file,len-4,4)
  if (ending ne 'fits') then begin
    successarr[i] = 0
    printlog,logfile,file,' DOES NOT END IN >>fits<<'
    goto,BOMB
  endif

  ; Check that the A.ALS and A.AP files exist
  ;------------------------------------------
  base = FILE_BASENAME(file,'.fits')
  alsfile = base+'a.als'
  apfile = base+'a.ap'
  alstest = FILE_TEST(alsfile)
  aptest = FILE_TEST(apfile)
  if (alstest eq 0) or (aptest eq 0) then begin
    successarr[i] = 0
    if alstest eq 0 then printlog,logfile,alsfile,' NOT FOUND'
    if aptest eq 0 then printlog,logfile,apfile,' NOT FOUND'
    goto,BOMB
  endif

  ; Check that the A.ALS and A.AP files have SOURCES in them
  ;----------------------------------------------------------
  nals_sources = FILE_LINES(alsfile)-3
  nap_sources = FILE_LINES(apfile)-3
  if (nals_sources lt 1) or (nap_sources lt 1) then begin
    successarr[i] = 0
    if nals_sources lt 1 eq 0 then printlog,logfile,alsfile,' has NO SOURCES'
    if nap_sources lt 1 eq 0 then printlog,logfile,apfile,' has NO SOURCES'
    goto,BOMB
  endif

  ; Check A.ALS header
  ;-----------------
  ;This is what it should look like:
  ; NL    NX    NY  LOWBAD HIGHBAD  THRESH     AP1  PH/ADU  RNOISE    FRAD
  ; 1  2040  2047    79.9 30000.0   32.07    3.00    2.20    2.86    2.18
  ; 
  ;   59  957.205   12.228   13.961   0.0151  141.334       4.    0.547   -0.020
  line1='' & line2='' & line3=''
  openr,unit,/get_lun,alsfile
  readf,unit,line1
  readf,unit,line2
  readf,unit,line3
  close,unit
  free_lun,unit

  dum1 = strtrim(strsplit(line1,' ',/extract),2)
  dum2 = strtrim(strsplit(line2,' ',/extract),2)
  line3 = strtrim(line3,2)
  if (dum1[0] ne 'NL') or (dum2[0] ne '1') or (line3 ne '') then begin
    successarr[i]=0
    printlog,logfile,alsfile,' DOES HOT HAVE A PROPER ALS HEADER'
    goto,BOMB
  endif

  ; Check A.AP header
  ;----------------
  ;This is what it should look like:
  ; NL    NX    NY  LOWBAD HIGHBAD  THRESH     AP1  PH/ADU  RNOISE    FRAD
  ;  2  2040  2047    79.9 30000.0   32.07    3.00    2.20    2.86    2.00
  line1='' & line2='' & line3=''
  openr,unit,/get_lun,apfile
  readf,unit,line1
  readf,unit,line2
  readf,unit,line3
  close,unit
  free_lun,unit

  dum1 = strtrim(strsplit(line1,' ',/extract),2)
  dum2 = strtrim(strsplit(line2,' ',/extract),2)
  line3 = strtrim(line3,2)
  if (dum1[0] ne 'NL') or (dum2[0] ne '2') or (line3 ne '') then begin
    successarr[i]=0
    printlog,logfile,apfile,' DOES HOT HAVE A PROPER AP HEADER'
    goto,BOMB
  endif


  ; Copy the files to daogrow/
  ;----------------------------
  FILE_COPY,alsfile,daogrowdir,/overwrite,/allow_same
  FILE_COPY,apfile,daogrowdir,/overwrite,/allow_same

  ; This file was okay
  successarr[i] = 1

  BOMB:

  CD,curdir

end

ncopy = long(total(successarr))
printlog,logfile,strtrim(ncopy,2),'/',strtrim(ninputlines,2),' successfully copied to daogrow/'


;--------------------------------
; Making the DAOGROW input files
;--------------------------------
printlog,logfile,'Making DAOGROW input files: daogrow/daogrow.inf and daogrow/daogrow.ext'

OPENW,infunit,/get_lun,daogrowdir+'daogrow.inf'
OPENW,extunit,/get_lun,daogrowdir+'daogrow.ext'

; Loop through the FITS files
For i=0,ninputlines-1 do begin

  ; Only use "good" files
  if (successarr[i] eq 1) then begin

    filedir = FILE_DIRNAME(inputlines[i])
    ;fitsfile = FILE_BASENAME(inputlines[i])
    fitsfile = inputlines[i]
    base = FILE_BASENAME(fitsfile,'.fits')

    ; Getting information from header
    head = HEADFITS(fitsfile)

    ; Getting FILTER
    filt = PHOTRED_GETFILTER(fitsfile,/numeric)

    ; Getting UT
    time = PHOTRED_GETUTTIME(fitsfile)
    if (time eq '') then begin
      successarr[i] = 0
      printlog,logfile,'NO UT-TIME'
      goto,BOMB2
    endif
    timarr = strsplit(time,':',/extract)
    uthr = timarr[0]
    utmin = timarr[1]

    ; Getting AIRMASS
    am = PHOTRED_GETAIRMASS(fitsfile,obs=observatory,/update)
    ;am = sxpar(head,'AIRMASS')
    ; no airmass, get from UT and RA/DEC
    if am lt 0.9 then begin
      successarr[i] = 0
      printlog,logfile,'NO AIRMASS'
      goto,BOMB2
    end

    ; Getting EXPTIME
    exp = PHOTRED_GETEXPTIME(fitsfile)
    if strtrim(exp,2) eq '-1' then begin
      successarr[i] = 0
      printlog,logfile,'NO EXPTIME'
      goto,BOMB2
    endif


    ; How long is the name
    len = strlen(base+'a')
    if len gt 23 then begin
      printlog,logfile,base+'a',' IS TOO LONG.  MAX 23 characters'
      return
    endif


    ; Outputting to the INF file
    ;format='(A1,A-11,A22,A4,A3,F7.3,F10.3)'
    format='(A1,A-23,A10,A4,A3,F7.3,F10.3)'
    printf,infunit,format=format,'',base+'a',filt,uthr,utmin,am,exp

    ; Outputting to the EXT file
    printf,extunit,base+'a.ap'

  endif  ; "good" file

  BOMB2:

Endfor

CLOSE,infunit
FREE_LUN,infunit
CLOSE,extunit
FREE_LUN,extunit


;-----------------------
; Run daogrow.sh script
;-----------------------
CD,'daogrow'

FILE_COPY,scriptsdir+'/daogrow.sh','.',/overwrite
FILE_COPY,scriptsdir+'/photo.opt','.',/overwrite

print,'RUNNING DAOGROW'
SPAWN,'./daogrow.sh daogrow > daogrow.log',out,errout


; DAOGROW Converged
;----------------------
; Read in the log file
READLINE,'daogrow.log',lines
ind = first_el(where(stregex(lines,'Converged',/boolean) eq 1,nind),/last)
if nind gt 0 then begin

  printlog,logfile,lines[ind[0]:*]
  printlog,logfile,'DAOGROW CONVERGED'


  ;--------------------
  ; Running MKDEL.PRO
  ;--------------------
  printlog,logfile,''
  printlog,logfile,'RUNNING MKDEL.PRO'
  printlog,logfile,''
  MKDEL,'daogrow.inf'


  ;--------------------
  ; Running APCOR.PRO
  ;--------------------
  printlog,logfile,''
  printlog,logfile,'RUNNING APCOR.PRO'
  printlog,logfile,''
  APCOR,'*.del','apcor.lst'

  ; Back to the original directory
  CD,curdir


  ;---------------------
  ; Checking OUTPUTS
  ;---------------------

  ; Check that apcor.lst file exists
  ;------------------------------------
  test = FILE_TEST('daogrow/apcor.lst')
  if (test eq 1) then begin

    ; Copy final apcor.lst file to the original directory
    FILE_COPY,'daogrow/apcor.lst','.',/overwrite,/allow

    printlog,logfile,''
    printlog,logfile,'FINAL Aperture Correction File = >>apcor.lst<<'


    ; Checking files in apcor.lst
    ;-----------------------------
    READCOL,'daogrow/apcor.lst',apcnames,apcvalue,format='A,F',/silent
    ; Names have a.del endings

    ; Loop through all files in fitsbaselist
    for i=0,ninputlines-1 do begin

      ; Only check "good" files
      if (successarr[i] eq 1) then begin

        fil = FILE_BASENAME(inputlines[i])
        base = FILE_BASENAME(fil,'.fits')
        delfile = base+'a.del'

        ; Check if the del file is in the apcor list
        gd = where(apcnames eq delfile,ngd)

        if (ngd eq 0) then begin
          successarr[i]=0
          printlog,logfile,file,' NOT FOUND IN apcor.lst'
        endif  

      endif  ; good files

    end  ; loop through all fits files


  ; No aperture correction file
  endif else begin
    successarr[*] = 0             ; all failed
    printlog,logfile,'NO apcor.lst FILE'
  endelse

; DAOGROW did NOT Converge
;-------------------------
endif else begin
  successarr[*] = 0               ; all failed
  printlog,logfile,lines,/logonly
  printlog,logfile,'DAOGROW DID NOT CONVERGE'
  CD,curdir
endelse



;##########################################
;#  UPDATING LIST FILES
;##########################################
undefine,outlist,successlist,failurelist

; We want APCOR.output and APCOR.success to be
; newly REWRITTEN with only the CURRENTLY successful
; files.  Previous success are unimportant since the
; apcor.lst gets remade each time.

; lists2 will make photred_updatelists think that there
; were not outputlines/successlines.
lists2 = lists
lists2.noutputlines = 0
lists2.nsuccesslines = 0

; Success List
ind = where(successarr eq 1,nind)
if nind gt 0 then successlist = inputlines[ind] else UNDEFINE,successlist

; Should this be APPENDED or REWRITTEN?
; Previous successes are unimportant since the apcor.lst gets
; remade each time anyway.


; Output List
; Creating the new output array, ALS files
if (nind gt 0) then begin
  outlist = inputlines[ind]
  ;outlist = fitsdirlist[ind]+'/'+fitsbaselist[ind]
endif else UNDEFINE,outlist

; Failure List
bd = where(successarr eq 0,nbd)
if (nbd gt 0) then begin
  failurelist = inputlines[bd]
endif else UNDEFINE,failurelist

PHOTRED_UPDATELISTS,lists2,outlist=outlist,successlist=successlist,$
                    failurelist=failurelist



printlog,logfile,'PHOTRED_APCOR Finished  ',systime(0)

if keyword_set(stp) then stop

end
