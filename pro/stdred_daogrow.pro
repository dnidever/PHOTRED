pro stdred_daogrow,redo=redo,stp=stp

;+
;
; STDRED_DAOGROW
;
; This runs DAOGROW on all of the aperture photometry files
; for each night to get the "total" magnitude for each star.
;
; INPUTS:
;  /redo Redo files that were already done.
;  /stp  Stop at the end of the program.
;
; OUTPUTS:
;  For each .ap a .tot file is created.
;
; By D.Nidever  May 2008
;-

COMMON photred,setup


print,''
print,'########################'
print,'RUNNING STDRED_DAOGROW'
print,'########################'
print,''

CD,current=curdir

; Does the logs/ directory exist?
testlogs = file_test('logs',/directory)
if testlogs eq 0 then FILE_MKDIR,'logs'

; Log files
;----------
thisprog = 'DAOGROW'
logfile = 'logs/'+thisprog+'.log'
logfile = FILE_EXPAND_PATH(logfile)  ; want absolute filename
if FILE_TEST(logfile) eq 0 then SPAWN,'touch '+logfile,out


; Print date and time to the logfile
printlog,logfile,''
printlog,logfile,'Starting STDRED_'+thisprog+'  ',systime(0)


; Check that all of the required programs are available
progs = ['readline','readlist','readpar','photred_getfilter','photred_getuttime',$
         'photred_getexptime','photred_getinput','photred_updatelists','photred_loadsetup',$
         'undefine','push','printlog','loadinput','photred_getairmass','photred_getdate',$
         'badpar','airmass','apcorrect','loadaper','writecol','writeline','first_el','minloc',$
         'mktemp','strsplitter','touchzero','sexig2ten']
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
precursor = 'APERPHOT'
lists = PHOTRED_GETINPUT(thisprog,precursor,redo=redo,ext='ap')
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

; Initialzing some arrays
UNDEFINE,outlist,successlist,failurelist

; Get the scripts directory from setup
scriptsdir = READPAR(setup,'SCRIPTSDIR')
if scriptsdir eq '' then begin
  printlog,logfile,'NO SCRIPTS DIRECTORY'
  return
endif


; How many nights are there?
base = FILE_BASENAME(inputlines,'.ap')
nightarr = intarr(ninputlines)
for i=0,ninputlines-1 do begin
  prefix = first_el(strsplit(base[i],'-',/extract))
  prefixarr = strsplit(prefix,'n',/extract)
  nightarr[i] = fix(prefixarr[1])
end

ui = uniq(nightarr,sort(nightarr))
nights = nightarr[ui]
nnights = n_elements(nights)


; Loop through the nights
;-------------------------
FOR i=0,nnights-1 do begin

  thisnight = nights[i]
  nightname = 'daogrow-n'+strtrim(thisnight,2)
  nightdir = FILE_EXPAND_PATH(nightname)
  ; Should I change directory name to 'daogrow-n#' ??

  printlog,logfile,''
  printlog,logfile,'----------------------------------------------------------'
  printlog,logfile,'Getting TOTAL APERTURE PHOTOMETRY for all NIGHT=',strtrim(thisnight,2),' frames'
  printlog,logfile,'----------------------------------------------------------'
  printlog,logfile,''

  ; Make the directory
  test = FILE_TEST(nightdir,/directory)
  ; directory exists already, erase and start again
  if test eq 1 then FILE_DELETE,nightdir,/recursive
  FILE_MKDIR,nightdir


  ; Copy the scripts to the night directory
  FILE_COPY,scriptsdir+'/daogrow.sh',nightdir,/overwrite
  FILE_COPY,scriptsdir+'/photo.opt',nightdir,/overwrite

  ; What files belong this night?
  nightind = where(nightarr eq thisnight,nnightind)


  ;-----------------------------------------------------
  ; Copying the .AP file to n#/
  ; and Checking that the FITS file exists
  ;-----------------------------------------------------
  printlog,logfile,'Copying AP files to '+nightname+'/'
  successarr = intarr(nnightind)-1          ; 0-bad, 1-good
  for j=0,nnightind-1 do begin

    longfile = inputlines[nightind[j]]
    file = FILE_BASENAME(longfile)
    filedir = FILE_DIRNAME(longfile)
    base = FILE_BASENAME(longfile,'.ap')
    fitsfile = base+'.fits'
    apfile = file
    aapfile = base+'a.ap'

    CD,filedir


    ; Check that the FITS file exists
    ;--------------------------------
    fitstest = FILE_TEST(fitsfile)
    if (fitstest eq 0) then begin
      successarr[j]=0
      PUSH,failurelist,longfile
      printlog,logfile,file,' NOT FOUND'
      goto,BOMB
    endif

    ; Check that the .AP and A.AP files exists
    ;-----------------------------------------
    aptest = FILE_TEST(apfile)
    aaptest = FILE_TEST(aapfile)
    if (aptest eq 0) or (aaptest eq 0) then begin
      successarr[i] = 0
      PUSH,failurelist,longfile
      if aptest eq 0 then printlog,logfile,apfile,' NOT FOUND'
      if aaptest eq 0 then printlog,logfile,aapfile,' NOT FOUND'
      goto,BOMB
    endif

    ; Check .AP header
    ;-----------------
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
      successarr[j]=0
      PUSH,failurelist,longfile
      printlog,logfile,apfile,' DOES HOT HAVE A PROPER AP HEADER'
      goto,BOMB
    endif


    ; Check A.AP header
    ;------------------
    ;This is what it should look like:
    ; NL    NX    NY  LOWBAD HIGHBAD  THRESH     AP1  PH/ADU  RNOISE    FRAD
    ;  2  2040  2047    79.9 30000.0   32.07    3.00    2.20    2.86    2.00
    line1='' & line2='' & line3=''   
    openr,unit,/get_lun,aapfile
    readf,unit,line1
    readf,unit,line2 
    readf,unit,line3
    close,unit
    free_lun,unit

    dum1 = strtrim(strsplit(line1,' ',/extract),2)
    dum2 = strtrim(strsplit(line2,' ',/extract),2)
    line3 = strtrim(line3,2)
    if (dum1[0] ne 'NL') or (dum2[0] ne '2') or (line3 ne '') then begin
      successarr[j]=0
      PUSH,failurelist,longfile
      printlog,logfile,aapfile,' DOES HOT HAVE A PROPER AP HEADER'
      goto,BOMB
    endif


    ; Copy the A.AP files to daogrow/
    ;--------------------------------
    FILE_COPY,aapfile,nightdir,/overwrite,/allow_same

    ; This file was okay
    successarr[j] = 1

    BOMB:

    CD,curdir

  end

  ncopy = long(total(successarr))
  printlog,logfile,strtrim(ncopy,2),'/',strtrim(nnightind,2),' successfully copied to '+nightname+'/'


  ; UPDATE the Lists
  PHOTRED_UPDATELISTS,lists,outlist=outlist,successlist=successlist,$
                      failurelist=failurelist,/silent



  ;--------------------------------
  ; Making the DAOGROW input files
  ;--------------------------------
  printlog,logfile,'Making DAOGROW input files: '+nightname+'/'+nightname+'.inf and '+$
                   nightname+'/'+nightname+'.ext'

  OPENW,infunit,/get_lun,nightdir+'/'+nightname+'.inf'
  OPENW,extunit,/get_lun,nightdir+'/'+nightname+'.ext'

  ; Loop through the AP files
  For j=0,nnightind-1 do begin

    ; Only use "good" files
    if (successarr[j] eq 1) then begin

      longfile = inputlines[nightind[j]]
      file = FILE_BASENAME(longfile)
      filedir = FILE_DIRNAME(longfile)
      base = FILE_BASENAME(longfile,'.ap')
      fitsfile = base+'.fits'
      apfile = file

        ; Getting information from header
      head = HEADFITS(fitsfile)

      ; Getting FILTER
      filt = PHOTRED_GETFILTER(fitsfile,/numeric)

      ; Getting UT
      time = PHOTRED_GETUTTIME(fitsfile)
      if (time eq '') then begin
        successarr[j] = 0
        PUSH,failurelist,longfile
        printlog,logfile,'NO UT-TIME'
        goto,BOMB2
      endif
      timarr = strsplit(time,':',/extract)
      uthr = timarr[0]
      utmin = timarr[1]

      ; Getting AIRMASS
      am = PHOTRED_GETAIRMASS(fitsfile,obs=observatory,/update,/recalculate)
      if am lt 0.9 then begin
        successarr[j] = 0
        PUSH,failurelist,longfile
        printlog,logfile,'NO AIRMASS'
        goto,BOMB2
      end

      ; Getting EXPTIME
      exp = PHOTRED_GETEXPTIME(fitsfile)
      if strtrim(exp,2) eq '-1' then begin
        successarr[j] = 0
        PUSH,failurelist,longfile
        printlog,logfile,'NO EXPTIME'
        goto,BOMB2
      endif


      ; How long is the name
      len = strlen(base+'a')
      if len gt 23 then begin
        printlog,logfile,base,' IS TOO LONG.  MAX 23 characters'
        return
      endif


      ; Outputting to the INF file
      format='(A1,A-23,A10,A4,A3,F7.3,F10.3)'
      printf,infunit,format=format,'',base+'a',filt,uthr,utmin,am,exp

      ; Outputting to the EXT file
      printf,extunit,base+'a.ap'


    endif  ; "good" file


    BOMB2:

  Endfor  ; AP files loop

  CLOSE,infunit
  FREE_LUN,infunit
  CLOSE,extunit
  FREE_LUN,extunit


  ; UPDATE the Lists
  PHOTRED_UPDATELISTS,lists,outlist=outlist,successlist=successlist,$
                      failurelist=failurelist,/silent


  ;-----------------------
  ; Run daogrow.sh script
  ;-----------------------
  CD,nightdir

  print,'RUNNING DAOGROW in '+nightname+'/'
  ; daogrow.sh only uses stars with err<0.2
  SPAWN,'./daogrow.sh '+nightname+' > '+nightname+'.log',out,errout


  ; Run daogrow on all the a.ap files, that have the stars picked by
  ; PSFPIXK.
  ; Use the aperture correction to correct all the stars in the .ap
  ; files.
  ; Get the cumulative aperture correction and its error from the
  ; dagrow.gro file.
  ; Use the aperture with the minimum error,
  ; err = sqrt( aperr^2.0 + cumerr^2.)
  ; since the aperture goes UP with radius and the correction error
  ; goes down with radius, it might not be the first aperture.


  ; DAOGROW Converged
  ;----------------------
  ; Read in the log file
  READLINE,nightname+'.log',lines
  ind = first_el(where(stregex(lines,'converged',/boolean,/fold_case) eq 1,nind),/last)
  if nind gt 0 then begin

    printlog,logfile,lines[ind[0]:*]
    printlog,logfile,'DAOGROW CONVERGED'


    ;-------------------------------------------
    ; Making aperture corrected photometry files
    ;-------------------------------------------

    ; Check the .GRO and .INF files
    inftest = FILE_TEST(nightname+'.inf')
    grotest = FILE_TESt(nightname+'.gro')
    if (inftest eq 0 or grotest eq 0) then begin
      if inftest eq 0 then printlog,logfile,nightname+'.inf NOT FOUND'
      if grotest eq 0 then printlog,logfile,nightname+'.gro NOT FOUND'
      successarr[*] = 0               ; all failed
      PUSH,failurelist,inputlines[nightind]
      goto,BOMB_NIGHT
    endif

    ; Load the .INF file
    READCOL,nightname+'.inf',infnames,format='A',/silent

    ; Load the .GRO file
    READLINE,nightname+'.gro',grolines


    ; Loop through all AP files
    for j=0,nnightind-1 do begin

      ; Only check "good" files
      if (successarr[j] eq 1) then begin

        longfile = inputlines[nightind[j]]
        file = FILE_BASENAME(longfile)
        filedir = FILE_DIRNAME(longfile)
        base = FILE_BASENAME(longfile,'.ap')
        fitsfile = base+'.fits'
        apfile = file
        totfile = filedir+'/'+base+'.tot'
        totafile = base+'a.tot'

        ; What number are we in the .INF file
        infind = where(infnames eq base+'a',ninfind)
        num = string(infind[0]+1,format='(I04)')

        ; Check the A.TOT file
        totatest = FILE_TEST(totafile)
        if totatest eq 1 then totalines = FILE_LINES(totafile) else totalines=0
      
        ; Load the .AP file
        LOADAPER,longfile,aper,aperhead
        aper_orig = aper

        ; Get only stars with "good" photometry
        gd = where(aper.err[0] lt 0.3,ngd)
        aper = aper[gd]
        nstars = n_elements(aper)

        ; Apply the aperture correction
        printlog,logfile,'Applying aperture Correction to >>'+base+'.ap<<'
        undefine,final,error
        APCORRECT,aper,final,grofile=nightname+'.gro',gronum=infind[0]+1,error=error
        nfinal = n_elements(final)
        printlog,logfile,''


        ; We have final photometry
        if (nfinal gt 1 and n_elements(error) eq 0 and totalines ge 2) then begin

          ; Write to file
          WRITECOL,totfile,final.id,final.x,final.y,final.mag,final.err,final.sky,$
                   final.magfap,final.apcorr,final.finalap,fmt='(I7,2F9.3,2F9.4,F9.2,F9.3,F9.4,I9)'
                   ;final.magfap,final.apcorr,final.finalap,fmt='(I7,2F9.3,2F9.4,2F9.3,F9.4,I9)'

          ; Load the A.TOT header
          line1='' & line2=''
          openr,unit,/get_lun,totafile
          readf,unit,line1
          readf,unit,line2
          close,unit
          free_lun,unit
          tothead = [line1,line2,'']

          ; Prepend to the TOT file
          WRITELINE,totfile,tothead,/prepend

          ; Final check of TOT file
          tottest = FILE_TEST(totfile)
          if tottest eq 1 then totlines=FILE_LINES(totfile) else totlines=0

          ; Successful
          if (totlines ge 4) then begin

            PUSH,successlist,longfile
            PUSH,outlist,totfile

          ; Failure
          endif else begin
            printlog,logfile,totfile+' PROBLEMS'
            PUSH,failurelist,longfile
          endelse

        ; Problems
        endif else begin
          if nfinal eq 0 or n_elements(error) gt 0 then printlog,logfile,base+' aperture correction problems'
          if totalines lt 2 then printlog,logfile,totafile,' PROBLEMS'
          PUSH,failurelist,longfile
        endelse

      endif                      ; good files

    endfor ; loop through all AP files


  ; DAOGROW did NOT Converge
  ;-------------------------
  endif else begin
    successarr[*] = 0               ; all failed
    PUSH,failurelist,inputlines[nightind]
    printlog,logfile,lines,/logonly
    printlog,logfile,'DAOGROW DID NOT CONVERGE.  All files failed.'
  endelse

  BOMB_NIGHT:


  ; CD back to original directory
  CD,curdir


  ;#####################
  ; UPDATE the Lists
  ;#####################
  PHOTRED_UPDATELISTS,lists,outlist=outlist,successlist=successlist,$
                      failurelist=failurelist,/silent

END   ; night loop


;#####################
; SUMMARY of the Lists
;#####################
PHOTRED_UPDATELISTS,lists,outlist=outlist,successlist=successlist,$
                    failurelist=failurelist


printlog,logfile,'STDRED_DAOGROW Finished  ',systime(0)

if keyword_set(stp) then stop

end
