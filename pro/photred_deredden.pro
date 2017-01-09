;+
;
; PHOTRED_DEREDDEN
;
; This dereddens photometry for photred
; See deredden.pro
;
; INPUTS:
;  /redo Redo files that were already done.
;  /stp  Stop at the end of the program.
;
; OUTPUTS:
;  The final calibrated and dereddened photometry files
;
; By D.Nidever  Mar 2008
;-
pro photred_deredden,redo=redo,stp=stp


COMMON photred,setup

print,''
print,'########################'
print,'RUNNING PHOTRED_DEREDDEN'
print,'########################'
print,''

CD,current=curdir

; Does the logs/ directory exist?
testlogs = file_test('logs',/directory)
if testlogs eq 0 then FILE_MKDIR,'logs'

; Log files
;----------
thisprog = 'DEREDDEN'
logfile = 'logs/'+thisprog+'.log'
logfile = FILE_EXPAND_PATH(logfile)  ; want absolute filename
if file_test(logfile) eq 0 then SPAWN,'touch '+logfile,out


; Print date and time to the logfile
printlog,logfile,''
printlog,logfile,'Starting PHOTRED_'+thisprog+'  ',systime(0)


; Check that all of the required programs are available
progs = ['readline','readlist','readpar','importascii','first_el','strsplitter','combine_structs',$
         'photred_getinput','photred_updatelists','photred_loadsetup','dust_getval','touchzero',$
         'wcs_getval','wcs_coord2pix','push','undefine','printlog','add_tag','printstr','writeline',$
         'bh_rdfort','djs_int2bin','mktemp','stress','wcs_coord2pix','djs_angpos','djs_ceil','strep']
test = PROG_TEST(progs)
if min(test) eq 0 then begin
  bd = where(test eq 0,nbd)
  printlog,logfile,'SOME NECESSARY PROGRAMS ARE MISSING:'
  printlog,logfile,progs[bd]
  return
endif


; Check that DUST_DIR is properly defined in the shell
dustdir = GETENV('DUST_DIR')
if (dustdir eq '') then begin
  print,'DUST_DIR TO SCHLEGEL MAPS NOT DEFINED IN SHELL'
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



;###################
; GETTING INPUTLIST
;###################
; INLIST         CMB files
; OUTLIST        DERED files
; SUCCESSLIST    CMB files

; Get input
;-----------
precursor = 'COMBINE'
lists = PHOTRED_GETINPUT(thisprog,precursor,redo=redo,ext='cmb')
ninputlines = lists.ninputlines


; No files to process
;---------------------
if ninputlines eq 0 then begin
  printlog,logfile,'NO FILES TO PROCESS'
  return
endif

inputlines = lists.inputlines



; Load the "extinction" file
;--------------------------
printlog,logfile,'Loading the >>extinction<< file'
test = FILE_TEST('extinction')
if (test eq 0) then begin
  scriptsdir = READPAR(setup,'SCRIPTSDIR')
  if scriptsdir eq '-1' or scriptsdir eq '0' then begin
    printlog,logfile,'NO SCRIPTSDIR'
    return
  endif
  test2 = FILE_TEST(scriptsdir+'/extinction')
  if (test2 eq 0) then begin
    printlog,logfile,'NO >>extinction<< file in ',scriptsdir
    return
  endif
  FILE_COPY,scriptsdir+'/extinction','.',/overwrite
endif

; Load the extinction values
; EXTRATIO is A(filter)/E(B-V), so A(filter) = EXTRATIO * E(B-V)
extstr = IMPORTASCII('extinction',fieldname=['filter','extratio'],comment='#',/noprint)
nextstr = n_elements(extstr)
if (nextstr eq 0) then begin
  printlog,logfile,'>>extinction<< file is EMPTY'
  return
endif


; What magnitudes and colors are we dereddening
;----------------------------------------------
printlog,logfile,''
printlog,logfile,'Getting the MAGNITUDES and COLORS to deredden'
todered = READPAR(setup,'TODERED')
if todered eq '-1' or todered eq '0' then undefine,todered
ntodered = n_elements(todered)
; Some mags/colors to deredden
if (ntodered gt 0) then begin
  todered = strcompress(todered,/remove_all)    ; remove all whitespace
  toderedarr = strsplit(todered,',',/extract)   ; split
  ntoderedarr = n_elements(toderedarr)

  ; Structure type
  dum = {type:0,input:'',mag1:'',mag2:'',extratio1:0.0,extratio2:0.0}


  ; Loop through the mags/colors to deredden
  for i=0,ntoderedarr-1 do begin

    ideredstr = dum
    ideredstr.input = toderedarr[i]

    ; Is it a color or magnitude
    minus = strpos(toderedarr[i],'-')
    ; Magnitude
    if (minus eq -1) then begin
      ideredstr.type = 1
      ideredstr.mag1 = strtrim(toderedarr[i],2)

      ; Get extinction value
      gdext = where(extstr.filter eq ideredstr.mag1,ngdext)
      if (ngdext eq 0) then begin
        printlog,logfile,'TODERED - ',toderedarr[i],'  NO EXTINCTION VALUE FOR FILTER ',ideredstr.mag1
        goto,EXT_BOMB
      endif
      ideredstr.extratio1 = extstr[gdext[0]].extratio
      
      ; Add to structure
      PUSH,toderedstr,ideredstr


    ; Color
    endif else begin
      ideredstr.type = 0

      magarr = strsplit(toderedarr[i],'-',/extract)
      magarr = strtrim(magarr,2)
      nmagarr = n_elements(magarr)

      ; Can only have two magnitudes in the color
      if (nmagarr gt 2) then begin
        printlog,logfile,'TODERED - ',ideredstr.input,'  HAS MORE THAN 2 MAGNITUDES IN THE COLOR'
        goto,EXT_BOMB
      endif

      ideredstr.mag1 = magarr[0]
      ideredstr.mag2 = magarr[1]

      ; Get extinction values, Magnitude1
      gdext1 = where(extstr.filter eq ideredstr.mag1,ngdext1)
      if (ngdext1 eq 0) then begin
        printlog,logfile,'TODERED - ',ideredstr.input,'  NO EXTINCTION VALUE FOR FILTER ',ideredstr.mag1
        goto,EXT_BOMB
      endif
      ideredstr.extratio1 = extstr[gdext1[0]].extratio

      ; Get extinction values, Magnitude2
      gdext2 = where(extstr.filter eq ideredstr.mag2,ngdext2)
      if (ngdext2 eq 0) then begin
        printlog,logfile,'TODERED - ',ideredstr.input,'  NO EXTINCTION VALUE FOR FILTER ',ideredstr.mag2
        goto,EXT_BOMB
      endif
      ideredstr.extratio2 = extstr[gdext2[0]].extratio

      ; Add to structure
      PUSH,toderedstr,ideredstr

    endelse

    EXT_BOMB:

  end  ; todered loop


  ; Print out magnitudes/colors to deredden
  ntoderedstr = n_elements(toderedstr)
  for i=0,ntoderedstr-1 do begin

    ; Print header
    if (i eq 0) then begin
      printlog,logfile,''
      printlog,logfile,'Magnitudes/Colors to Deredden'
      printlog,logfile,'------------------------'
      printlog,logfile,'MAG/COLOR  A(M)/E(M1-M2)'
      printlog,logfile,'------------------------'
    endif

    tstr = toderedstr[i]

    ; Magnitude
    if (tstr.type eq 1) then begin
      printlog,logfile,'',tstr.mag1,tstr.extratio1,format='(A3,A-8,F8.3)'
      
    ; Color
    endif else begin
      red = tstr.extratio1-tstr.extratio2
      printlog,logfile,'',tstr.input,red,format='(A3,A-8,F8.3)'

    endelse

    ; Print footer
    if (i eq (ntoderedstr-1)) then begin
      printlog,logfile,'------------------------'
    endif

  end


; No mags/colors to deredden
endif else begin
  printlog,logfile,'NO MAGNITUDES/COLORS to DEREDDEN'
endelse


;stop

; Extinction values to add
;----------------------------
printlog,logfile,''
printlog,logfile,'Getting EXTINCTION values to add'
toextadd = READPAR(setup,'TOEXTADD')
if toextadd eq '-1' or toextadd eq '0' then undefine,toextadd
ntoextadd = n_elements(toextadd)
; Some extinction values to add
if (ntoextadd gt 0) then begin
  toextadd = strcompress(toextadd,/remove_all)    ; remove all whitespace
  toextaddarr = strsplit(toextadd,',',/extract)   ; split
  ntoextaddarr = n_elements(toextaddarr)

  ; Structure type
  dum = {type:0,input:'',mag1:'',mag2:'',extratio1:0.0,extratio2:0.0}


  ; Loop through the mags/colors to deredden
  for i=0,ntoextaddarr-1 do begin

    iextaddstr = dum
    iextaddstr.input = toextaddarr[i]

    ; Is it a color or magnitude
    minus = strpos(toextaddarr[i],'-')
    ; Magnitude
    if (minus eq -1) then begin
      iextaddstr.type = 1
      iextaddstr.mag1 = strtrim(toextaddarr[i],2)

      ; Get extinction value
      gdext = where(extstr.filter eq iextaddstr.mag1,ngdext)
      if (ngdext eq 0) then begin
        printlog,logfile,'TODERED - ',toextaddarr[i],'  NO EXTINCTION VALUE FOR FILTER ',iextaddstr.mag1
        goto,EXT_BOMB2
      endif
      iextaddstr.extratio1 = extstr[gdext[0]].extratio
      
      ; Add to structure
      PUSH,toextaddstr,iextaddstr


    ; Color
    endif else begin
      iextaddstr.type = 0

      magarr = strsplit(toextaddarr[i],'-',/extract)
      magarr = strtrim(magarr,2)
      nmagarr = n_elements(magarr)

      ; Can only have two magnitudes in the color
      if (nmagarr gt 2) then begin
        printlog,logfile,'TODERED - ',iextaddstr.input,'  HAS MORE THAN 2 MAGNITUDES IN THE COLOR'
        goto,EXT_BOMB2
      endif

      iextaddstr.mag1 = magarr[0]
      iextaddstr.mag2 = magarr[1]

      ; Get extinction values, Magnitude1
      gdext1 = where(extstr.filter eq iextaddstr.mag1,ngdext1)
      if (ngdext1 eq 0) then begin
        printlog,logfile,'TODERED - ',iextaddstr.input,'  NO EXTINCTION VALUE FOR FILTER ',iextaddstr.mag1
        goto,EXT_BOMB2
      endif
      iextaddstr.extratio1 = extstr[gdext1[0]].extratio

      ; Get extinction values, Magnitude2
      gdext2 = where(extstr.filter eq iextaddstr.mag2,ngdext2)
      if (ngdext2 eq 0) then begin
        printlog,logfile,'TODERED - ',iextaddstr.input,'  NO EXTINCTION VALUE FOR FILTER ',iextaddstr.mag2
        goto,EXT_BOMB2
      endif
      iextaddstr.extratio2 = extstr[gdext2[0]].extratio

      ; Add to structure
      PUSH,toextaddstr,iextaddstr

    endelse

    EXT_BOMB2:

  end  ; todered loop


  ; Print out magnitudes/colors to deredden
  ntoextaddstr = n_elements(toextaddstr)
  for i=0,ntoextaddstr-1 do begin

    ; Print header
    if (i eq 0) then begin
      printlog,logfile,''
      printlog,logfile,'Extinctions to add'
      printlog,logfile,'------------------------'
      printlog,logfile,'MAG/COLOR  A(M)/E(M1-M2)'
      printlog,logfile,'------------------------'
    endif

    tstr = toextaddstr[i]

    ; Magnitude
    if (tstr.type eq 1) then begin
      printlog,logfile,'',tstr.mag1,tstr.extratio1,format='(A3,A-8,F8.3)'
      
    ; Color
    endif else begin
      red = tstr.extratio1-tstr.extratio2
      printlog,logfile,'',tstr.input,red,format='(A3,A-8,F8.3)'

    endelse

    ; Print footer
    if (i eq (ntoextaddstr-1)) then begin
      printlog,logfile,'------------------------'
    endif

  end


; No extinctions to add
endif else begin
  printlog,logfile,'NO EXTINCTIONS TO ADD'
endelse


; Only adding E(B-V)
ntoderedstr = n_elements(toderedstr)
ntoextaddstr = n_elements(toextaddstr)
if (ntoderedstr eq 0) and (ntoextaddstr eq 0) then begin
  printlog,logfile,'ONLY ADDING E(B-V)'
endif

;stop

;##################################################
;#  PROCESSING THE FILES
;##################################################
printlog,logfile,''
printlog,logfile,'-----------------------'
printlog,logfile,'PROCESSING THE FILES'
printlog,logfile,'-----------------------'
printlog,logfile,systime(0)

;successarr = intarr(ninputlines)-1     ; 0-bad, 1-good
;undefine,outputarr
undefine,outlist,successlist,failurelist

; Looping through the input files
FOR i=0,ninputlines-1 do begin

  longfile = inputlines[i]
  file = FILE_BASENAME(longfile)
  filedir = FILE_DIRNAME(longfile)
  ending = first_el(strsplit(file,'.',/extract),/last)
  base = FILE_BASENAME(file,'.'+ending)  

  printlog,logfile,''
  printlog,logfile,'============================='
  printlog,logfile,'DEREDDENING ',file
  printlog,logfile,'============================='
  printlog,logfile,''
  printlog,logfile,systime(0)
  
  CD,filedir

  ; Make sure it exists
  test = FILE_TEST(file)
  if test eq 1 then nlines = FILE_LINES(file)
  if (test eq 0) or (nlines eq 0) then begin
    ;successarr[i] = 0
    if test eq 0 then printlog,logfile,file,' NOT FOUND'
    if test eq 1 and nlines eq 0 then printlog,logfile,file,' HAS 0 LINES'
    PUSH,failurelist,longfile
    goto,BOMB
  endif

  ; Load the data file
  ;--------------------
  printlog,logfile,'Loading file'
  str = IMPORTASCII(file,/header,/noprint)

  ; Checking that we've got coordinates
  tags = TAG_NAMES(str)
  raind = where(tags eq 'RA',nraind)
  decind = where(tags eq 'DEC',ndecind)
  if (nraind eq 0) or (ndecind) eq 0 then begin
    ;successarr[i] = 0
    printlog,logfile,file,' DOES NOT HAVE RA/DEC'
    PUSH,failurelist,longfile
    goto,BOMB
  endif

  ; Add an EBV tag to the structure
  gdebv = where(tags eq 'EBV',ngdebv)
  if ngdebv eq 0 then $
    ADD_TAG,str,'EBV',0.0,str

  nstars = n_elements(str)
  print,'Dereddening ',strtrim(nstars,2),' sources'


  ; Getting E(B-V) from the Schlegel maps
  ;----------------------------------------
  ra = str.ra
  dec = str.dec

  ; Convert to galactic coordinates
  GLACTC,ra,dec,2000.,glon,glat,1,/degree

  ; Get the E(B-V) reddening from the Schlegel maps
  ebv = DUST_GETVAL(glon,glat,/interp,/noloop)
  str.EBV = EBV



  ; There can be problems with the dereddening if there
  ; are multiple observations per band.  These will have
  ; names of M1, M2, etc. and won't match the filter names.

  ;---------------------------------------
  ; ADDING Dereddened Magnitudes/Colors
  ;---------------------------------------
  ; Loop through the toderedstr structure
  for j=0,ntoderedstr-1 do begin

    tstr = toderedstr[j]

    ; Magnitude
    ;----------
    if (tstr.type eq 1) then begin

      ; Check that the magnitude exists
      ;magind = where(tags eq strupcase(tstr.mag1),nmagind)
      magind = where(tags eq strupcase(tstr.mag1)+'MAG',nmagind)
      if (nmagind eq 0) then begin
        printlog,logfile,'NO Magnitude ',tstr.mag1
        goto,TODERED_BOMB
      endif

      ; Getting the magnitude
      mag = str.(magind[0])

      ; Adding the tag to the structure
      newtagname = strupcase(tstr.mag1)+'0'
      newtagname = IDL_VALIDNAME(newtagname,/convert_all)    ; make sure it's a valid IDL name
      gdnewtagname = where(tags eq newtagname,ngdnewtagname)
      if (ngdnewtagname eq 0) then begin
        ADD_TAG,str,newtagname,0.0,str

      ; Tag already exists, may be overwriting something
      endif else begin
        printlog,logfile,'TAG ',newtagname,' ALREADY EXISTS.  MAY BE OVERWRITTING SOME DATA'
      endelse
      tags = TAG_NAMES(str)
      deredmagind = where(tags eq newtagname,nderedmagind)

      ; Dereddening
      red = tstr.extratio1 * str.ebv
      deredmag = mag - red
      ; Only deredden stars with decent magnitudes
      bdmag = where(mag gt 50.,nbdmag)
      if (nbdmag gt 0) then deredmag[bdmag] = 99.99
      ; Adding to the structure
      str.(deredmagind) = deredmag
      

    ; Color
    ;------
    endif else begin

      ; Check that magnitude1 exists
      ;mag1ind = where(tags eq strupcase(tstr.mag1),nmag1ind)
      mag1ind = where(tags eq strupcase(tstr.mag1)+'MAG',nmag1ind)
      if (nmag1ind eq 0) then begin
        printlog,logfile,'NO Magnitude ',tstr.mag1
        goto,TODERED_BOMB
      endif

      ; Check that magnitude2 exists
      ;mag2ind = where(tags eq strupcase(tstr.mag2),nmag2ind)
      mag2ind = where(tags eq strupcase(tstr.mag2)+'MAG',nmag2ind)
      if (nmag2ind eq 0) then begin
        printlog,logfile,'NO Magnitude ',tstr.mag2
        goto,TODERED_BOMB
      endif

      ; Getting the magnitude
      mag1 = str.(mag1ind[0])
      mag2 = str.(mag2ind[0])


  ; If it is E(B-V) then DO NOT redo it.!!


      ; Adding the tag to the structure
      newtagname = strupcase(tstr.mag1)+strupcase(tstr.mag2)+'0'
      newtagname = IDL_VALIDNAME(newtagname,/convert_all)    ; make sure it's a valid IDL name
      gdnewtagname = where(tags eq newtagname,ngdnewtagname)
      if (ngdnewtagname eq 0) then begin
        ADD_TAG,str,newtagname,0.0,str

      ; Tag already exists, may be overwriting something
      endif else begin
        printlog,logfile,'TAG ',newtagname,' ALREADY EXISTS.  MAY BE OVERWRITTING SOME DATA'
      endelse
      tags = TAG_NAMES(str)
      deredcolind = where(tags eq newtagname,nderedmagind)

      ; Dereddening
      red1 = tstr.extratio1 * str.ebv
      red2 = tstr.extratio2 * str.ebv
      deredcol = (mag1 - red1) - (mag2 - red2)
      ; Only deredden stars with decent magnitudes
      bdmag = where(mag1 gt 50. or mag2 gt 50.,nbdmag)
      if (nbdmag gt 0) then deredcol[bdmag] = 99.99
      ; Adding to the structure
      str.(deredcolind) = deredcol

    endelse

    TODERED_BOMB:

  endfor


  ;---------------------------------------
  ; ADDING EXTINCTIONS
  ;---------------------------------------
  ; Loop through the toderedstr structure
  for j=0,ntoextaddstr-1 do begin

    tstr = toextaddstr[j]

    ; Magnitude
    ;----------
    if (tstr.type eq 1) then begin

      ; Check that the magnitude exists
      ;magind = where(tags eq strupcase(tstr.mag1),nmagind)
      magind = where(tags eq strupcase(tstr.mag1)+'MAG',nmagind)
      if (nmagind eq 0) then begin
        printlog,logfile,'NO Magnitude ',tstr.mag1
        goto,TODERED_BOMB2
      endif

      ; Adding the tag to the structure
      newtagname = 'A'+strupcase(tstr.mag1)
      newtagname = IDL_VALIDNAME(newtagname,/convert_all)    ; make sure it's a valid IDL name
      gdnewtagname = where(tags eq newtagname,ngdnewtagname)
      if (ngdnewtagname eq 0) then begin
        ADD_TAG,str,newtagname,0.0,str

      ; Tag already exists, may be overwriting something
      endif else begin
        printlog,logfile,'TAG ',newtagname,' ALREADY EXISTS.  MAY BE OVERWRITTING SOME DATA'
      endelse
      tags = TAG_NAMES(str)
      extind = where(tags eq newtagname,nextind)

      ; Extinction
      red = tstr.extratio1 * str.ebv
      ; Adding to the structure
      str.(extind) = red
 

    ; Color
    ;------
    endif else begin

      ; Check that magnitude1 exists
      ;mag1ind = where(tags eq strupcase(tstr.mag1),nmag1ind)
      mag1ind = where(tags eq strupcase(tstr.mag1)+'MAG',nmag1ind)
      if (nmag1ind eq 0) then begin
        printlog,logfile,'NO Magnitude ',tstr.mag1
        goto,TODERED_BOMB2
      endif

      ; Check that magnitude2 exists
      ;mag2ind = where(tags eq strupcase(tstr.mag2),nmag2ind)
      mag2ind = where(tags eq strupcase(tstr.mag2)+'MAG',nmag2ind)
      if (nmag2ind eq 0) then begin
        printlog,logfile,'NO Magnitude ',tstr.mag2
        goto,TODERED_BOMB2
      endif

      ; Adding the tag to the structure
      newtagname = 'E'+strupcase(tstr.mag1)+strupcase(tstr.mag2)
      newtagname = IDL_VALIDNAME(newtagname,/convert_all)    ; make sure it's a valid IDL name
      gdnewtagname = where(tags eq newtagname,ngdnewtagname)
      if (ngdnewtagname eq 0) then begin
        ADD_TAG,str,newtagname,0.0,str

      ; Tag already exists, may be overwriting something
      endif else begin
        printlog,logfile,'TAG ',newtagname,' ALREADY EXISTS.  MAY BE OVERWRITTING SOME DATA'
      endelse      
      tags = TAG_NAMES(str)
      extind = where(tags eq newtagname,nextind)

      ; Dereddening
      red1 = tstr.extratio1 * str.ebv
      red2 = tstr.extratio2 * str.ebv
      red = (red1 - red2)
      ; Adding to the structure
      str.(extind) = red

    endelse

    TODERED_BOMB2:

  endfor


  ; WRITE new file
  ;------------------
  ; Output the structure to the DERED file
  deredfile = base+'.dered'
  printlog,logfile,'New file with extinctions is: ',deredfile
  PRINTSTR,str,deredfile

  ; Check that the file DERED file is there
  deredtest = FILE_TEST(deredfile)
  if (deredtest eq 0) then begin
    PUSH,failurelist,longfile
    printlog,logfile,deredfile,' NOT FOUND'
  endif


  ; Add to outlist/successlist
  PUSH,outlist,filedir+'/'+deredfile
  PUSH,successlist,longfile


  BOMB:

  CD,curdir


  ;#####################
  ; UPDATE the Lists
  ;#####################
  PHOTRED_UPDATELISTS,lists,outlist=outlist,successlist=successlist,$
                      failurelist=failurelist,/silent

ENDFOR


;#####################
; SUMMARY of the Lists
;#####################
PHOTRED_UPDATELISTS,lists,outlist=outlist,successlist=successlist,$
                    failurelist=failurelist


printlog,logfile,'PHOTRED_DEREDDEN Finished  ',systime(0)

if keyword_set(stp) then stop

end
