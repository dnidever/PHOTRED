;+
;
; PHOTRED_CALIB
;
; This does the calibration (i.e runs photcalib.pro) for photred
; See photcalib_prep.pro and photcalib.pro
;
; INPUTS:
;  /redo Redo files that were already done.
;  /stp  Stop at the end of the program.
;
; OUTPUTS:
;  The calibrated photometry PHOT files from PHOTCALIB
;
; By D.Nidever  Mar 2008
;-

pro photred_calib,nmulti=nmulti,redo=redo,stp=stp

COMMON photred,setup

print,''
print,'########################'
print,'RUNNING PHOTRED_CALIB'
print,'########################'
print,''

CD,current=curdir

; Does the logs/ directory exist?
testlogs = file_test('logs',/directory)
if testlogs eq 0 then FILE_MKDIR,'logs'

; Log files
;----------
thisprog = 'CALIB'
logfile = 'logs/'+thisprog+'.log'
logfile = FILE_EXPAND_PATH(logfile)  ; want absolute filename
if file_test(logfile) eq 0 then SPAWN,'touch '+logfile,out


; Print date and time to the logfile
printlog,logfile,''
printlog,logfile,'Starting PHOTRED_'+thisprog+'  ',systime(0)


; Check that all of the required programs are available
progs = ['readline','readlist','readpar','photred_photcalib_prep','photred_getfilter',$
         'photred_getexptime','photcalib','photred_getinput','photred_updatelists','photred_loadsetup',$
         'undefine','push','printlog','importascii','loadmch','photred_photcalib_prep','photcalib',$
         'photred_getairmass','photred_getdate','badpar','airmass','first_el','strsplitter',$
         'touchzero','writeline','mktemp','photred_getuttime','sexig2ten','stress','strep',$
         'rotsph','rotsphcen','robust_mean','mad','wmeanerr','read_trans']
test = PROG_TEST(progs)
if min(test) eq 0 then begin
  bd = where(test eq 0,nbd)
  printlog,logfile,'SOME NECESSARY PROGRAMS ARE MISSING:'
  printlog,logfile,progs[bd]
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

; Catalog format to use
catformat = READPAR(setup,'catformat')
if catformat eq '0' or catformat eq '' or catformat eq '-1' then catformat='ASCII'
if catformat ne 'ASCII' and catformat ne 'FITS' then catformat='ASCII'

; MCHUSETILES
mchusetiles = READPAR(setup,'MCHUSETILES')
if mchusetiles eq '0' or mchusetiles eq '' or mchusetiles eq '-1' then undefine,mchusetiles
tilesep = '+'

; Are we redoing?
doredo = READPAR(setup,'REDO')
if keyword_set(redo) or (doredo ne '-1' and doredo ne '0') then redo=1

; Hyperthread?
hyperthread = READPAR(setup,'hyperthread')
if hyperthread ne '0' and hyperthread ne '' and hyperthread ne '-1' then hyperthread=1
if strtrim(hyperthread,2) eq '0' then hyperthread=0
if n_elements(hyperthread) eq 0 then hyperthread=0

; Getting NMULTI
if n_elements(nmulti) eq 0 then begin
  nmulti = READPAR(setup,'NMULTI')
  nmulti = long(nmulti)
  ; Use NMULTI_WCS if set
  nmultiwcs = READPAR(setup,'NMULTI_WCS')
  if nmultiwcs ne '0' and nmultiwcs ne '' and nmultiwcs ne '-1' then nmulti=long(nmultiwcs)
endif
nmulti = nmulti > 1  ; must be >=1   

; Telescope, Instrument
telescope = READPAR(setup,'TELESCOPE')
telescope = strupcase(telescope)
instrument = READPAR(setup,'INSTRUMENT')
instrument = strupcase(instrument)
observatory = READPAR(setup,'OBSERVATORY')
if observatory eq '0' or observatory eq '-1' or observatory eq '' then undefine,observatory
if strlowcase(telescope) eq 'blanco' then observatory='ctio'
if strlowcase(telescope) eq 'swope' then observatory='lco'
if strlowcase(telescope) eq 'magellan' then observatory='lco'
if strlowcase(telescope) eq 'lbt' then observatory='mgio'
ddo51radoffset = READPAR(setup,'DDO51RADOFFSET')
ddo51radoffset = strtrim(ddo51radoffset,2)
if ddo51radoffset ne '1' then undefine,ddo51radoffset

; Get the scripts directory from setup
scriptsdir = READPAR(setup,'SCRIPTSDIR')
if scriptsdir eq '' then begin
  printlog,logfile,'NO SCRIPTS DIRECTORY'
  return
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


;###################
; GETTING INPUTLIST
;###################
; INLIST         AST files
; OUTLIST        PHOT files
; SUCCESSLIST    AST files

; Get input
;-----------
precursor = 'ASTROM'
lists = PHOTRED_GETINPUT(thisprog,precursor,redo=redo,ext='ast')
ninputlines = lists.ninputlines


; No files to process
;---------------------
if ninputlines eq 0 then begin
  printlog,logfile,'NO FILES TO PROCESS'
  return
endif

inputlines = lists.inputlines



; Getting the transformation filename
transfile = READPAR(setup,'TRANS')
if transfile eq '0' or transfile eq '-1' then begin
  printlog,logfile,'There is NO TRANS file in the PHOTRED.SETUP file'
  return
endif
transfile = FILE_SEARCH(transfile,/fully,count=ntransfile)
if ntransfile eq 0 then begin
  transfile = READPAR(setup,'TRANS')
  printlog,logfile,'Transformation file >>',transfile,'<< NOT FOUND'
  return
endif
; Load the transformation equations
READ_TRANS,transfile,trans,/silent


; Check that the apcor.lst file exists
apcortest = FILE_TEST('apcor.lst')
if (apcortest eq 0) then begin
  printlog,logfile,'apcor.lst NOT FOUND'
  return
endif
apcorfile = 'apcor.lst'
apcorfile = FILE_SEARCH(apcorfile,/fully)

; Getting the aperture correction structure
apcor = IMPORTASCII(apcorfile,fieldnames=['name','value'],/noprint)
; Remove the 'a.del' endings for the names
apcor_orig = apcor
apcor.name = repstr(apcor.name,'a.del','')  ; base names


; Are we keeping the INSTRUMENTAL magnitudes 
dokeepinstr=0
keepinstr = READPAR(setup,'KEEPINSTR')
if (keepinstr ne '-1' and keepinstr ne '0') then dokeepinstr=1
; Are we AVERAGING multiple observations
doavgmag=0 & doavgonlymag=0
avgmag = READPAR(setup,'AVGMAG')
if (avgmag ne '-1' and avgmag ne '0') then doavgmag=1
avgonlymag = READPAR(setup,'AVGONLYMAG')
if (avgonlymag ne '-1' and avgonlymag ne '0') then doavgonlymag=1
; If BOTH are set then do AVGONLY
if keyword_set(doavgmag) and keyword_set(doavgonlymag) then doavgmag=0
; Printing PHOTCALIB Settings
if keyword_set(dokeepinstr) then printlog,logfile,$
  'KEEPING INSTRUMENTAL MAGNITUDES IN OUTPUT FILES'
if keyword_set(doavgmag) then printlog,logfile,$
  'PUTTING INDIVIDUAL MAGNITUDES AND AVERAGE MAGNITUDES IN OUTPUT FILES'
if keyword_set(doavgonlymag) then printlog,logfile,$
  'PUTTING *ONLY* AVERAGE MAGNITUDES IN OUTPUT FILES'


; APPLYING DDO51 Radial OFFSET
; GETTING FIELD CENTERS
;---------------------------------
If keyword_set(ddo51radoffset) and (telescope eq 'BLANCO') and (instrument eq 'MOSAIC') then begin

  printlog,logfile,''
  printlog,logfile,'Applying DDO51 Radial OFFSET'
  printlog,logfile,'GETTING FIELD CENTERS'

  ;; Make sure we have the DDO51 band
  ;bases = FILE_BASENAME(inputlines,'.ast')
  ;filedirs = FILE_DIRNAME(inputlines)
  ;fitsfiles = filedirs+'/'+bases+'.fits'
  ;filters = PHOTRED_GETFILTER(fitsfiles)

  ; Get the fields
  shortfieldarr = strarr(ninputlines)
  For i=0,ninputlines-1 do begin
    sep = '-'  ; '.'
    file = FILE_BASENAME(inputlines[i])
    shortfield = first_el(strsplit(file,sep,/extract))
    shortfieldarr[i] = shortfield
  Endfor

  ui = uniq(shortfieldarr,sort(shortfieldarr))
  ui = ui[sort(ui)]
  sfields = shortfieldarr[ui]
  nsfields = n_elements(sfields)
  printlog,logfile,strtrim(nsfields,2)+' unique fields found'

  fieldcenters = REPLICATE({field:'',ra:999999.0d0,dec:999999.0d0},nsfields)

  ; Loop through each unique field
  For i=0,nsfields-1 do begin

    isfield = sfields[i]

    ;printlog,logfile,'Getting Field Center for '+isfield

    ; Get all the AST files for this field
    astfiles = FILE_SEARCH(isfield+'-*.ast',count=nastfiles)

    ; Load in all the AST files
    undefine,all
    For j=0,nastfiles-1 do begin

      iastfile = astfiles[j]
      iastbase = FILE_BASENAME(iastfile,'.ast')
      ext = first_el(strsplit(iastbase,'_',/extract),/last)
      ext = long(ext)

      ast = IMPORTASCII(astfiles[j],/header,/noprint)
      nast = n_elements(ast)
      temp1 = REPLICATE({ext:0,ra:0.0d0,dec:0.0d0},nast)
      temp1.ext = ext
      temp1.ra = ast.ra
      temp1.dec = ast.dec
      PUSH,all,temp1
    End  ; ast files
    nall = n_elements(all)

    ; No files found
    if (nall eq 0) then begin
      printlog,logfile,'NO AST Files for '+isfield
      goto,FIELDCENTER_BOMB
    endif

    ; Now find the field center
    ; This is not that simple because with SWARP images the center of the images
    ;  is probably not the center of the chip/amp.
    
    ; Half-way between extremes of middle 8 amplifiers
    ; This is the best method if there are enough amplifiers to use
      
    ; RA
    ; middle 8 amps
    graleft = where(all.ext ge 3 and all.ext le 6,ngraleft)
    ; all
    if ngraleft eq 0 then $
      graleft = where(all.ext ge 1 and all.ext le 8,ngraleft)
    graright = where(all.ext ge 11 and all.ext le 14,ngraright)
    ; all
    if ngraright eq 0 then $
      graright = where(all.ext ge 9 and all.ext le 16,ngraright)
    ; DEC
    ; 1st amps from center
    gdecbot = where(all.ext eq 4 or all.ext eq 12,ngdecbot)
    gdectop = where(all.ext eq 5 or all.ext eq 13,ngdectop)
    ; 2nd amps from center
    if (ngdecbot eq 0 or ngdectop eq 0) then begin
      gdecbot = where(all.ext eq 3 or all.ext eq 11,ngdecbot)
      gdectop = where(all.ext eq 6 or all.ext eq 14,ngdectop)
    endif
    ; 3rd amps from center
    if (ngdecbot eq 0 or ngdectop eq 0) then begin
      gdecbot = where(all.ext eq 2 or all.ext eq 10,ngdecbot)
      gdectop = where(all.ext eq 7 or all.ext eq 15,ngdectop)
    endif
    ; 4th amps from center
    if (ngdecbot eq 0 or ngdectop eq 0) then begin
      gdecbot = where(all.ext eq 1 or all.ext eq 9,ngdecbot)
      gdectop = where(all.ext eq 8 or all.ext eq 16,ngdectop)
    endif
      
    if (ngraleft gt 0 and ngraright gt 0 and ngdecbot gt 0 and ngdectop gt 0) then begin
      cra = mean( [ min(all[graleft].ra), max(all[graright].ra) ] )
      cdec = mean( [ min(all[gdecbot].dec), max(all[gdectop].dec) ] )
      
    ; Calculating the field center from each amplifier
    endif else begin
    
      ; Get unique amps
      uiamps = uniq(all.ext,sort(all.ext))
      amps = all[uiamps].ext
      namps = n_elements(amps)
      
      ; This is (ampmnra-cra)*cos(cdec/!radeg)
      amp_rashift = [ -0.15089910, -0.15110530, -0.15431237, -0.15238810,$
                      -0.15440239, -0.15363146, -0.15023996, -0.15083753,$
                       0.15594973,  0.16516746,  0.15622857,  0.15950704,$
                       0.15472859,  0.15436622,  0.15452388,  0.14968334 ]
      ; This is (ampmndec-cdec)
      amp_decshift = [ -0.269097, -0.195312, -0.115852, -0.0405336,$
                        0.0405336, 0.115852, 0.195312, 0.269097,$
                       -0.269097, -0.195312, -0.115852, -0.0405336,$
                        0.0405336, 0.115852, 0.195312, 0.269097 ]
        
      ampmidra = dblarr(namps)
      ampmiddec = dblarr(namps)
      ampcra = dblarr(namps)
      ampcdec = dblarr(namps)
      for j=0,namps-1 do begin
        gdamp = where(all.ext eq amps[j],ngdamp)
        midra = mean( [ min(all[gdamp].ra), max(all[gdamp].ra) ] )
        middec = mean( [ min(all[gdamp].dec), max(all[gdamp].dec) ] )
        ampmidra[j] = midra   
        ampmiddec[j] = middec
        
        ampcdec[j] = middec - amp_decshift[amps[j]-1]
        ampcra[j] = midra - amp_rashift[amps[j]-1]/cos(ampcdec[i]/!radeg)
        ;ampcra[j] = midra - amp_rashift[amps[j]-1]/cos(middec/!radeg)
        ;ampcra[j] = midra - amp_rashift[amps[j]-1]/cos(cdec/!radeg)
      end
                         
      ROBUST_MEAN,ampcra,cra_ampmn
      ROBUST_MEAN,ampcdec,cdec_ampmn
                          
      cra = cra_ampmn
      cdec = cdec_ampmn
        
    Endelse
        
    ; Put in structure
    fieldcenters[i].field = isfield
    fieldcenters[i].ra = cra
    fieldcenters[i].dec = cdec

    ; Print
    printlog,logfile,isfield+' Field Center = RA:'+strtrim(cra,2)+'  DEC:'+strtrim(cdec,2)

    FIELDCENTER_BOMB:

  End ; each field

Endif  ; ddo51 radial offset, getting field centers




;##################################################
;#  PROCESSING THE FILES
;##################################################
printlog,logfile,''
printlog,logfile,'-----------------------'
printlog,logfile,'PROCESSING THE FILES'
printlog,logfile,'-----------------------'
printlog,logfile,''
printlog,logfile,systime(0)

undefine,outlist,successlist,failurelist
undefine,cmd,cmddir

; Looping through the input files
FOR i=0,ninputlines-1 do begin

  longfile = inputlines[i]
  file = FILE_BASENAME(longfile)
  filedir = FILE_DIRNAME(longfile)
  base = FILE_BASENAME(file,'.ast')
  mchfile = base+'.mch'

  cmd1 = "PHOTRED_CALIB_SINGLE,'"+longfile+"','"+scriptsdir+"/imagers','"+transfile+"'"
  cmd1 += ",'"+apcorfile+"',catformat='"+catformat+"',telescope='"+telescope+"'"
  cmd1 += ",instrument='"+instrument+"',observatory='"+observatory+"'"
  if keyword_set(ddo51radialoffset) then begin
    cmd1 += ',/ddo51radialoffset'
    cmd1 += ',fieldra='+strtrim(fieldcenters[i].ra,2)
    cmd1 += ',fielddec='+strtrim(fieldcenters[i].dec,2)
  endif
  if keyword_set(mchusetiles) then cmd1 += ',/mchusetiles'
  if keyword_set(dokeepstr) then cmd1 += ',/keepinstr'
  if keyword_set(doavgmag) then cmd1 += ',/avgmag'
  if keyword_set(doavgonly) then cmd1 += ',/avgonlymag'

  PUSH,cmd,cmd1
  PUSH,cmddir,filedir

ENDFOR
ncmd = n_elements(cmd)

if ncmd gt 0 then begin
  cmd = "cd,'"+cmddir+"' & "+cmd ; go to the directory
  ; Submit the jobs to the daemon
  PBS_DAEMON,cmd,cmddir,nmulti=nmulti,prefix='calib',hyperthread=hyperthread,/idle,$
             waittime=1,/cdtodir,scriptsdir=scriptsdir
endif


;; Check the output files
for i=0,ninputlines-1 do begin
  longfile = inputlines[i]
  file = FILE_BASENAME(longfile)
  filedir = FILE_DIRNAME(longfile)
  base = FILE_BASENAME(file,'.ast')

  ; Check that the PHOT file exists
  ;-----------------------------------
  photfile = filedir+'/'+base+'.phot'
  phottest = FILE_TEST(photfile)
  if phottest eq 1 then nlines = FILE_LINES(photfile) else nlines=0
  if (phottest eq 1) and (nlines gt 0) then begin
    PUSH,outlist,photfile        ; add to outputarr
    PUSH,successlist,longfile
  endif else begin
    PUSH,failurelist,longfile
    if phottest eq 0 then printlog,logfile,photfile,' NOT FOUND'
    if phottest eq 1 and nlines eq 0 then printlog,logfile,photfile,' HAS 0 LINES'
  endelse
endfor


;#####################
; SUMMARY of the Lists
;#####################
PHOTRED_UPDATELISTS,lists,outlist=outlist,successlist=successlist,$
                    failurelist=failurelist,setupdir=curdir


printlog,logfile,'PHOTRED_CALIB Finished  ',systime(0)

if keyword_set(stp) then stop

end
