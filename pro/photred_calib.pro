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

pro photred_calib,redo=redo,stp=stp

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


; Are we redoing?
doredo = READPAR(setup,'REDO')
if keyword_set(redo) or (doredo ne '-1' and doredo ne '0') then redo=1

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
  End

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

; Looping through the input files
FOR i=0,ninputlines-1 do begin

  longfile = inputlines[i]
  file = FILE_BASENAME(longfile)
  filedir = FILE_DIRNAME(longfile)
  base = FILE_BASENAME(file,'.ast')
  mchfile = base+'.mch'

  printlog,logfile,'CALIBRATING ',file
  printlog,logfile,systime(0)

  CD,filedir

  ; Check that the AST and MCH files exist
  ;---------------------------------------
  asttest = FILE_TEST(longfile)
  if asttest eq 1 then nastlines = FILE_LINES(longfile) else nastlines=0
  mchtest = FILE_TEST(mchfile)
  if mchtest eq 1 then nmchlines = FILE_LINES(mchfile) else nmchlines=0
  if (nmchlines eq 0) or (nastlines eq 0) then begin
    PUSH,failurelist,longfile
    if asttest eq 0 then printlog,logfile,file+' NOT FOUND'
    if asttest eq 1 and nastlines eq 0 then printlog,logfile,file+' HAS 0 LINES'
    if mchtest eq 0 then printlog,logfile,mchfile,' NOT FOUND'
    if mchtest eq 1 and nmchlines eq 0 then printlog,logfile,mchfile,' HAS 0 LINES'
    goto,BOMB
  endif


  ; Check that the individual FITS files exist
  ; and appear in the apcor.lst
  ; the 'a.del' endings were already removed
  LOADMCH,mchfile,alsfiles
  nalsfiles = n_elements(alsfiles)
  for j=0,nalsfiles-1 do begin

    ialsfile = alsfiles[j]
    ibase = FILE_BASENAME(ialsfile,'.als')
    ifitsfile = ibase+'.fits'

    ifitstest = FILE_TEST(ifitsfile)
    igd = where(apcor.name eq ibase,nigd)
    if (ifitstest eq 0) or (nigd eq 0) then begin
      PUSH,failurelist,longfile
      if ifitstest eq 0 then printlog,logfile,ifitsfile,' NOT FOUND'
      if nigd eq 0 then printlog,logfile,ialsfile,' NOT in apcor.lst'
      goto,BOMB
    endif
  end


  ; DDO51 Radial Offset Correction
  ;-------------------------------
  undefine,photfile
  If keyword_set(ddo51radoffset) and (telescope eq 'BLANCO') and (instrument eq 'MOSAIC') then begin

    ; Check that we have DDO51
    alsbases = FILE_BASENAME(alsfiles,'.als')
    alsfitsfiles = alsbases+'.fits'
    filters = PHOTRED_GETFILTER(alsfitsfiles)

    ; We have a DDO51 filter
    gd_ddo51 = where(filters eq 'D',ngd_ddo51)
    if (ngd_ddo51 gt 0) then begin

      ; Load the AST file
      ast = IMPORTASCII(longfile,/header,/noprint)
      tags = tag_names(ast)

      ; Make sure we have RA/DEC in the AST structure
      if TAG_EXIST(ast,'RA') eq 0 then begin
        printlog,logfile,file+' DOES NOT HAVE RA.  CANNOT DO DDO51 Radial Offset Correction'
        PUSH,failurelist,longfile
        goto,BOMB
      endif
      if TAG_EXIST(ast,'DEC') eq 0 then begin
        printlog,logfile,file+' DOES NOT HAVE DEC.  CANNOT DO DDO51 Radial Offset Correction'
        PUSH,failurelist,longfile
        goto,BOMB
      endif

      ; Get field center
      ishortname = first_el(strsplit(base,'-',/extract))
      gdshort = where(fieldcenters.field eq ishortname,ngdshort)
      if (ngdshort eq 0) then begin
        printlog,logfile,'NO FIELD CENTER FOR '+ishortname
        PUSH,failurelist,longfile
        goto,BOMB
      endif
      ifieldcenters = fieldcenters[gdshort[0]]

      printlog,logfile,'Applying DDO51 Radial Offset Correction'

      ; Convert from RA/DEC to X/Y
      ROTSPHCEN,ast.ra,ast.dec,ifieldcenters.ra,ifieldcenters.dec,xi,eta,/gnomic
      ; The MOSAIC camera is oriented so that N is to the right.
      ; So xi -> YB, and eta -> XB.
      yb = xi*3600./0.26   ; convert from deg to pixels
      xb = eta*3600./0.26
      rad = sqrt(xb^2.0 + yb^2.0)
      xb = xb + 4096
      yb = yb + 4096

      ; Calculate the offset
      ; expr = 'P[0]*exp(-0.5*(X-P[1])^2.0/P[2]^2.0)+P[3]+P[4]*X'
      ; 0.0632479 593.665 600.155 0.0186573 -9.44991e-06
      ddo51_radoffset = 0.0632479*exp(-0.5*(rad-593.665)^2.0/(600.155^2.0)) + 0.0186573 - 9.44991e-06*rad
      ddo51_radoffset = -ddo51_radoffset     ; we want to remove the offset by addition

      ; Add to the structure
      ADD_TAG,ast,'RPIX',0.0,ast
      ADD_TAG,ast,'XB',0.0,ast
      ADD_TAG,ast,'YB',0.0,ast
      ADD_TAG,ast,'DDO51_RADOFFSET',0.0,ast
      ast.rpix = rad
      ast.xb = xb
      ast.yb = yb
      ast.ddo51_radoffset = ddo51_radoffset


      ; Now correct the magnitudes
      for j=0,ngd_ddo51-1 do begin
        ; The columns are: ID, X, Y, MAG1, MAG1ERR, MAG2, MAG2ERR, ...
        thismag = gd_ddo51[j]
        thiscol = 3+2*thismag
        ;thiscol = where(tags eq 'MAG'+strtrim(thismag+1,2),nthiscol)

        ; Applying the correction
        ast.(thiscol) = ast.(thiscol) + ddo51_radoffset

        ; bad values are still bad
        bdval = where(ast.(thiscol) gt 90,nbdval)
        if nbdval gt 0 then ast[bdval].(thiscol)=99.9999

      end ; DDO51 filter loop

      ; Write the new file  
      photfile = base+'.temp'
      PRINTSTR,ast,photfile,/silent

    endif else begin
      printlog,logfile,'NO DDO51 Filter.  CANNOT APPLY DDO51 Radial Offset'
      undefine,photfile
    endelse

  Endif  ; DDO51 radial offset


  ; Make the input file with PHOTRED_PHOTCALIB_PREP
  ;------------------------------------------------
  inpfile = base+'.input'
  printlog,logfile,'Making input file: ',inpfile
  PHOTRED_PHOTCALIB_PREP,mchfile,apcor,inpfile,error=error,imager=thisimager,$
                         observatory=observatory,photfile=photfile

  ; Make sure the input file exists
  inptest = FILE_TEST(inpfile)
  if (inptest eq 0) or (n_elements(error) ne 0) then begin
    PUSH,failurelist,longfile
    if (inptest eq 0) then printlog,logfile,inpfile,' NOT FOUND' else $
        printlog,logfile,inpfile,' ERROR'
    goto,BOMB
  endif


  ;*****************************************************
  ; DOES THE TRANS FILE HAVE THE MAGNITUDES WE NEED???
  ;*****************************************************
  ; I think PHOTCALIB deals with this

  ; Check info there's CHIP information in the trans eqns.
  ;  CHIP=-1 means there's no info
  chipinfo = 0
  if tag_exist(trans,'CHIP') then begin
    gchip = where(trans.chip ge 1,ngchip)
    if ngchip gt 0 then chipinfo=1
  endif

  ; CHIP-SPECIFIC transformation equations
  if chipinfo eq 1 then begin
    inptransfile = ''
    ext = first_el(strsplit(base,thisimager.separator,/extract),/last)
    chip = long(ext)
    ind = where(trans.chip eq chip,nind)

    ; Nothing for this chip
    if nind eq 0 then begin
      PUSH,failurelist,longfile
      printlog,logfile,'NO transformation information for CHIP=',strtrim(chip,2)
      goto,BOMB
    endif
    inptrans = trans[ind]

  ; Global transformation equations
  endif else begin
    undefine,inptrans
    inptransfile = transfile
  endelse

  ; Run PHOTCALIB
  ;---------------
  printlog,logfile,'Calibrating photometry for ',file
  PHOTCALIB,inpfile,inptransfile,inptrans=inptrans,average=doavgmag,keepinstrumental=dokeepinstr,$
            onlyaverage=doavgonlymag,logfile=logfile,/header


  ; Check that the PHOT file exists
  ;-----------------------------------
  photfile = base+'.phot'
  phottest = FILE_TEST(photfile)
  if phottest eq 1 then nlines = FILE_LINES(photfile) else nlines=0
  if (phottest eq 1) and (nlines gt 0) then begin
    PUSH,outlist,filedir+'/'+photfile     ; add to outputarr
    PUSH,successlist,longfile
  endif else begin
    PUSH,failurelist,longfile
    if phottest eq 0 then printlog,logfile,photfile,' NOT FOUND'
    if phottest eq 1 and nlines eq 0 then printlog,logfile,photfile,' HAS 0 LINES'
  endelse


  BOMB:

  CD,curdir


  ;##########################################
  ;#  UPDATING LIST FILES
  ;##########################################
  PHOTRED_UPDATELISTS,lists,outlist=outlist,successlist=successlist,$
                      failurelist=failurelist,/silent

  ;stop

ENDFOR



;#####################
; SUMMARY of the Lists
;#####################
PHOTRED_UPDATELISTS,lists,outlist=outlist,successlist=successlist,$
                    failurelist=failurelist


printlog,logfile,'PHOTRED_CALIB Finished  ',systime(0)

if keyword_set(stp) then stop

end
