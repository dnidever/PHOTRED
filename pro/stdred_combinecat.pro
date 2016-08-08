pro stdred_combinecat,redo=redo,stp=stp

;+
;
; STDRED_COMBINECAT
;
; This uses the "matched" catalogs from STDRED_MATCHCAT
; and combines all of the observations for a given filter
; into a final catalog that can be used to derive
; photometric transformation equations.
;
; INPUTS:
;  /redo Redo files that were already done.
;  /stp  Stop at the end of the program.
;
; OUTPUTS:
;  CAT files for each frame that contains some standard
;  stars.  The CAT files has the same columns as the
;  AST files plus the calibrated data.
;
; By D.Nidever  May 2008
;-

COMMON photred,setup

print,''
print,'#########################'
print,'RUNNING STDRED_COMBINECAT'
print,'#########################'
print,''

CD,current=curdir

; Does the logs/ directory exist?
testlogs = file_test('logs',/directory)
if testlogs eq 0 then FILE_MKDIR,'logs'

; Log files
;----------
thisprog = 'COMBINECAT'
logfile = 'logs/'+thisprog+'.log'
logfile = FILE_EXPAND_PATH(logfile)  ; want absolute filename
if file_test(logfile) eq 0 then SPAWN,'touch '+logfile,out


; Print date and time to the logfile
printlog,logfile,''
printlog,logfile,'Starting STDRED_'+thisprog+'  ',systime(0)


; Check that all of the required programs are available
progs = ['readline','readlist','readpar','importascii','add_tag','printstr','photred_getinput',$
         'photred_updatelists','photred_loadsetup','push','undefine','printlog','writecol',$
         'photred_getairmass','photred_getfilter','photred_getuttime','strsplitter','touchzero',$
         'writeline','airmass','first_el','mktemp','photred_getdate','sexig2ten','stress','badpar',$
         'strep','robust_mean','mad','wmeanerr']
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
  PHOTRED_LOADSETUP,setup,count=count,/std
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
ddo51radoffset = READPAR(setup,'DDO51RADOFFSET')
ddo51radoffset = strtrim(ddo51radoffset,2)
if ddo51radoffset ne '1' then undefine,ddo51radoffset



;###################
; GETTING INPUTLIST
;###################
; INLIST         PHOT files
; OUTLIST        AST files
; SUCCESSLIST    PHOT files

; Get input
;-----------
precursor = 'MATCHCAT'
lists = PHOTRED_GETINPUT(thisprog,precursor,redo=redo,ext='cat')
ninputlines = lists.ninputlines


; No files to process
;---------------------
if ninputlines eq 0 then begin
  printlog,logfile,'NO FILES TO PROCESS'
  return
endif

inputlines = lists.inputlines




; Get the scripts directory from setup
scriptsdir = READPAR(setup,'SCRIPTSDIR')
if (scriptsdir eq '') then begin
  printlog,logfile,'NO SCRIPTS DIRECTORY'
  return
endif


; How many unique filters are there
filtarr = strarr(ninputlines)
nightarr = intarr(ninputlines)
framearr = strarr(ninputlines)
for i=0,ninputlines-1 do begin
  file = FILE_BASENAME(inputlines[i],'.cat')
  filarr = strsplit(file,'-',/extract)
  prefix = filarr[0]
  prefixarr = strsplit(prefix,'n',/extract)
  filtarr[i] = prefixarr[0]
  nightarr[i] = prefixarr[1]
  framearr[i] = filarr[1]
end

ui = uniq(filtarr,sort(filtarr))
filters = filtarr[ui]
nfilters = n_elements(filters)


; Load the TRANS.SETUP file
printlog,logfile,'Loading the parameters for the transformation equations, >>trans.setup<<'
testtrans = FILE_TEST('trans.setup')
if (testtrans eq 0) then begin
  printlog,logfile,'NO local >>trans.setup<< file found.  Checking scriptsdir'
  testtrans2 = FILE_TEST(scriptsdir+'/trans.setup')
  if (testtrans2 eq 0) then begin
    printlog,logfile,'NO >>trans.setup<< file found in scriptsdir'
    return
  endif else begin
    FILE_COPY,scriptsdir+'/trans.setup','.',/overwrite,/allow
    READCOL,'trans.setup',tmagnames,tcolnames,format='A,A'
  endelse

; Load local trans.setup file
endif else begin
  READCOL,'trans.setup',tmagnames,tcolnames,format='A,A'
endelse
ntmagnames = n_elements(tmagnames)

; Get "shortnames" for the magnitude names
tmagnames2 = PHOTRED_GETFILTER(filtname=tmagnames,/silent,/noupdate)

; Match the magnitude names in the "trans.setup" to the
; color to use
usemag = strarr(nfilters)
usecol = strarr(nfilters,2)
for i=0,nfilters-1 do begin

  ifilter = filters[i]
  gdfilt = where(tmagnames2 eq ifilter,ngdfilt)
  if (ngdfilt gt 0) then begin
    twocol = strsplit(tcolnames[gdfilt[0]],'-',/extract)
    ntwocol = n_elements(twocol)
    if (ntwocol eq 2) then begin
      usemag[i] = tmagnames[gdfilt[0]]
      usecol[i,*] = twocol
    endif else begin
      printlog,logfile,'Color ',tcolnames[gdfilt[0]],' NOT understood'
    endelse

  endif else begin
    printlog,logfile,'NO entry for FILTER=',ifilter,' in >>trans.setup<<'
  endelse

end


; APPLYING DDO51 Radial OFFSET
; GETTING FIELD CENTERS
;---------------------------------
If keyword_set(ddo51radoffset) and (telescope eq 'BLANCO') and (instrument eq 'MOSAIC') then begin

  printlog,logfile,''
  printlog,logfile,'Applying DDO51 Radial OFFSET'
  printlog,logfile,'GETTING FIELD CENTERS'

  ; Make sure we have the DDO51 band
  bases = FILE_BASENAME(inputlines,'.cat')
  gd_ddo51 = where(filters eq 'D',ngd_ddo51)
  
  ; Some DDO51 filter
  If (ngd_ddo51 gt 0) then begin

    ; Get the frame names
    framearr = strsplitter(bases,'_',/extract)
    framearr = reform(framearr[0,*])
    
    filternames = strmid(strtrim(framearr,2),0,1)
    gd_dframes = where(filternames eq 'D',ngd_dframes)
    dframearr = framearr[gd_dframes]

    ui = uniq(dframearr,sort(dframearr))
    ui = ui[sort(ui)]
    dframes = dframearr[ui]
    ndframes = n_elements(dframes)
    printlog,logfile,strtrim(ndframes,2)+' unique D frames found'

    fieldcenters = REPLICATE({frame:'',ra:0.0d0,dec:0.0d0},ndframes)

    ; Loop through each unique frame
    For i=0,ndframes-1 do begin

      iframe = dframes[i]

      ; Get all the AST files for this field
      astfiles = FILE_SEARCH(iframe+'_*.ast',count=nastfiles)

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
          ampcra[j] = midra - amp_rashift[amps[j]-1]/cos(ampcdec[j]/!radeg)
          ;ampcra[j] = midra - amp_rashift[amps[j]-1]/cos(cdec/!radeg)
        endfor

        ROBUST_MEAN,ampcra,cra_ampmn
        ROBUST_MEAN,ampcdec,cdec_ampmn

        cra = cra_ampmn
        cdec = cdec_ampmn

      Endelse

      ; Put in structure
      fieldcenters[i].frame = iframe
      fieldcenters[i].ra = cra
      fieldcenters[i].dec = cdec

      ; Print
      printlog,logfile,iframe+' Frame Center = RA:'+strtrim(cra,2)+'  DEC:'+strtrim(cdec,2)

    End ; each frame

  ; NO DDO51 filters
  Endif else begin
    printlog,logfile,'NO DDO51 Observations'
  Endelse

Endif  ; ddo51 radial offset, getting field centers



;##################################################
;#  PROCESSING THE FILES
;##################################################
printlog,logfile,''
printlog,logfile,'-----------------------'
printlog,logfile,'PROCESSING THE FILES'
printlog,logfile,'-----------------------'

undefine,outlist,successlist,failurelist

; Loop through the filters
FOR i=0,nfilters-1 do begin

  undefine,combcat

  ifilter = filters[i]
  filtind = where(filtarr eq ifilter,nfiltind)
  filtfiles = inputlines[filtind]

  printlog,logfile,''
  printlog,logfile,'================================================='
  printlog,logfile,' Combining standard star data for FILTER = ',ifilter
  printlog,logfile,'================================================='
  printlog,logfile,strtrim(nfiltind,2),' CAT files to combine'
  printlog,logfile,''


  ; Check that we have a color to use for this magnitude
  if (usecol[i,0] eq '') then begin
    printlog,logfile,'NO COLOR to use for FILTER=',ifilter
    goto,BOMB2
  endif

  ; What magnitude and color to use for the transformation equation
  ; for this filter.
  iusemag = usemag[i]
  iusecol = reform(usecol[i,*])
  printlog,logfile,'Transformation magnitude: ',iusemag
  printlog,logfile,'Transformation color: ',iusecol[0],'-',iusecol[1]


  ; How many nights are there
  uinight = uniq(nightarr[filtind],sort(nightarr[filtind]))
  nights = nightarr[filtind[uinight]]
  nnights = n_elements(nights)

  printlog,logfile,strtrim(nnights,2),' nights of data: ',strjoin('n'+strtrim(nights,2),', ')
  printlog,logfile,''


  ; Loop through the nights
  For n=0,nnights-1 do begin

    inight = nights[n]
    filtnightind = where(filtarr eq ifilter and nightarr eq inight,nfiltnightind)
    files = inputlines[filtnightind]


    ; Loop through the files for this filter/night pair
    For j=0,nfiltnightind-1 do begin

      longfile = files[j]
      file = FILE_BASENAME(longfile)
      base = FILE_BASENAME(file,'.cat')
      filedir = FILE_DIRNAME(longfile)
      fitsfile = base+'.fits'

      iframe = framearr[filtnightind[j]]
      ichip = long( (strsplit(base,'_',/extract))[1] )

      printlog,logfile,'Adding ',file


      CD,filedir

      ; Test that the CAT file exists
      test = FILE_TEST(file)
      if test eq 1 then nlines = FILE_LINES(file)
      if (test eq 0) or (nlines eq 0) then begin
        PUSH,failurelist,longfile
        if test eq 0 then printlog,logfile,file,' NOT FOUND'
        if test eq 1 and nlines eq 0 then printlog,logfile,file,' HAS 0 LINES'
        goto,BOMB
      endif

      ; Test that the FITS file exists
      fitstest = FILE_TEST(fitsfile)
      if (fitstest eq 0) then begin
        PUSH,failurelist,longfile
        printlog,logfile,fitsfile,' NOT FOUND'
        goto,BOMB
      endif

      ; Get MJD night number
      mjd = PHOTRED_GETMJD(fitsfile,observatory)
      if (mjd le 0) then begin
        PUSH,failurelist,longfile
        printlog,logfile,'NO MJD for ',fitsfile
        goto,BOMB
      endif
      
      ; Get UT-TIME
      uttime = PHOTRED_GETUTTIME(fitsfile)
      if (uttime eq '') then begin
        PUSH,failurelist,longfile
        printlog,logfile,'NO UT-TIME for ',fitsfile
        goto,BOMB
      endif

      ; Get EXPTIME
      exptime = PHOTRED_GETEXPTIME(fitsfile)
      if (exptime le 0.0) then begin
        PUSH,failurelist,longfile
        printlog,logfile,'NO EXPTIME for ',fitsfile
        goto,BOMB
      endif

      ; Get AIRMASS
      airmass = PHOTRED_GETAIRMASS(fitsfile,obs=observatory)
      if (airmass lt 0.9) then begin
        PUSH,failurelist,longfile
        printlog,logfile,'NO AIRMASS for ',fitsfile
        goto,BOMB
      endif

      ; Loading the CAT file
      ;---------------------
      cat = IMPORTASCII(file,/header,/noprint)
      ncat = n_elements(cat)
      printlog,logfile,'Nstars = ',strtrim(ncat,2)

      ; DDO51 Radial Offset Correction
      ;-------------------------------
      undefine,photfile
      If (ifilter eq 'D') and keyword_set(ddo51radoffset) and (telescope eq 'BLANCO') and (instrument eq 'MOSAIC') then begin

        ; Get field center
        iframename = first_el(strsplit(base,'_',/extract))
        gdframe = where(fieldcenters.frame eq iframename,ngdframe)
        if (ngdframe eq 0) then begin
          printlog,logfile,'NO FIELD CENTER FOR '+iframename
          PUSH,failurelist,longfile
          goto,BOMB
        endif
        ifieldcenters = fieldcenters[gdframe[0]]

        printlog,logfile,'Applying DDO51 Radial Offset Correction'

        ; Convert from RA/DEC to X/Y
        ROTSPHCEN,cat.ra,cat.dec,ifieldcenters.ra,ifieldcenters.dec,xi,eta,/gnomic
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
        ADD_TAG,cat,'RPIX',0.0,cat
        ADD_TAG,cat,'XB',0.0,cat
        ADD_TAG,cat,'YB',0.0,cat
        ADD_TAG,cat,'DDO51_RADOFFSET',0.0,cat
        cat.rpix = rad
        cat.xb = xb
        cat.yb = yb
        cat.ddo51_radoffset = ddo51_radoffset

        ; Now correct the magnitudes
        cat.mag = cat.mag + ddo51_radoffset

        ; bad values are still bad
        bdval = where(cat.mag gt 90,nbdval)
        if nbdval gt 0 then cat[bdval].mag=99.9999

      Endif  ; DDO51 radial offset


      ; Figure out what calibrated magnitudes are in the CAT file
      ;----------------------------------------------------------
      ; Use the error tags
      cat_tags = TAG_NAMES(cat)
      cerrind = where(stregex(cat_tags,'^C_',/boolean) eq 1 and $
                      (stregex(cat_tags,'ERR$',/boolean) eq 1 or $
                      stregex(cat_tags,'_ERR_',/boolean) eq 1),ncerr)
      cerr_tags = cat_tags[cerrind]
      cerrnames = strmid(cerr_tags,2,20)
      cmagnames = strarr(ncerr)
      len = strlen(cerrnames)
      for k=0,ncerr-1 do begin
        ; C_RERR, CMERR, style name
        if stregex(cerrnames[k],'ERR$',/boolean) then begin
          cmagnames[k] = strmid(cerrnames[k],0,len[k]-3)
        ; C_ERR_R, C_ERR_U style name
        endif else begin
          dum = strsplit(cerrnames[k],'_',/extract)
          cmagnames[k] = strtrim(dum[1],2)
        endelse
      endfor
      tmag_tags = 'C_'+cmagnames

      ; Check if the magnitudes corresponding to the error tags exist
      ; and have short filternames
      tcatmagshnames = strarr(ncerr)
      tcatmagind = intarr(ncerr)-1
      for k=0,ncerr-1 do begin
        tmagind = where(cat_tags eq tmag_tags[k],ntmagind)
        if (ntmagind gt 0) then begin
          tcatmagshnames[k] = PHOTRED_GETFILTER(filtname=cmagnames[k],/noupdate,/silent,/fold_case)
          tcatmagind[k] = tmagind[0]       ; This is the tag index
        endif
      end
      tcatmagshnames = strupcase(tcatmagshnames)
      ; Get the magnitudes that we can use
      gdcatmags = where(tcatmagshnames ne '' and tcatmagind ge 0,ngdcatmags)
      if (ngdcatmags gt 0) then begin
        catmagnames = cmagnames[gdcatmags]
        catmagshnames = tcatmagshnames[gdcatmags]
        catmagind = tcatmagind[gdcatmags]
        caterrind = cerrind[gdcatmags]

      ; No useful magnitudes found in this cat file
      endif else begin
        printlog,logfile,'NO useful magnitudes found in the ',file,' CAT file'
        PUSH,failurelist,longfile
        goto,BOMB
      endelse


      ; Which calibrated magnitude and color do we need for this FILTER
      ;----------------------------------------------------------------
      havemagind = where(catmagshnames eq strupcase(iusemag),nhavemagind)
      havecol1ind = where(catmagshnames eq strupcase(iusecol[0]),nhavecol1ind)
      havecol2ind = where(catmagshnames eq strupcase(iusecol[1]),nhavecol2ind)

      if (nhavemagind eq 0 or nhavecol1ind eq 0 or nhavecol2ind eq 0) then begin
        if (nhavemagind eq 0) then printlog,logfile,'Magnitude ',iusemag,' NOT FOUND in CAT file'
        if (nhavecol1ind eq 0) then printlog,logfile,'Magnitude ',iusecol[0],' NOT FOUND in CAT file'
        if (nhavecol2ind eq 0) then printlog,logfile,'Magnitude ',iusecol[1],' NOT FOUND in CAT file'
        PUSH,failurelist,longfile
       goto,BOMB
      endif

      ; These are the tag indices to use
      usemagind = catmagind[havemagind[0]]
      usemagerrind = caterrind[havemagind[0]]
      usecol1ind = catmagind[havecol1ind[0]]
      usecol1errind = caterrind[havecol1ind[0]]
      usecol2ind = catmagind[havecol2ind[0]]
      usecol2errind = caterrind[havecol2ind[0]]


      ; What names to use for the calibrated magnitude and color
      ; use the short names
      ;usemagname = catmagshnames[havemagind[0]]
      ;usemagerrname = usemagname+'ERR'
      ;usecolname = catmagshnames[havecol1ind[0]]+catmagshnames[havecol2ind[0]]
      ;usecolerrname = usecolname+'ERR'
      usemagname = iusemag[0]
      usemagerrname = usemagname+'ERR'
      usecolname = iusecol[0]+'_'+iusecol[1]
      usecolerrname = usecolname+'ERR'



      ; Make the combined standard star data structure
      ;
      ; Need to pick out the calibrated magnitude, color and errors
      ; that we need, and IDs.
      ;
      ; frame+cat.id, cat.c_id, cat.mag, cat.err, cat.c_mag,
      ; cat.c_magerr, cat.c_color, cat.c_colerr
      ; night

      dum = {frame:'',frameid:'',id:'',night:0L,mag:0.0,err:0.0}
      dum = CREATE_STRUCT(dum,usemagname,0.0)
      dum = CREATE_STRUCT(dum,usemagerrname,0.0)
      dum = CREATE_STRUCT(dum,usecolname,0.0)
      dum = CREATE_STRUCT(dum,usecolerrname,0.0)
      dum = CREATE_STRUCT(dum,'airmass',0.0,'nightcount',0L,'mjd',0L,'ut','','weight',0.0,'chip',0L,'exptime',0.0)
      dum = CREATE_STRUCT(dum,'ra',0.0d0,'dec',0.0d0,'stdfield','')
      if TAG_EXIST(cat,'RPIX') then $
        dum=CREATE_STRUCT(dum,'RPIX',0.0,'XB',0.0,'YB',0.0,'DDO51_RADOFFSET',0.0)
      temp = REPLICATE(dum,ncat)

      ; ID ind
      cidind = where(cat_tags eq 'C_ID',ncidind)
      if ncidind eq 0 then cidind=where(cat_tags eq 'C_OBJID',ncidind)
      if ncidind eq 0 then begin
        print,'NO ID found in CAT file'
        PUSH,failurelist,longfile
        goto,BOMB
      endif
      cidind = cidind[0]

      ; Stuff in the information
      temp.frame = iframe
      temp.frameid = iframe+'-'+strtrim(cat.id,2)
      temp.id = strtrim(cat.(cidind),2)
      temp.night = inight
      temp.mag = cat.mag + 2.5*alog10(exptime)   ; CORRECT for exposure time!!!
      temp.err = cat.err
      temp.chip = ichip
      temp.exptime = exptime
      temp.ra = cat.ra
      temp.dec = cat.dec
      temp.stdfield = cat.stdfield

      ; Get the calibrated photometry
      col1 = cat.(usecol1ind)
      col1err = cat.(usecol1errind)
      col2 = cat.(usecol2ind)
      col2err = cat.(usecol2errind)
      color = col1 - col2
      colerr = sqrt( col1err^2.0 + col2err^2.0 )

      ; Stuff into the structure
      temp.(6) = cat.(usemagind)      ; calib. magnitude
      temp.(7) = cat.(usemagerrind)   ; calib. magnitude error
      temp.(8) = color                ; calib. color
      temp.(9) = colerr               ; calib. color error


      ; Calculate the weights
      merr1 = temp.err > 0.001        ; instrumental magnitude error
      merr2 = temp.(7) > 0.001        ; calib. magnitude error
      weight = 0.0014142135624 / sqrt(merr1^2. + merr1^2.)
      ; From FILDAOPHOT.PRO
      ; Evaluate the weight as inverse of quadrature addition of measured error
      ; and quoted error in magnitudes.  The function is normalized so that the
      ; highest possible wt = 1.00
      ;if (merr1 lt 0.001) then merr1 = 0.001
      ;if (magerr1 lt 0.001) then magerr1 = 0.001
      ;wt = 0.0014142135624 / sqrt(merr1^2. + magerr1^2.)

      ; How many observations are there for this night already
      nthisnight = 0
      if n_elements(combcat) gt 0 then dumb = where(combcat.night eq inight,nthisnight)

      ; Extra info
      temp.airmass = airmass

      temp.nightcount = lindgen(ncat)+nthisnight+1
      temp.mjd = mjd
      temp.ut = uttime
      temp.weight = weight

      ; Add DDO51 radial offset correction information
      if TAG_EXIST(temp,'RPIX') then begin
        temp.rpix = cat.rpix
        temp.xb = cat.xb
        temp.yb = cat.yb
        temp.ddo51_radoffset = cat.ddo51_radoffset
      endif

      ; Only use good calibrate photometry
      gdphot = where(temp.mag gt 0 and temp.mag lt 50 and $    ; instrumental mag
                     temp.(6) gt 0 and temp.(6) lt 50 and $    ; calibrated mag
                     temp.(8) gt -5 and temp.(8) lt 5,ngdphot,comp=bdphot,ncomp=nbdphot)   ; calibrated color
      if ngdphot gt 0 and nbdphot gt 0 then begin
        printlog,logfile,'Removing '+strtrim(nbdphot,2)+' stars with BAD photometry'
        temp = temp[gdphot]
      endif
      if ngdphot eq 0 then begin
        printlog,logfile,'NO stars with GOOD photometry'
        PUSH,failurelist,longfile
        goto,BOMB        
      endif

      ; Add to the combined structure for this filter
      PUSH,combcat,temp


      ; Successful
      PUSH,successlist,longfile

      BOMB:

      CD,curdir


      ;##########################################
      ;#  UPDATING LIST FILES
      ;##########################################
      PHOTRED_UPDATELISTS,lists,outlist=outlist,successlist=successlist,$
                        failurelist=failurelist,/silent


    Endfor ; files loop

  Endfor ; night loop


  ; We have a combined CAT structure
  ncombcat = n_elements(combcat)
  if (ncombcat gt 0) then begin

    ; Save the combined filter CAT file
    ;----------------------------------
    printlog,logfile,''
    printlog,logfile,strtrim(ncombcat,2),' Observations for FILTER=',ifilter
    filtcatfile = ifilter+'.cat'
    printlog,logfile,'Combined CAT file for FILTER=',ifilter,' written to "',filtcatfile,'"'
    PRINTSTR,combcat,filtcatfile,/silent

    PUSH,outlist,curdir+'/'+filtcatfile

    ; Adding the night count information
    ; DN 1/22/2014, I DON'T UNDERSTAND WHAT THE POINT OF THESE LINES IS!!
    ;uinight = uniq(combcat.night,sort(combcat.night))
    ;nights = combcat[uinight].night
    ;temp.nightcount = lindgen(ncat)+n_elements(combcat)+1


    ; Make the "skawdphot" FILTER.data file
    ;--------------------------------------
    ;readcol,ddfile,starname,realmag,realcol,secz,instmag,night,num,UT,wt,$
    ;  format='a,d,d,d,d,i,i,a,d'
    ; frameid, cmag, col, airmass, mag, night, nightcount, ut, weight
    skawdfile = ifilter+'.data'
    printlog,logfile,'Making "skawdphot" file "',skawdfile,'"'
    WRITECOL,skawdfile,combcat.frameid,combcat.(6),combcat.(8),combcat.airmass,$
             combcat.mag,combcat.night,combcat.nightcount,combcat.ut,combcat.weight,$
             fmt='(A14,F13.5,F12.5,F12.5,F12.5,I5,I5,A12,F12.5)'

  ; No observations for this filter
  endif else begin
    printlog,logfile,''
    printlog,logfile,'NO Observations for FILTER=',ifilter
  endelse



  BOMB2:

  CD,curdir

  ;##########################################
  ;#  UPDATING LIST FILES
  ;##########################################
  PHOTRED_UPDATELISTS,lists,outlist=outlist,successlist=successlist,$
                    failurelist=failurelist,/silent


ENDFOR ; filter loop



;#####################
; SUMMARY of the Lists
;#####################
PHOTRED_UPDATELISTS,lists,outlist=outlist,successlist=successlist,$
                    failurelist=failurelist


printlog,logfile,'STDRED_COMBINECAT Finished  ',systime(0)

if keyword_set(stp) then stop


end
