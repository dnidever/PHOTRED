;+
;
; COMPLETENESS
;
; Figure out completeness from artificial starsn for a single PHOTRED "field".
;
; Put together list of original detected source + input artificial
; stars. Then match those up with the final list of real sources
; and artificial stars.  The data from all chips and all mocks are combined.
;
; INPUTS:
;  photfiles     List of AST phot files.
;  =imager       The imager information structure.
;  =logfile      A logfile to print output to.
;  =maindir      The main directory for this PHOTRED run.
;  /redo         Recreate AST files if they already exist.
;
; OUTPUTS:
;  The structure of artificial stars and information on whether
;  they were recovered and with what magnitudes is saved.
;  =error     The error message if there was one.
;
; USAGE:
;  IDL>completeness,photfiles,imager=imager,maindir=maindir
; 
; By D.Nidever  based on getcomplete.pro   July 2017
;-

pro completeness,photfiles,imager=imager,logfile=logfile,$
                 maindir=maindir,error=error,redo=redo

undefine,error
  
; Not enough inputs
nphotfiles = n_elements(photfiles)
if nphotfiles eq 0 or n_elements(imager) eq 0 or n_elements(maindir) eq 0 then begin
  error = 'Not enough inputs'
  print,'Syntax - completeness,photfiles,imager=imager,maindir=maindir,logfile=logfile,error=error,redo=redo'
  return
endif

; Logfile
if keyword_set(logfile) then logf=logfile else logf=-1

; Get the "short" field name, e.g. F4
;  The names should look like this: F4M6-00478573_22.phot
allfieldmockname = reform( (strsplitter(file_basename(photfiles),'-',/extract))[0,*] )
allfieldname = reform( (strsplitter(allfieldmockname,'M',/extract))[0,*] )
uifieldname = uniq(allfieldname,sort(allfieldname))
if n_elements(uifieldname) gt 1 then begin
  error = 'ERROR: files for more than one field input'
  printlog,logf,error
  return
endif
shfield = allfieldname[0]  ; they should all be the same

; Get the chip numbers
allchips = lonarr(nphotfiles)
for i=0,nphotfiles-1 do allchips[i]=PHOTRED_GETCHIPNUM(repstr(photfiles[i],'.phot','.fits'),imager)
uichips = uniq(allchips,sort(allchips))
uchips = allchips[uichips]
nchips = n_elements(uchips)


; Load "fields" file to get "global" field name
READCOL,maindir+'/fields',shnames,lnames,format='A,A',/silent
nameind = where(shnames eq shfield,nnameind)
if nnameind eq 0 then begin
  error = 'Short field name >>'+shfield+'<< not found in "fields" file'
  printlog,logf,error
  return
endif
globalfield = lnames[nameind[0]]

printlog,logf,'Getting completeness function for ',shfield,' = ',globalfield
printlog,logf,'-------------------------------------------------------------'


; Loop over the separate chips
;-----------------------------
undefine,bigast
For i=0,nchips-1 do begin

  ichip = uchips[i]
  chipind = where(allchips eq ichip,nmocks)
  if i gt 0 and nchips gt 1 then printlog,logf,' '
  printlog,logf,'  '+strtrim(i+1,2)+'/'+strtrim(nchips,2)+' CHIP='+strtrim(ichip,2)+' - '+strtrim(nmocks,2)+' mock(s)'
  printlog,logf,systime(0)

  ; Check if the chip output file already exists
  photdir = file_dirname(photfiles[chipind[0]])
  photbase = file_basename(photfiles[chipind[0]],'.phot')  
  refname = (strsplit(photbase,'-'+imager.separator,/extract))[1]
  outchipfile = photdir+'/'+shfield+'-'+refname+imager.separator+string(ichip,format='(i02)')+'_complete.fits'
  if file_test(outchipfile+'.gz') eq 1 and not keyword_set(redo) then begin
    print,outchipfile,' already EXISTS and /redo NOT set.  Just loading in the existing file.'
    chipast = MRDFITS(outchipfile+'.gz',1,/silent)
    goto,chipcombine
  endif
  
  ; Loop over the mocks for this chip
  ;----------------------------------
  undefine,chipast
  For j=0,nmocks-1 do begin
    mockphotfile = photfiles[chipind[j]]
    mockphotdir = file_dirname(mockphotfile)
    mockphotbase = file_basename(mockphotfile,'.phot')
    mockfitsfile = mockphotdir+'/'+mockphotbase+'.fits'
    ; Get the mock number for the filename, e.g. F2M5- gives mocknum=5
    fieldmockname = first_el(strsplit(mockphotbase,'-',/extract))
    mocknum = long(first_el(strsplit(fieldmockname,'M',/extract),/last))
    if j gt 0 and nmocks gt 1 then printlog,logf,' '
    printlog,logf,'  MOCK='+strtrim(mocknum,2)
    ; Load the FITS header
    head = headfits(mockfitsfile)
    
    ; Get the reference exposure number
    ;  00379732 in F4M4-00379732_22.phot
    refname = (strsplit(mockphotbase,'-'+imager.separator,/extract))[1]
    
    ; Check phot file header to see if ALLFRAME was run (see if PROB exists)
    READLINE,mockphotfile,phothead,nlineread=1
    photcols = strtrim(strsplit(phothead[0],' ',/extract),2)
    probind = where(strupcase(photcols) eq 'PROB',nprobind)
    if nprobind gt 0 then allframe=1 else allframe=0
    ; Load the PHOT file
    printlog,logf,'  Loading '+mockphotfile
    phot = IMPORTASCII(mockphotfile,/header,/silent)
    nphot = n_elements(phot)
    
    ; -- Add average magnitudes to PHOT file --
    ; We need average magnitues per band for each object
    phottags = tag_names(phot)
    photavgmag = where(stregex(phottags,'MAG$',/boolean) eq 1,nphotavgmag)
    if nphotavgmag eq 0 then begin
      ; Find all the magnitude columns
      photmagind = where(stregex(photcols,'MAG',/boolean) eq 1 and $
                         stregex(photcols,'^I_',/boolean) eq 0 and $
                         stregex(photcols,'ERR',/boolean) eq 0,nphotmagind)
      if nphotmagind eq 0 then begin
        printlog,logf,'NO photometric magnitudes found in '+mockphotfile
        goto,BOMB1 
      endif
      allphotmags = photcols[photmagind]  ; iMAG3, zMAG
      allphotmagsband = reform((strsplitter(allphotmags,'MAG',/extract))[0,*])  ; just the band name
      ; Unique magnitudes
      uiphotmagsband = uniq(allphotmagsband,sort(allphotmagsband))
      uphotmagsband = allphotmagsband[uiphotmagsband]
      nphotmagsband = n_elements(uphotmagsband)
      ; Add the columns to PHOT
      photschema = phot[0]
      STRUCT_ASSIGN,{dum:0},photschema  ; blank it out
      for k=0,nphotmagsband-1 do $
         photschema=CREATE_STRUCT(photschema,uphotmagsband[k]+'MAG',99.99,uphotmagsband[k]+'ERR',9.99,$
                                  'NDET'+uphotmagsband[k],0)
      temp = phot & undefine,phot
      phot = REPLICATE(photschema,nphot)  ; put in new structure format
      STRUCT_ASSIGN,temp,phot,/nozero
      phottags = tag_names(phot)  ; new phot tags
      
      ; Loop over unique bands and average
      for k=0,nphotmagsband-1 do begin
        indphmag = where(allphotmagsband eq uphotmagsband[k],nindphmag)  ; indices into ALLPHOTMAGS/PHOTCOLS
        phmagcols = photcols[photmagind[indphmag]]
        MATCH,phottags,phmagcols,tagindphmag,ind2,/sort  ; indices for 
        tempmag = fltarr(nphot,nindphmag)
        temperr = tempmag
        for l=0,nindphmag-1 do begin
          MATCH,phottags,phmagcols[l],tagindphmag1,ind2,/sort
          tempmag[*,l]=phot.(tagindphmag1)
          temperr[*,l]=phot.(tagindphmag1+1)  ; ERR is always the next one
        endfor
          stop
        ; Average the mags
        AVERAGEMAG,tempmag,temperr,newmag,newerr,/robust
        magind = where(phottags eq uphotmagsband[k]+'MAG',nmagind)
        errind = where(phottags eq uphotmagsband[k]+'ERR',nerrind)
        phot.(magind) = newmag
        phot.(errind) = newerr
        ndet = long(total(tempmag lt 50,2))
        detind = where(phottags eq 'NDET'+strupcase(uphotmagsband[k]),ndetind)
        phot.(detind) = ndet
        phot.ndet += ndet
     endfor                     ; unique band loop      

    ; Add NDETX columns
    endif else begin

      ; Find all the magnitude columns
      photmagind = where(stregex(photcols,'MAG',/boolean) eq 1 and $
                         stregex(photcols,'^I_',/boolean) eq 0 and $
                         stregex(photcols,'ERR',/boolean) eq 0,nphotmagind)
      if nphotmagind eq 0 then begin
        printlog,logf,'NO photometric magnitudes found in '+mockphotfile
        goto,BOMB1 
      endif
      allphotmags = photcols[photmagind]  ; iMAG3, zMAG
      allphotmagsband = reform((strsplitter(allphotmags,'MAG',/extract))[0,*])  ; just the band name
      ; Unique magnitudes
      uiphotmagsband = uniq(allphotmagsband,sort(allphotmagsband))
      uphotmagsband = allphotmagsband[uiphotmagsband]
      nphotmagsband = n_elements(uphotmagsband)
      ; Add the columns to PHOT
      photschema = phot[0]
      STRUCT_ASSIGN,{dum:0},photschema  ; blank it out
      for k=0,nphotmagsband-1 do $
         photschema=CREATE_STRUCT(photschema,'NDET'+uphotmagsband[k],0)
      photschema=CREATE_STRUCT(photschema,'NDET',0)
      
      temp = phot & undefine,phot
      phot = REPLICATE(photschema,nphot)  ; put in new structure format
      STRUCT_ASSIGN,temp,phot,/nozero
      phottags = tag_names(phot)  ; new phot tags

      ; Loop over unique bands and average
      for k=0,nphotmagsband-1 do begin
        indphmag = where(allphotmagsband eq uphotmagsband[k] and $   ; indices into ALLPHOTMAGS/PHOTCOLS
                         allphotmags ne uphotmagsband[k],nindphmag)        ; exclude avg mag
        phmagcols = photcols[photmagind[indphmag]]
        MATCH,phottags,phmagcols,tagindphmag,ind2,/sort  ; indices for 
        tempmag = fltarr(nphot,nindphmag)
        for l=0,nindphmag-1 do begin
          MATCH,phottags,strupcase(phmagcols[l]),tagindphmag1,ind2,/sort
          tempmag[*,l]=phot.(tagindphmag1)
        endfor
        ndet = long(total(tempmag lt 50,2))
        detind = where(phottags eq 'NDET'+strupcase(uphotmagsband[k]),ndetind)
        phot.(detind) = ndet
        phot.ndet += ndet
     endfor               ; unique band loop      
      
    endelse     ; no avg phot mags

    ; Load the final, recovered photometry file, but use the mag/raw
    ;  instrumental photometry.  Use this to match up recovered to
    ;  original objects.
    ;  Use .raw or .mag file depending if allframe was used.
    if allframe eq 1 then magext='.mag' else magext='.raw'
    finalfile = repstr(mockphotfile,'.phot','.mag')
    if (file_info(finalfile)).exists eq 0 then begin
      printlog,logf,finalfile+' NOT FOUND'
      goto,BOMB1
    endif
    if allframe eq 1 then LOADMAG,finalfile,final else LOADRAW,finalfile,final
    nfinal = n_elements(final)
    printlog,logf,'  NFINAL='+strtrim(nfinal,2)
    
    ; Load the original "real" star data file
    ;  Use .raw or .mag file depending if allframe was used.
    origfile = mockphotdir+'/'+shfield+'-'+refname+imager.separator+string(ichip,format='(i02)')+magext
    if (file_info(origfile)).exists eq 0 then begin
      printlog,logf,origfile+' NOT FOUND'
      goto,BOMB1
    endif
    if allframe eq 1 then LOADMAG,origfile,orig else LOADRAW,origfile,orig
    norig = n_elements(orig)
    printlog,logf,'  NORIG='+strtrim(norig,2)
    
    ; Load the artificial star data file, with information on
    ;   the injected stars
    synthfile = mockphotdir+'/'+fieldmockname+'-'+refname+'-add'+imager.separator+string(ichip,format='(i02)')+magext
    if (file_info(synthfile)).exists eq 0 then begin
      printlog,logf,synthfile+' NOT FOUND'
      goto,BOMB1
    endif
    if allframe eq 1 then LOADMAG,synthfile,synth else LOADRAW,synthfile,synth
    nsynth = n_elements(synth)
    printlog,logf,'  NSYNTH='+strtrim(nsynth,2)
    ; Remove stars with NO good magnitudes, were not actually used
    stags = tag_names(synth)
    smagind = where(stregex(stags,'^MAG',/boolean) eq 1 and stregex(stags,'ERR$',/boolean) eq 0,nsmagind)
    badmask = lonarr(n_elements(synth))+1
    for k=0,nsmagind-1 do badmask = badmask AND (synth.(smagind[k]) gt 50)
    sbdmag = where(badmask eq 1,nsbdmag)
    if nsbdmag gt 0 then begin
      printlog,logf,'  Removing '+strtrim(nsbdmag,2)+' synthetic stars with NO good magnitudes'
      REMOVE,sbdmag,synth
      nsynth = n_elements(synth)
    endif
    
    ; Get the calibrated photometry from the CALMAG synthetic star file.
    ;  This has calibrated, apparent magnitudes for the synthetic
    ;  stars for this chip+mock.
    calsynthfile = mockphotdir+'/'+fieldmockname+'-calmag'+imager.separator+string(ichip,format='(i02)')+'.mag'
    if (file_info(calsynthfile)).exists eq 0 then begin
      printlog,logf,synthfile+' NOT FOUND'
      goto,BOMB1
    endif
    calsynth = IMPORTASCII(calsynthfile,/header,count=ncalsynth,/silent)
    printlog,logf,'  NCALSYNTH='+strtrim(ncalsynth,2)

    ; Do the CROSS-MATCHING between the three lists (final/orig/synth)
    ;-----------------------------------------------------------------
    ; If two stars land right on top of each other and only one is recovered then use the magnitude
    ; to match it up to the correct star (real or artificial)
    ; We don't want to think we recovered an artificial star when it's
    ; actually a real star (likely bright one).
    ; We are using the instrumental photomery (.mag/.raw) to do this matching.
    printlog,logf,'  Matching'
    COMPLETENESS_CROSSMATCH,final,orig,synth,find,sind,nmatch=nmatch,logfile=logf,error=xmatcherror
    if n_elements(xmatcherror) gt 0 then goto,BOMB1
    
    ; Put all AST information into one structure
    ;-------------------------------------------

    ; --- Create the structure schema ---
    astschema = {astid:'',photid:'',recovered:-1,field:'',chip:0L,mock:0L,$
              inp_x:0.0d0,inp_y:0.0d0,inp_ra:0.0d0,inp_dec:0.0d0}
    ; Add columns for input synthetic calibrated photometry
    calsynthtags = tag_names(calsynth)
    calsynthphotcolind = where(calsynthtags ne 'ID',ncalsynthphotcolind)
    calsynthphotcols = calsynthtags[calsynthphotcolind]
    for k=0,ncalsynthphotcolind-1 do astschema=CREATE_STRUCT(astschema,'INP_'+calsynthphotcols[k],0.0)
    ; Add recovered information
    astschema = CREATE_STRUCT(astschema,'ID','','X',999999.0d0,'Y',999999.d0,'RA',999999.0d0,'DEC',999999.0d0,'NDET',0L)
    ; Add columms from PHOT file
    ;   get all of the unique filters, anything with MAG and not I_ and
    ;   not with ERR
    photmagind = where(stregex(photcols,'MAG',/boolean) eq 1 and $
                       stregex(photcols,'^I_',/boolean) eq 0 and $
                       stregex(photcols,'ERR',/boolean) eq 0,nphotmagind)
    if nphotmagind eq 0 then begin
      printlog,logf,'NO photometric magnitudes found in '+mockphotfile
      goto,BOMB1 
    endif
    allphotmags = photcols[photmagind]  ; iMAG3, zMAG
    allphotmags = reform((strsplitter(allphotmags,'MAG',/extract))[0,*])  ; just the mag name
    uiphotmags = uniq(allphotmags,sort(allphotmags))
    uphotmags = allphotmags[uiphotmags]
    nphotmags = n_elements(uphotmags)
    ;  Add three columns for each filter: magnitude, error, Ndetections
    for k=0,nphotmags-1 do $
       astschema=CREATE_STRUCT(astschema,uphotmags[k],99.99,uphotmags[k]+'ERR',9.99,'NDET'+uphotmags[k],0L)
    ; Add final morphology columns
    astschema = CREATE_STRUCT(astschema,'chi',999999.0,'sharp',999999.0,'flag',-1L,'prob',999999.0)
    asttags = tag_names(astschema)

    ; --- Make the structure and copy over the data ---
    ast = replicate(astschema,nsynth)
    ; Copy over the input synth data
    ;   ASTID is a unique name across all mocks and chips
    ast.astid = shfield+'.'+strtrim(ichip,2)+'.'+strtrim(mocknum,2)+'.'+strtrim(synth.id,2)
    ast.photid = strtrim(synth.id,2)
    ast.field = shfield
    ast.chip = ichip
    ast.mock = mocknum
    ast.inp_x = synth.x
    ast.inp_y = synth.y
    ; Add the RA/DEC values using the input X/Y coordinates
    ;  and the reference header
    HEAD_XYAD,head,ast.inp_x-1,ast.inp_y-1,inpra,inpdec,/deg
    ast.inp_ra = inpra
    ast.inp_dec = inpdec
    
    ; Stick in the calibration synthetic input photometry from CALSYNTH    
    MATCH2,calsynth.id,synth.id,ind1,ind2
    dum = where(ind2 gt -1,ncalmatch)
    if ncalmatch ne nsynth then begin
      printlog,logf,'Not all synth elements found matches in calsynth'
      goto,BOMB1
    endif
    ; Loop over calsynth columns and copy
    for k=0,ncalsynthphotcolind-1 do begin
      asttagind = where(asttags eq 'INP_'+calsynthphotcols[k],nasttagind)
      calsynthtagind = where(calsynthtags eq calsynthphotcols[k],ncalsynthtagind)
      ast.(asttagind) = calsynth[ind2].(calsynthtagind)
    endfor
    ; Copy over the recovered values using STRUCT_ASSIGN
    ;  LEFT was created from the PHOT structure
    ;temp = ast[aind2]
    ;STRUCT_ASSIGN,left[aind1],temp,/nozero
    ;ast[aind2] = temp
    ;ast[aind2].recovered = 1
    temp = ast[sind]
    STRUCT_ASSIGN,phot[find],temp,/nozero
    temp.id = strtrim(temp.id,2)
    for k=0,nphotmags-1 do begin
       tind1 = where(asttags eq strupcase(uphotmags[k]),ntind1)         ; U, G, R,..
       tind2 = where(phottags eq strupcase(uphotmags[k])+'MAG',ntind2)  ; UMAG, GMAG, RMAG, ...
       temp.(tind1) = phot[find].(tind2)
    endfor
    ast[sind] = temp
    ast.recovered = 0  ; nothing recovered to start
    ast[sind].recovered = 1

    ; Write out the ast file for this mock
    outmockfile = mockphotdir+'/'+mockphotbase+'_complete.fits'
    print,'  Writing AST catalog for this Field/Chip/Mock to ',outmockfile
    MWRFITS,ast,outmockfile,/create
    if file_test(outmockfile+'.gz') then file_delete,outmockfile+'.gz'
    SPAWN,['gzip',outmockfile],/noshell

    ; Combine the mock AST structure for this chip
    if n_elements(chipast) gt 0 then begin
      chipasttags = tag_names(chipast)
      bad = where(chipasttags ne asttags,nbad)
      if nbad gt 0 then begin
        printlog,logf,'  AST structures do NOT match between MOCKS'
        return 
      endif
      ; Combine
      nchipast = n_elements(chipast)
      chipastschema = chipast[0]
      struct_assign,{dum:''},chipastschema
      temp = chipast
      chipast = REPLICATE(chipastschema,nchipast+nsynth)
      chipast[0:nchipast-1] = temp
      chipast[nchipast:*] = ast
      undefine,temp
      nchipast = n_elements(chipast)
    endif else begin
      chipast = ast
      chipasttags = tag_names(chipast)
    endelse

    BOMB1:    
  Endfor  ; mock loop

  ; Write AST structure for this chip
  outchipfile = mockphotdir+'/'+shfield+'-'+refname+imager.separator+string(ichip,format='(i02)')+'_complete.fits'
  print,'  Writing AST catalog for this Field/Chip to ',outchipfile
  MWRFITS,chipast,outchipfile,/create
  if file_test(outchipfile+'.gz') then file_delete,outchipfile+'.gz'
  SPAWN,['gzip',outchipfile],/noshell

  ; Combine the chip AST structure for this field
  CHIPCOMBINE:
  if n_elements(bigast) gt 0 and n_elements(chipast) gt 0 then begin
    chipasttags = tag_names(chipast)
    bigasttags = tag_names(bigast)
    bad = where(bigasttags ne chipasttags,nbad)
    if nbad gt 0 then begin
      printlog,logf,'  AST structures do NOT match between CHIPS'
      return 
    endif
    ; Combine
    nbigast = n_elements(bigast)
    nchipast = n_elements(chipast)
    bigastschema = bigast[0]
    struct_assign,{dum:''},bigastschema
    temp = bigast
    bigast = REPLICATE(bigastschema,nbigast+nchipast)
    bigast[0:nbigast-1] = temp
    bigast[nbigast:*] = chipast
    undefine,temp
    nbigast = n_elements(bigast)
  endif else begin
    bigast = chipast
    bigasttags = tag_names(bigast)
   endelse
    
  ;stop
  
  BOMB2:
  
Endfor  ; chip loop

; Nothing to save
if n_elements(bigast) eq 0 then begin
  error = 'Nothing to save'
  printlog,logf,error
  return
endif

; Write AST structure for this chip
;fdir = mockphotdir
;; remove chipXX portion of directory if it's there
;if strmid(file_basename(dir),0,4) eq 'chip' then fdir=file_dirname(dir)
outastfile = maindir+'/'+globalfield+'_complete.fits'
print,'  Writing FINAL AST catalog to ',outastfile
MWRFITS,bigast,outastfile,/create
if file_test(outastfile+'.gz') then file_delete,outastfile+'.gz'
SPAWN,['gzip',outastfile],/noshell


; Make the figures
;-----------------
gdrecover = where(bigast.recovered eq 1,ngdrecover)
; Use g and i if they exist
if tag_exist(bigast,'INP_G') eq 1 and tag_exist(bigast,'INP_I') eq 1 then begin
  mag1ind = where(bigasttags eq 'INP_G',nmag1ind)
  mag2ind = where(bigasttags eq 'INP_I',nmag2ind) 
  xtit = 'g-i'
  ytit = 'g'
endif else begin
  ; Use 1st and 2nd mag columns
  magcols = where(strmid(bigasttags,0,4) eq 'INP_',nmagind)
  mag1ind = magcols[0]
  mag2ind = magcols[1]
  xtit = strmid(bigasttags[mag1ind],4)+'-'+strmid(bigasttags[mag2ind],4)
  ytit = strmid(bigasttags[mag1ind],4)
endelse
; Color and magnitude
mag = bigast.(mag1ind)
col = bigast.(mag1ind)-bigast.(mag2ind)
; Ranges and stepsizes
xr = [min(col),max(col)]
xr = [floor(xr[0]*20)/20., ceil(xr[1]*20)/20.]  ; round to nearest 0.05
yr = [min(mag),max(mag)]
yr = [floor(yr[0]*20)/20., ceil(yr[1]*20)/20.]  ; round to nearest 0.05
dx = 0.2
dy = 0.4
;xr = [-1,3.5]
;yr = [17.0,27.0]
; Bin the data in the CMD
undefine,dum
hess,col,mag,dum,imall,dx=dx,dy=dy,xr=xr,yr=yr,xarr=xarr,yarr=yarr,/noplot
if ngdrecover gt 0 then begin
  hess,col[gdrecover],mag[gdrecover],dum,imrec,dx=dx,dy=dy,xr=xr,yr=yr,xarr=xarr,yarr=yarr,/noplot
endif else imrec=imall*0

; Make some figures
if file_test(maindir+'plots/',/directory) eq 0 then file_mkdir,maindir+'/plots/'
setdisp
;loadcol,3
!p.font = 0
; Input ASTS
file = maindir+'/plots/'+globalfield+'_input'
ps_open,file,/color,thick=4,/encap
device,/inches,xsize=8.5,ysize=10.5
displayc,imall,xarr,yarr,/yflip,xtit=xtit,ytit=ytit,tit='Input ASTs for '+globalfield,charsize=1.3
ps_close
ps2png,file+'.eps',/eps
spawn,['epstopdf',file+'.eps'],/noshell
; Recovered
file = maindir+'/plots/'+globalfield+'_recovered'
ps_open,file,/color,thick=4,/encap
device,/inches,xsize=8.5,ysize=10.5
displayc,imrec,xarr,yarr,/yflip,xtit=xtit,ytit=ytit,tit='Recovered ASTs for '+globalfield,charsize=1.3
ps_close
ps2png,file+'.eps',/eps
spawn,['epstopdf',file+'.eps'],/noshell
; Completeness
file = maindir+'/plots/'+globalfield+'_completeness'
ps_open,file,/color,thick=4,/encap
device,/inches,xsize=8.5,ysize=10.5
displayc,float(imrec)/(imall>1),xarr,yarr,/yflip,xtit=xtit,ytit=ytit,tit='Completeness for '+globalfield,charsize=1.3
ps_close
ps2png,file+'.eps',/eps
spawn,['epstopdf',file+'.eps'],/noshell
; Combine the figures
pdffiles = maindir+'/plots/'+globalfield+'_'+['input','recovered','completeness']+'.pdf'
spawn,'gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile='+maindir+'/plots/'+globalfield+'_complete.pdf '+strjoin(pdffiles,' ')
printlog,logf,'Completeness plots are in '+maindir+'/plots/'+globalfield+'_complete.pdf'

;stop

end
