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
;  bigsynthfile  File with all of the calibrated photometry for the
;                   artificial stars.
;  =imager       The imager information structure.
;  =logfile      A logfile to print output to.
;  =maindir      The main directory for this PHOTRED run.
;
; OUTPUTS:
;  The structure of artificial stars and information on whether
;  they were recovered and with what magnitudes is saved.
;  =error     The error message if there was one.
;
; USAGE:
;  IDL>completeness,photfiles,bigsynthfile,imager=imager,maindir=maindir
; 
; By D.Nidever  based on getcomplete.pro   July 2017
;-

pro completeness,photfiles,bigsynthfile,imager=imager,logfile=logfile,$
                 maindir=maindir,error=error

undefine,error
  
; Not enough inputs
nphotfiles = n_elements(photfiles)
if nphotfiles eq 0 or n_elements(bigsynthfile) eq 0 or n_elements(imager) eq 0 or n_elements(maindir) eq 0 then begin
  error = 'Not enough inputs'
  print,'Syntax - completeness,photfiles,bigsynthfile,imager=imager,maindir=maindir,logfile=logfile,error=error'
  return
endif

; Logfile
if keyword_set(logfile) then logf=logfile else logf=-1

; Load the file with all the artificial star information
;   and the calibrated photometry
if (file_info(bigsynthfile)).exists eq 0 then begin
  error = 'BIGSYNTHFILE='+bigsynthfile+' NOT FOUND'
  printlog,logfile,error
  return
endif
if first_el(strsplit(bigsynthfile,'.',/extract),/last) eq 'fits' then $
  bigsynthstr = MRDFITS(bigsynthfile,1,/silent) else $
  bigsynthstr = IMPORTASCII(bigsynthfile,/header,/silent)

; Get the "short" field name, e.g. F4
;  The names should look like this: F4M6-00478573_22.phot
allfieldmockname = reform( (strsplitter(file_basename(photfiles),'-',/extract))[0,*] )
allfieldname = reform( (strsplitter(allfieldmockname,'M',/extract))[0,*] )
uifieldname = uniq(allfieldname,sort(allfieldname))
if n_elements(uifieldname) gt 1 then begin
  error = 'ERROR: files for more than one field input'
  printlog,logfile,error
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
  printlog,logfile,error
  return
endif
globalfield = lnames[nameind[0]]

printlog,logfile,'Getting completeness function for ',shfield,' = ',globalfield
printlog,logfile,'-------------------------------------------------------------'


; Loop over the separate chips
;-----------------------------
undefine,bigast
For i=0,nchips-1 do begin

  ichip = uchips[i]
  chipind = where(allchips eq ichip,nmocks)
  printlog,logfile,'  '+strtrim(i+1,2)+'/'+strtrim(nchips,2)+' CHIP='+strtrim(ichip,2)+' - '+strtrim(nmocks,2)+' mocks'

  ; Loop over the mocks for this chip
  ;----------------------------------
  undefine,chipast
  For j=0,nmocks-1 do begin
    mockphotfile = photfiles[chipind[j]]
    mockphotdir = file_dirname(mockphotfile)
    mockphotbase = file_basename(mockphotfile)
    ; Get the mock number for the filename, e.g. F2M5- gives mocknum=5
    fieldmockname = first_el(strsplit(mockphotbase,'-',/extract))
    mocknum = long(first_el(strsplit(fieldmockname,'M',/extract),/last))
    printlog,logfile,'  MOCK='+strtrim(mocknum,2)

    ; Get the reference exposure number
    ;  00379732 in F4M4-00379732_22.phot
    refname = (strsplit(mockphotbase,'-'+imager.separator,/extract))[1]
    
    ; Check phot file header to see if ALLFRAME was run (see if PROB exists)
    READLINE,mockphotfile,phothead,nlineread=1
    photcols = strtrim(strsplit(phothead[0],' ',/extract),2)
    probind = where(strupcase(photcols) eq 'PROB',nprobind)
    if nprobind gt 0 then allframe=1 else allframe=0
    ; Load the PHOT file
    phot = IMPORTASCII(mockphotfile,/header,/silent)
    nphot = n_elements(phot)
    
    ; -- Add average magnitudes to PHOT file --
    ; We need average magnitues per band for reach object
    phottags = tag_names(phot)
    photavgmag = where(stregex(phottags,'MAG$',/boolean) eq 1,nphotavgmag)
    if nphotavgmag eq 0 then begin
      ; Find all the magnitude columns
      photmagind = where(stregex(photcols,'MAG',/boolean) eq 1 and $
                         stregex(photcols,'^I_',/boolean) eq 0 and $
                         stregex(photcols,'ERR',/boolean) eq 0,nphotmagind)
      if nphotmagind eq 0 then begin
        printlog,logfile,'NO photometric magnitudes found in '+mockphotfile
        goto,BOMB1 
      endif
      allphotmags = photcols[photmagind]  ; iMAG3, zMAG
      allphotmagsband = reform((strsplitter(allphotmags,'MAG',/extract))[0,*])  ; just the band name
      ; Unique magnitudes
      uiphotmagsband = uniq(allphotmagsband,sort(allphotmagsband))
      uphotmagsband = allphotmags[uiphotmagsband]
      nphotmagsband = n_elements(photmagsband)
      ; Add the columns to PHOT
      photschema = phot[0]
      STRUCT_ASSIGN,photschema,{dum:0}  ; blank it out
      for k=0,nphotmagsband-1 do $
        photschema=CREATE_STRUCT(photschema,uphotmagsband[k]+'MAG',0.0,uphotmagsband[k]+'ERR',0.0)
      temp = phot & undefine,phot
      phot = REPLICATE(photschema,nphot)  ; put in new structure format
      STRUCT_ASSIGN,temp,phot,/nozero
      phottag = tag_names(phot)  ; new phot tags
      
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
          
        ; Average the mags
        AVERAGEMAG,tempmag,temperr,newmag,newerr,/robust
        magind = where(phottags eq uphotmagsband[k]+'MAG',nmagind)
        errind = where(phottags eq uphotmagsband[k]+'ERR',nerrind)
        phot.(magind) = newmag
        phot.(errind) = newerr
      endfor  ; unique band loop      
    endif  ; adding average photometry per band      

    
    ; Load the final, recovered photometry file, but use the mag/raw
    ;  instrumental photometry.  Use this to match up recovered to
    ;  original objects.
    ;  Use .raw or .mag file depending if allframe was used.
    if allframe eq 1 then magext='.mag' else magext='.raw'
    finalfile = repstr(mockphotfile,'.phot','.mag')
    if (file_info(finalfile)).exists eq 0 then begin
      printlog,logfile,finalfile+' NOT FOUND'
      goto,BOMB1
    endif
    if allframe eq 1 then LOADMAG,finalfile,final else LOADRAW,finalfile,final
    nfinal = n_elements(final)
    printlog,logfile,'  NFINAL='+strtrim(nfinal,2)
    
    ; Load the original "real" star data file
    ;  Use .raw or .mag file depending if allframe was used.
    origfile = mockphotdir+'/'+shfield+'-'+refname+imager.separator+string(ichip,format='(i02)')+magext
    if (file_info(origfile)).exists eq 0 then begin
      printlog,logfile,origfile+' NOT FOUND'
      goto,BOMB1
    endif
    if allframe eq 1 then LOADMAG,origfile,orig else LOADRAW,origfile,orig
    norig = n_elements(orig)
    printlog,logfile,'  NORIG='+strtrim(norig,2)
    
    ; Load the artificial star data file, with information on
    ;   the injected stars
    synthfile = mockphotdir+'/'+fieldmockname+'-'+refname+'-add'+imager.separator+string(ichip,format='(i02)')+magext
    if (file_info(synthfile)).exists eq 0 then begin
      printlog,logfile,synthfile+' NOT FOUND'
      goto,BOMB1
    endif
    if allframe eq 1 then LOADMAG,synthfile,synth else LOADRAW,synthfile,synth
    nsynth = n_elements(synth)
    printlog,logfile,'  NSYNTH='+strtrim(nsynth,2)

    ; Do the CROSS-MATCHING between the three lists (final/orig/synth)
    ;-----------------------------------------------------------------
    ; If two stars land right on top of each other and only one is recovered then use the magnitude
    ; to match it up to the correct star (real or artificial)
    ; We don't want to think we recovered an artificial star when it's
    ; actually a real star (likely bright one).
    ; We are using the instrumental photomery (.mag/.raw) to do this matching.
    COMPLETENESS_CROSSMATCH,final,orig,synth,find,sind,nmatch=nmatch,logfile=logfile,error=xmatcherror
    if n_elements(xmatcherror) gt 0 then goto,BOMB1
    
    ;; Crossmatch final to orig and final to synth
    ;SRCMATCH,final.x,final.y,orig.x,orig.y,2,oind1,oind2,count=nomatch
    ;SRCMATCH,final.x,final.y,synth.x,synth.y,2,aind1,aind2,count=nastmatch
    ;
    ;; Deal with duplicates
    ;;  all we care about is real objects falsely
    ;;  identified as ASTs, so bad matches in AIND
    ;finaldblind = doubles([oind1,aind1],count=nfinaldbl)
    ;if nfinaldbl gt 0 then begin
    ;  printlog,logfile,'  Resolving '+strtrim(nfinaldbl,2)+' duplicates'
    ;  finaldbl = ([oind1,aind1])(finaldblind)
    ;  flag = lonarr(nfinaldbl)  ; 0-real, 1-ast
    ;  origdbl = lonarr(nfinaldbl)
    ;  synthdbl = lonarr(nfinaldbl)
    ;  ; Get magnitude column indices for the three structures
    ;  finaltags = tag_names(final)
    ;  finalmagind = where(stregex(finaltags,'^MAG',/boolean) eq 1,nfinalmagind)
    ;  origtags = tag_names(orig)
    ;  origmagind = where(stregex(origtags,'^MAG',/boolean) eq 1,norigmagind)
    ;  synthtags = tag_names(synth)
    ;  synthmagind = where(stregex(synthtags,'^MAG',/boolean) eq 1,nsynthmagind)
    ;  if nfinalmagind ne norigmagind or nfinalmagind ne nsynthmagind then begin
    ;    printlog,logfile,'  FINAL/ORIG/SYNTH photometry files have different magnitude columns'
    ;    goto,BOMB1
    ;  endif
    ;  ; Loop over duplicates
    ;  for k=0,nfinaldbl-1 do begin
    ;    finaldbl1 = finaldbl[k]  ; the FINAL index
    ;    final1 = final[finaldbl1]
    ;    finalmag = fltarr(nfinalmagind)
    ;    for l=0,nfinalmagind-1 do finalmag[l]=final1.(finalmagind[l])
    ;
    ;    ; ORIG
    ;    MATCH,oind1,finaldbl1,ind1,/sort
    ;    origdbl[k] = ind1  ; index into OIND1/2
    ;    orig1 = orig[oind2[ind1]]
    ;    odist = sqrt( (final1.x-orig1.x)^2 + (final1.y-orig1.y)^2 )
    ;    origmag = fltarr(norigmagind)
    ;    for l=0,norigmagind-1 do origmag[l]=orig1.(origmagind[l])
    ;    gdorig = where(finalmag lt 50 and origmag lt 50,ngdorig)
    ;    if ngdorig gt 0 then omagdiff = mean(abs(origmag[gdorig]-finalmag[gdorig])) else omagdiff=99.99
    ;    ofinaldiff = sqrt(odist^2 + omagdiff^2)
    ;
    ;    ; SYNTH
    ;    MATCH,aind1,finaldbl1,ind2,/sort
    ;    synthdbl[k] = ind2  ; index into AIND1/2
    ;    synth1 = synth[aind2[ind2]]
    ;    adist = sqrt( (final1.x-synth1.x)^2 + (final1.y-synth1.y)^2 )
    ;    synthmag = fltarr(nsynthmagind)
    ;    for l=0,nsynthmagind-1 do synthmag[l]=synth1.(synthmagind[l])
    ;    gdsynth = where(finalmag lt 50 and synthmag lt 50,ngdsynth)
    ;    if ngdsynth gt 0 then amagdiff = mean(abs(synthmag[gdsynth]-finalmag[gdsynth])) else amagdiff=99.99
    ;    afinaldiff = sqrt(adist^2 + amagdiff^2)
    ;
    ;    ; Which one is the match
    ;    if ofinaldiff lt afinaldiff then begin
    ;      com='  REAL'
    ;      flag[k] = 0
    ;    endif else begin
    ;      com='  AST'
    ;      flag[k] = 1
    ;    endelse
    ;    if keyword_set(verbose) then $
    ;      printlog,logfile,'  '+strtrim(i+1,2),odist,omagdiff,ofinaldiff,' ',adist,amagdiff,afinaldiff,com
    ;    ;if amagdiff lt omagdiff and adist gt odist then stop
    ;  endfor  ; duplicates loop
    ;
    ;  ; Remove ASTs from the "ORIG" list
    ;  bdomatch = where(flag eq 1,nbdomatch,ncomp=ngdomatch)  
    ;  print,strtrim(nbdomatch,2),' are ASTs and ',strtrim(ngdomatch,2),' are REAL sources'
    ;  if nbdomatch gt 0 then begin
    ;    if nbdomatch lt nomatch then begin
    ;      bdorigdbl = origdbl[bdomatch]
    ;      REMOVE,bdorigdbl,oind1,oind2
    ;      nomatch = n_elements(oind1)
    ;    endif else begin
    ;      undefine,oind1,oind2
    ;      nomatch = 0
    ;    endelse
    ;  endif
    ;endif ; duplicates
    ;
    ;; "Prune" the real sources from the list of recovered sources
    ;left = phot  ; use the calibrated photometry, same sources/order as FINAL
    ;if nomatch gt 0 then remove,oind1,left
    ;; Now rematch SYNTH to the leftover sources
    ;SRCMATCH,left.x,left.y,synth.x,synth.y,2,aind1,aind2,count=nastmatch
    ;print,strtrim(nastmatch,2),' ASTs recovered'

    
    ; Put all AST information into one structure
    ;-------------------------------------------

    ; --- Create the structure schema ---
    astschema = {astid:'',photid:'',recovered:-1,field:'',chip:0L,mock:0L,$
              inp_x:0.0d0,inp_y:0.0d0}
    ; Add columns for input synthetic calibrated photometry
    bigsynthtags = tag_names(bigsynthstr)
    bigsynthphotcolind = where(bigsynthtags ne 'ID',nbigsynthphotcolind)
    bigsynthphotcols = bigsynthtags[bigsynthphotcolind]
    for k=0,nbigsynthphotcolind-1 do astschema=CREATE_STRUCT(astschema,'INP_'+bigsynthphotcols[k],0.0)
    ; Add recovered information
    astschema = CREATE_STRUCT(astschema,'ID','','X',0.0d0,'Y',0.d0,'NDET',0L)
    ; Add columms from PHOT file
    ;   get all of the unique filters, anything with MAG and not I_ and
    ;   not with ERR
    photmagind = where(stregex(photcols,'MAG',/boolean) eq 1 and $
                       stregex(photcols,'^I_',/boolean) eq 0 and $
                       stregex(photcols,'ERR',/boolean) eq 0,nphotmagind)
    if nphotmagind eq 0 then begin
      printlog,logfile,'NO photometric magnitudes found in '+mockphotfile
      goto,BOMB1 
    endif
    allphotmags = photcols[photmagind]  ; iMAG3, zMAG
    allphotmags = reform((strsplitter(allphotmags,'MAG',/extract))[0,*])  ; just the mag name
    uiphotmags = uniq(allphotmags,sort(allphotmags))
    uphotmags = allphotmags[uiphotmags]
    nphotmags = n_elements(photmags)
    ;  Add three columns for each filter: magnitude, error, Ndetections
    for k=0,nphotmags-1 do $
       astschema=CREATE_STRUCT(astschema,photmags[k],0.0,photmags[k]+'ERR',0.0,'NDET'+photmags[k],0L)
    ; Add final morphology columns
    astschema = CREATE_STRUCT(astschema,'chi',0.0,'sharp',0.0,'flag',0L,'prob',0.0)
    asttags = tag_names(astschema)
    
    ; --- Make the structure and copy over the data ---
    ast = replicate(astschema,nsynth)
    ; Copy over the input synth data
    ;   ASTID is a unique name across all mocks and chips
    ast.astid = shfield+'.'+strtrim(ichip,2)+'.'+strtrim(mocknum,2)+'.'+strtrim(synth.id,2)
    ast.photid = synth.id
    ast.field = shfield
    ast.chip = ichip
    ast.mock = mocknum
    ast.inp_x = synth.x
    ast.inp_y = synth.y
    ; Get the calibrated photometry from BIGSYNTHSTR
    ;  match2.pro deals with duplicates better
    ;  The IDs in the add.mag file (SYNTH) and BIGSYNTHSTR should match
    MATCH2,bigsynthstr.id,synth.id,ind1,ind2
    dum = where(ind2 gt -1,nbigmatch)
    if nbigmatch ne nsynth then begin
      printlog,logfile,'Not all synth elements found matches in bigsynthstr'
      goto,BOMB1
    endif
    ; loop over bigsynth columns and copy
    for k=0,nbigsynthcolind-1 do begin
      asttagind = where(asttags eq 'INP_'+bigsynthphotcols[k],nasttagind)
      bigsynttagind = where(bigsynthtags eq bigsynthcphotcols[k],nbigsynthtagind)
      ast.(asttagind) = bigsynthstr[ind1].(bigsynthtagind)
    endfor
    ; Copy over the recovered values using STRUCT_ASSIGN
    ;  LEFT was created from the PHOT structure
    ;temp = ast[aind2]
    ;STRUCT_ASSIGN,left[aind1],temp,/nozero
    ;ast[aind2] = temp
    ;ast[aind2].recovered = 1
    temp = ast[sind]
    STRUCT_ASSIGN,phot[find],temp,/nozero
    ast[sind] = temp
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
        printlog,logfile,'  AST structures do NOT match between MOCKS'
        return 
      endif
      ; Combine
      nchipast = n_elements(chipast)
      temp = chipast
      chipast = REPLICATE(chipasttags,nchipast+nsynth)
      chipast[0:nchipast-1] = temp
      chipast[nchipast:*] = ast
      undefine,temp
      nchipast = n_elements(chipast)
    endif else chipast=ast
    
    stop

    BOMB1:
    
  Endfor  ; mock loop

  ; Write AST structure for this chip
  outchipfile = mockphotdir+'/'+shfield+'-'+refname+imager.separator+string(ichip,format='(i02)')+'_complete.fits'
  print,'  Writing AST catalog for this Field/Chip to ',outchipfile
  MWRFITS,chipast,outchipfile,/create
  if file_test(outchipfile+'.gz') then file_delete,outchipfile+'.gz'
  SPAWN,['gzip',outchipfile],/noshell

  ; Combine the chip AST structure for this field
  if n_elements(bigast) gt 0 and n_elements(chipast) gt 0 then begin
    bigasttags = tag_names(bigast)
    bad = where(bigasttags ne chipasttags,nbad)
    if nbad gt 0 then begin
      printlog,logfile,'  AST structures do NOT match between CHIPS'
      return 
    endif
    ; Combine
    nbigast = n_elements(bigast)
    nchipast = n_elements(chipast)
    temp = bigast
    bigast = REPLICATE(bigasttags,nbigast+nchipast)
    bigast[0:nbigast-1] = temp
    bigast[nbigast:*] = chipast
    undefine,temp
    nbigast = n_elements(bigast)
  endif else bigast=chipast
  
  stop
  
  BOMB2:
  
Endfor  ; chip loop

; Nothing to save
if n_elements(bigast) eq 0 then begin
  error = 'Nothing to save'
  printlog,logfile,error
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


; Figure out the completeness
gdrecover = where(bigast.recovered eq 1,ngdrecover)
dx = 0.2
dy = 0.4
xr = [-1,3.5]
yr = [17.0,27.0]
hess,bigast.inp_g-bigast.inp_i,bigast.inp_g,dum,imall,dx=dx,dy=dy,xr=xr,yr=yr,xarr=xarr,yarr=yarr,/noplot
hess,bigast[gdrecover].inp_g-bigast[gdrecover].inp_i,bigast[gdrecover].inp_g,dum,imrec,dx=dx,dy=dy,xr=xr,yr=yr,xarr=xarr,yarr=yarr,/noplot


; Make some figures
if file_test(maindir+'plots/',/directory) eq 0 then file_mkdir,maindir+'plots/'
setdisp
;loadcol,3
!p.font = 0
; Input ASTS
file = maindir+'plots/'+globalfield+'_input'
ps_open,file,/color,thick=4,/encap
device,/inches,xsize=8.5,ysize=10.5
displayc,imall,xarr,yarr,/yflip,xtit='g-i',ytit='g',tit='Input ASTs for '+globalfield,charsize=1.3
ps_close
ps2png,file+'.eps',/eps
spawn,['epstopdf',file+'.eps'],/noshell
; Recovered
file = maindir+'plots/'+globalfield+'_recovered'
ps_open,file,/color,thick=4,/encap
device,/inches,xsize=8.5,ysize=10.5
displayc,imrec,xarr,yarr,/yflip,xtit='g-i',ytit='g',tit='Recovered ASTs for '+globalfield,charsize=1.3
ps_close
ps2png,file+'.eps',/eps
spawn,['epstopdf',file+'.eps'],/noshell
; Completeness
file = maindir+'plots/'+globalfield+'_completeness'
ps_open,file,/color,thick=4,/encap
device,/inches,xsize=8.5,ysize=10.5
displayc,float(imrec)/(imall>1),xarr,yarr,/yflip,xtit='g-i',ytit='g',tit='Completeness for '+globalfield,charsize=1.3
ps_close
ps2png,file+'.eps',/eps
spawn,['epstopdf',file+'.eps'],/noshell
; Combine the figures
pdffiles = maindir+'/plots/'+globalfield+'_'+['input','recovered','completeness']+'.pdf'
spawn,'gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile='+maindir+'/plots/'+globalfield+'_complete.pdf '+strjoin(pdffiles,' ')

stop

end
