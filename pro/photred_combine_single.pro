;+
;
; PHOTRED_COMBINE_SINGLE
;
; This program combines the data from all the chips
; See comb_mscphot.pro
;
; INPUTS:
;  files         Files to combine.
;  imagerfile    Filename of imager information
;  fieldname     Field name
;  =catformat    Type of output catalog format to use.
;  =telescope    Telescope name.
;  =instrument   Instrument name.
;  /force        Force combination even if not all amps are there.
;  /posonly    
;  /mchusetiles  Use tiles.
;  /redo         Redo files that were already done.
;  /stp          Stop at the end of the program.
;
; OUTPUTS:
;  The combined photometry files for each pointing
;
; By D.Nidever  Mar 2008
;-

pro photred_combine_single,files,imagerfile,fieldname,telescope=telescope,instrument=instrument,$
                           force=force,posonly=posonly,catformat=catformat,$
                           mchusetiles=mchusetiles,sepchipdir=sepchipdir,$
                           redo=redo,stp=stp

  CD,current=curdir

  if n_elements(catformat) eq 0 then catformat='FITS'

  ; LOAD THE "imagers" FILE
  ;----------------------------
  print,'Loading imager information'
  ; The columns need to be: Telescope, Instrument, Naps, separator
  imagers_fieldnames = ['telescope','instrument','observatory','namps','separator']
  imagers_fieldtpes = [7,7,7,3,7]
  imagers = IMPORTASCII(imagerfile,fieldnames=imagers_fieldnames,$
                        fieldtypes=imagers_fieldtypes,comment='#')
  imagers.telescope = strupcase(strtrim(imagers.telescope,2))
  imagers.instrument = strupcase(strtrim(imagers.instrument,2))
  imagers.observatory = strupcase(strtrim(imagers.observatory,2))
  singleind = where(imagers.namps eq 1,nsingle)
  if nsingle gt 0 then imagers[singleind].separator = ''
  if (n_tags(imagers) eq 0) then begin
    print,'NO imagers in '+imagerfile
    return
  endif

  ; What IMAGER are we using??
  ;---------------------------
  ind_imager = where(imagers.telescope eq strupcase(telescope) and $
                     imagers.instrument eq strupcase(instrument),nind_imager)
  if nind_imager eq 0 then begin
    print,'TELESCOPE='+telescope+' INSTRUMENT='+instrument+' NOT FOUND in >>imagers<< file'
    return
  endif
  thisimager = imagers[ind_imager[0]]
  ; print out imager info
  print,''
  print,'USING IMAGER:'
  print,'Telescope = '+thisimager.telescope
  print,'Instrument = '+thisimager.instrument
  print,'Namps = '+strtrim(thisimager.namps,2)
  print,"Separator = '"+thisimager.separator+"'"
  print,''


  print,''
  print,'================================================'
  print,'Combining photometry files for '+fieldname
  print,'================================================'
  print,''
  print,systime(0)

  fieldlines = files
  nfieldlines = n_elements(fieldlines)
  base = FILE_BASENAME(fieldlines,'.phot')
  basedir = FILE_DIRNAME(fieldlines[0])
  shortfield = (strsplit(base[0],'-',/extract))[0]
  print,strtrim(nfieldlines,2)+' File(s) for FIELD='+shortfield
  
  ;; Not enough files for this imager
  if (thisimager.namps gt 1) and (nfieldlines ne thisimager.namps) and not keyword_set(force) then begin
    print,'Only '+strtrim(nfieldlines,2)+' for '+fieldname+'.  Need '+strtrim(thisimager.namps,2)+' for '+$
                     thisimager.telescope+'+'+thisimager.instrument
    return
  endif

  print,'COMBINING '+strtrim(nfieldlines,2)+' files for '+fieldname

  ;; Getting the main name
  firstname = FILE_BASENAME(fieldlines[0],'.phot')
  basedir = FILE_DIRNAME(fieldlines[0])
  if keyword_set(mchusetiles) then basedir=FILE_DIRNAME(basedir)  ;; F1/F1-T10/F1-00507801+T10.phot
  if not keyword_set(mchusetiles) and keyword_set(sepchipdir) then basedir=FILE_DIRNAME(basedir)  ;; F1/chip01/
  ;; Getting basename
  basename = firstname   ;; namps=1 and not using tiles
  if keyword_set(mchusetiles) then basename = (strsplit(firstname,'\'+tilesep+'T',/extract))[0]
  if not keyword_set(mchusetiles) and (thisimager.namps gt 1) then begin
    ending = first_el(strsplit(firstname,thisimager.separator,/extract),/last)
    endlen = strlen(ending)
    len = strlen(firstname)
    basename = strmid(firstname,0,len-endlen-1)
  endif

  ;; Get file information and find unique exposures
  undefine,allfiles
  for j=0,nfieldlines-1 do begin
    dir1 = file_dirname(fieldlines[j])
    base1 = file_basename(fieldlines[j],'.phot')
    LOADMCH,dir1+'/'+base1+'.mch',indivfiles
    push,allfiles,dir1+'/'+file_basename(indivfiles,'.als')+'.fits'
  endfor
  bd = where(file_test(allfiles) eq 0,nbd)
  if nbd gt 0 then allfiles[bd]+='.fz'
  ;; Get exposure names
  ;; single-amp F1-12340044.als
  ;; multi-amp  F1-12340044_10.als
  arr1 = strsplitter(file_basename(allfiles),'-',/extract)
  allexpnum = reform(arr1[1,*])         ;; single-amp
  if thisimager.namps gt 1 then begin   ;; multi-amp
    arr2 = strsplitter(reform(arr1[1,*]),thisimager.separator,/extract)
    allexpnum = reform(arr2[0,*])
  endif
  uiexp = uniq(allexpnum,sort(allexpnum))
  uexpnum = allexpnum[uiexp]
  nexp = n_elements(uexpnum)
  ;; Get exposure information
  expstr = replicate({expnum:'',filter:'',dateobs:'',mjd:0.0d0,exptime:0.0,cenra:0.0d0,cendec:0.0d0},nexp)
  PHOTRED_GATHERFILEINFO,allfiles[uiexp],filestr
  struct_assign,filestr,expstr
  expstr.expnum = uexpnum
  for j=0,nexp-1 do expstr[j].mjd = date2jd(expstr[j].dateobs,/mjd)
  si = sort(expstr.mjd)  ; put in chronological order
  expstr = expstr[si]
  print,strtrim(nexp,2)+' exposures'

  ;; Check that all of the files have the same reference exposure
  fbase = file_basename(fieldlines,'.phot')
  fbase = reform((strsplitter(fbase,'-',/extract))[1,*])
  fbase = reform((strsplitter(fbase,'_',/extract))[0,*])
  fui = uniq(fbase,sort(fbase))
  if n_elements(fui) gt 1 then begin
    ufbase = fbase[fui]
    findex = create_index(fbase)
    print,'More than one REFERENCE exposure used for this field: ',strjoin(findex.value+' ('+strtrim(findex.num,2)+' files)',', ')
    bestind = first_el(maxloc(findex.num))
    ind = findex.index[findex.lo[bestind]:findex.hi[bestind]]
    fieldlines = fieldlines[ind]
    nfieldlines = n_elements(fieldlines)
    print,'ASSUMING ',findex.value[bestind],' is the correct REFERENCE exposure - ',strtrim(nfieldlines,2),' PHOT files left'
  endif  

  ;-------------------------------------------------
  ; LOOP through all the PHOT files for this field
  ;-------------------------------------------------
  ncombined = 0
  undefine,str,all
  undefine,fieldnames0,fieldtypes0
  For j=0,nfieldlines-1 do begin
    file = fieldlines[j]
    filebase = FILE_BASENAME(file,'.phot')
    filedir = FILE_DIRNAME(file)

    ; Check that the PHOT file exists
    test = FILE_TEST(file)
    if test eq 1 then photlines = FILE_LINES(file) else photlines=0
    if (photlines eq 0) then begin
      PUSH,failurelist,file
      if test eq 0 then print,file+' NOT FOUND'
      if test eq 1 and nlines eq 0 then print,file+' HAS 0 LINES'
      goto,BOMB2
    endif

    ;; Load the PHOT file
    str = PHOTRED_READFILE(file,count=nstr)

    ;; Get the column names and types
    fieldnames = tag_names(str)
    fieldtypes = lonarr(n_elements(fieldnames))
    for k=0,n_elements(fieldnames)-1 do fieldtypes[k]=size(str[0].(k),/type)
    ;; ID must be a string
    idind = where(fieldnames eq 'ID',nidind)
    if fieldtypes[idind] ne 7 then begin
      fieldtypes[idind] = 7
      str0 = str[0] & struct_assign,{dum:''},str0
      schema = create_struct(fieldnames[0],fix(str0.(0),type=fieldtypes[0]))
      for k=1,n_elements(fieldnames)-1 do schema=create_struct(schema,fieldnames[k],fix(str0.(k),type=fieldtypes[k]))
      str_orig = str
      str = replicate(schema,nstr)
      struct_assign,str_orig,str
      str.id = strtrim(str.id,2)
      undefine,str_orig
    endif

    ;; Print out the file info
    print,''
    print,'ADDING '+filebase+' Nstars='+strtrim(nstr,2)
    print,''

    ;; Adding EXT tag for amp number
    tags = TAG_NAMES(str)
    if (thisimager.namps gt 1) or keyword_set(mchusetiles) then begin
      extgd = where(tags eq 'EXT',nextgd)
      if nextgd eq 0 then ADD_TAG,str,'EXT',0,str
    endif

    ;; Single-amp and NOT using tiles
    ;;-------------------------------
    If (not keyword_set(mchusetiles)) and (thisimager.namps eq 1) then begin
      ;; Updating the IDs
      ;; FIELD.IDNUMBER, i.e. 190L182a.17366
      ;;---------------------------------------------
      id2 = fieldname+'.'+strtrim(str.id,2)
      str.id = id2
    Endif
       
    ;; Multi-amps and NOT using tiles
    ;;-------------------------------
    If (not keyword_set(mchusetiles)) and (thisimager.namps gt 1) then begin
      ;; FORCE the format to be the same for ALL chip files, otherwise we sometimes
      ;;  get "type mismatch" errors with double/floats
      if j eq 0 then begin
        fieldnames0 = fieldnames
        fieldtypes0 = fieldtypes
        file0 = file
      endif
      ;; Check that the fieldnames are the same
      ;; Just a warning now.  Formats are fixed below with
      ;; photred_combine_reformatphot.pro even if there are missing files
      if n_elements(fieldnames0) ne n_elements(fieldnames) then begin
        print,''
        print,'Warning: '+file+' format does NOT agree with '+file0
      endif
      if total(strcmp(fieldnames0,fieldnames)) ne n_elements(fieldnames0) then begin
        print,''
        print,'Warning: '+file+' format does NOT agree with '+file0
      endif      
      ; Getting the amplifier number (extension)
      ;-----------------------------------------
      ext = photred_getchipnum(filebase,thisimager)
      str.ext = ext
      ; Updating the IDs
      ; FIELD_EXT.IDNUMBER, i.e. 190L182a_5.17366
      ;---------------------------------------------
      id2 = fieldname+'_'+strtrim(ext,2)+'.'+strtrim(str.id,2)
      str.id = id2
    Endif
      
    ;; Using TILES, need to reformat photometry structure
    If keyword_set(mchusetiles) then begin
      ;; The chip/amp information is not in the name
      ;; each image could have a different chip
      ;; F1-00507801+T2, no chip extension
      ;; Use the TILE Number as the EXT
      ext = first_el(strsplit(filebase,tilesep+'T',/extract),/last)
      str.ext = ext
      ; Updating the IDs
      ; FIELD_TILE.IDNUMBER, i.e. 190L182a_5.17366
      ;---------------------------------------------
      id2 = fieldname+'_'+strtrim(ext,2)+'.'+strtrim(str.id,2)
      str.id = id2
    Endif
      
    ;; Reformat the photometry structure
    ;; Put the photometry in exposure columns
    ;;---------------------------------------
    ;; Get information on the individual images
    LOADMCH,filedir+'/'+filebase+'.mch',alsfiles
    fitsfiles = filedir+'/'+file_basename(alsfiles,'.als')+'.fits'
    bdfits = where(file_test(fitsfiles) eq 0,nbdfits)
    if nbdfits gt 0 then fitsfiles[bdfits]+='.fz'
    PHOTRED_GATHERFILEINFO,fitsfiles,filestr
    ;; Get exposure names
    add_tag,filestr,'expnum','',filestr
    arr1 = strsplitter(file_basename(filestr.file),'-',/extract)
    allexpnum = reform(arr1[1,*])         ;; single-amp
    if thisimager.namps gt 1 then begin   ;; multi-amp
      arr2 = strsplitter(reform(arr1[1,*]),thisimager.separator,/extract)
      allexpnum = reform(arr2[0,*])
    endif
    filestr.expnum = allexpnum
    ;; Reformat the photometry structure
    str_orig = str & undefine,str

    PHOTRED_COMBINE_REFORMATPHOT,str_orig,filestr,expstr,str

    ; Check that the structure has RA/DEC
    ;------------------------------------
    gdra = where(tags eq 'RA',ngdra)
    gddec = where(tags eq 'DEC',ngddec)
    if (ngdra eq 0) or (ngddec eq 0) then begin
      print,file+' DOES NOT HAVE RA/DEC FIELDS'
      goto,BOMB2
    endif

    ; COMBINING
    ;-----------
    ; Adding together with PHOT_OVERLAP.PRO to deal with overlaps
    ; at the edges of chips/amplifiers
    ; This uses RA/DEC so if these aren't in the structures it won't
    ; "combine" stars properly.
    if (n_elements(all) gt 0) then begin
      PHOT_OVERLAP,all,str,outstr,silent=silent,posonly=cmbposonly
      noutstr = n_elements(outstr)
      ;; There was a problem
      if (noutstr eq 0) then begin
        print,'Problem in PHOT_OVERLAP for '+file
        goto,BOMB2
      endif        
      all = outstr
      ncombined++

    ; First time, initialize the "all" structure
    endif else begin
      all = str
      ncombined++
    endelse

    BOMB2:
  Endfor  ; loop through the individual PHOT files

  ;; At least one file failed
  ;; Field and files are successful only if *ALL* succeeded
  ;;--------------------------------------------------------
  if (ncombined ne nfieldlines) then begin
    print,fieldname,' FAILED because ',strtrim(ngd-ncombined,2),' Files FAILED'
    return
  endif

  ;; Are there any stars?
  ;;----------------------
  nall = n_elements(all)
  if (nall eq 0) then begin
    print,'NO STARS for FIELD='+fieldname
    return
  endif
  
  ;; We were successful!
  ;;-----------------------
  print,''
  print,'Combined '+strtrim(ncombined,2)+' files'
  print,'Nstars = '+strtrim(nall,2)+' for FIELD='+fieldname

  ;; Output the data
  ;;----------------
  outname = basedir+'/'+basename+'.cmb'
  print,'OUTPUTTING data to '+outname
  if catformat eq 'ASCII' then begin
    PRINTSTR,all,outname
    PRINTSTR,expstr,outname+'.meta'
  endif else begin
    MWRFITS,all,outname,/create
    MWRFITS,expstr,outname,/silent   ;; add meta-data on exposure information 
  endelse

  ;; Make sure that it exists, and add to OUTLIST
  outtest = FILE_TEST(outname)
  if outtest eq 1 then outlines = FILE_LINES(outname) else outlines=0
  if (outlines eq 0) then begin
    if outtest eq 0 then print,outname+' NOT FOUND'
    if outtest eq 1 and outlines eq 0 then print,outname+' HAS 0 LINES'
  endif

  if keyword_set(stp) then stop

end
