;+
;
; PHOTRED_PICKREFERENCEFRAME
;
; Pick the reference frame for a set of images.
;
; INPUTS:
;  base       Set of FITS file base names.  This should be a unique
;               set of exposure names, i.e. if there are multiple
;               chips per exposure then only give the file names
;               for one chip (e.g., *_01).
;  filtref    The reference filter name.  Can be an array of filter
;               names in decreasing order of preference.
;  thisimager The structure with information on this imager.
;  =logfile   The log file name.
;  /fake      This is running in FAKERED, artificial star tests.
;
; OUTPUTS:
;  refstr     Structure with the reference base name, filter and exptime.
;  =error     The error message if one occurred.
;
; USAGE:
;  IDL> photred_pickreferenceframe,base,refimbase
;
; By D.Nidever  Jan 2017
;-

pro photred_pickreferenceframe,base0,filtref,thisimager,refstr,logfile=logfile,fake=fake,error=error

undefine,refimbase,error
  
; Not enough inputs
if n_elements(base0) eq 0 or n_elements(filtref) eq 0 then begin
  print,'Syntax - photred_pickreferenceframe,base,filtref,thisimager,refstr,logfile=logfile,fake=fake'
  error = 'Not enough inputs'
  return
endif

if n_elements(logfile) eq 0 then logfile=-1
  
; Getting the REFERENCE Image
;----------------------------

;; Using chip 1 for multi-chip imagers
nbase0 = n_elements(base0)
if thisimager.namps gt 1 then begin
  amp = strarr(nbase0)
  for k=0,nbase0-1 do begin
    dum = strsplit(base0[k],thisimager.separator,/extract)
    amp[k] = first_el(dum,/last)
  endfor
  uiamp = uniq(amp,sort(amp))
  amps = amp[uiamp]
  MATCH,amp,amps[0],chip1ind,ind2,/sort
  base = base0[chip1ind]
endif else base=base0

nbase = n_elements(base)
nfiltref = n_elements(filtref)
if not keyword_set(fake) then begin
  ; Get filters for first amp
  filters = PHOTRED_GETFILTER(base+'.fits')
  exptime = PHOTRED_GETEXPTIME(base+'.fits')
  rexptime = round(exptime*10)/10.  ; rounded to nearest 0.1s
  utdate = PHOTRED_GETDATE(base+'.fits')
  uttime = PHOTRED_GETUTTIME(base+'.fits')
  dateobs = utdate+'T'+uttime
  jd = dblarr(nbase)
  for l=0,nbase-1 do jd[l]=DATE2JD(dateobs[l])
      
  ; Find matches to the reference filter in priority order
  ngdref=0 & refind=-1
  repeat begin
    refind++
    gdref = where(filters eq filtref[refind],ngdref)
  endrep until (ngdref gt 0) or (refind eq nfiltref-1)
  if ngdref gt 0 then usefiltref=filtref[refind]
  ;gdref = where(filters eq filtref,ngdref)
  ; No reference filters
  if ngdref eq 0 then begin
    printlog,logfile,'NO IMAGES IN REFERENCE FILTER - '+filtref
    printlog,logfile,'MODIFY >>photred.setup<< file parameter FILTREF'
    printlog,logfile,'FILTERS AVAILABLE: '+filters[uniq(filters,sort(filters))]
    error = 'No images in reference filter = '+filtref
    ;printlog,logfile,'Failing field '+thisfield+' and going to the next'
    ;PUSH,failurelist,dirs[i]+'/'+base+'.als'
    ;goto,BOMB
    return
  endif

  ; More than one exposure in reference filter
  ; Use image with LONGEST exptime
  if ngdref gt 1 then begin
    ; Getting image with longest exptime
    refbase = base[gdref]
    maxind = maxloc(rexptime[gdref])
    if n_elements(maxind) gt 1 then begin  ; pick chronological first image
      si = sort(jd[gdref[maxind]])
      maxind = maxind[si[0]]  
    endif 
    ;exptime2 = exptime[gdref]
    ;maxind = first_el(maxloc(exptime2))
    refimbase = refbase[maxind[0]]
    refexptime = exptime[gdref[maxind[0]]]

    printlog,logfile,'Two images in reference filter.'
    printlog,logfile,'Picking the image with the longest exposure time'

  ; Single frame
  endif else begin
    refimbase = base[gdref[0]]
    refexptime = exptime[gdref[0]]
  endelse

  ; Getting just the base, without the extension, e.g. "_1"
  arr = strsplit(refimbase,thisimager.separator,/extract)
  amp = first_el(arr,/last)
  len = strlen(refimbase)
  lenend = strlen(thisimager.separator+amp)
  refimbase = strmid(refimbase,0,len-lenend)
  refstr = {base:refimbase,filter:usefiltref,exptime:float(refexptime)}
  
; FAKE, pick reference image of existing MCH file
;  This ensures that we use exactly the same reference frame.
endif else begin
  filters = PHOTRED_GETFILTER(base+'.fits')
  exptime = PHOTRED_GETEXPTIME(base+'.fits')
  gdref = where(file_test(base+'.mch') eq 1,ngdref)
  if ngdref eq 0 then begin
    error = '/FAKE, no existing MCH file.'
    printlog,logfile,error
    return
  endif
  if ngdref gt 1 then begin
     error = '/FAKE, '+strtrim(ngdref,2)+' MCH files. Too many!'
     printlog,logfile,error
    return
  endif
  refimbase = base[gdref[0]]
  usefiltref = filters[gdref[0]]
  refexptime = exptime[gdref[0]]
  refstr = {base:refimbase,filter:usefiltref,exptime:float(refexptime)}
endelse
        
; Reference image information
printlog,logfile,'REFERENCE IMAGE = '+refstr.base+' Filter='+refstr.filter+' Exptime='+strtrim(refstr.exptime,2)

;stop

end
