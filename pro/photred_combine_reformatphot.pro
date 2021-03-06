;+
;
; PHOTRED_COMBINE_REFORMATPHOT
;
; Reformat the photometry structure for tile/groups making
; sure that there is a column for each unique exposure.
;
; INTPUTS:
;  phot     The photometry structure output by CALIB.
;  filestr  Structure giving information on all the images
;             included in PHOT.
;  expstr   Structure giving information on all the exposures
;             included in all the tiles/groups for this field.
;
; OUTPUTS:
;  newphot  The new photometr structure with magnitude columns for
;             each uniqe exposure.
;  =error   The error message if one occurred.
;
; USAGE:
;  IDL>photred_combine_reformatphot,phot,filestr,expstr,newphot
;
; By D. Nidever  Jan 2019
;-

pro photred_combine_reformatphot,phot,filestr,expstr,newphot,error=error

;; Create a structure that has a separate column for each EXPOSURE

nphot = n_elements(phot)
nfiles = n_elements(filestr)
nexp = n_elements(expstr)
; Not enough inputs
if nphot eq 0 or nfiles eq 0 or nexp eq 0 then begin
  error = 'Not Enough Inputs'
  print,'Syntax - photred_combine_reformatphot,phot,filestr,expstr,newphot,error=error'
  return
endif

;; Deal with
;;  1) instrumental individual, e.g. I_R1, I_R1ERR, I_I2, I_I2ERR
;;  2) calibrated individual magnitudes, e.g. RMAG1, R1ERR, GMAG2,
;;       G2ERR.  We always have these.
;; 3) calibrated average magnitudes, e.g. ZMAG, ZERR, UMAG, UMAG, UERR

phtags = tag_names(phot)
;; Do we have instrumental magnitudes?
instmag = (total(stregex(phtags,'^I_',/boolean) eq 1) gt 0)
;; Do we have calibrated AVERAGE magnitudes?
avgmag = (total(stregex(phtags,'MAG$',/boolean) eq 1) gt 0)

;; Get the unique EXPOSURE filters
expfiltindex = create_index(expstr.filter)
nexpfilt = n_elements(expfiltindex.value)
;; Unique FILE filters
filefiltindex = create_index(filestr.filter)
nfilefilt = n_elements(filefiltindex.value)

;; Instrumental magnitudes
;;------------------------
;;  e.g. I_R1, I_R1ERR, I_I2, I_I2ERR 
if instmag eq 1 then begin
  ;; Create the schema, ALL exposures in chronological order
  ;;   even if they are not represented in this tile
  ;; EXPSTR is already in chronological order
  ;; Loop over filters
  magnames = strarr(nexp)
  errnames = strarr(nexp)
  for i=0,nexpfilt-1 do begin
    filt = strupcase(expfiltindex.value[i])
    ind = expfiltindex.index[expfiltindex.lo[i]:expfiltindex.hi[i]]
    nind = expfiltindex.num[i]
    magnames[ind] = 'I_'+filt+strtrim(lindgen(nind)+1,2)
    errnames[ind] = 'I_'+filt+strtrim(lindgen(nind)+1,2)+'ERR'
  endfor
  for i=0,nexp-1 do begin
    if i eq 0 then inst_schema = create_struct(magnames[0],0.0) else $
      inst_schema = create_struct(inst_schema,magnames[i],0.0)
    inst_schema = create_struct(inst_schema,errnames[i],0.0)    
  endfor
  inst_phot = replicate(inst_schema,nphot)
  ;; All PHOT instrumental MAG/ERR column indices in order
  instphmagind = where(stregex(phtags,'^I_',/boolean) eq 1 and stregex(phtags,'ERR',/boolean) eq 0,ninstphmagind)
  instpherrind = where(stregex(phtags,'^I_',/boolean) eq 1 and stregex(phtags,'ERR',/boolean) eq 1,ninstpherrind)
  ;; Stuff in the instrumental photometry
  totalfluxwt = dblarr(nphot,nexp)  ; temporary array for calculations
  totalwt = dblarr(nphot,nexp)
  for i=0,nfiles-1 do begin
    ; Which column are we using
    expind = where(expstr.expnum eq filestr[i].expnum,nexpind)
    magind = instphmagind[i]  ; input indices
    errind = instpherrind[i]
    omagind = expind*2        ; output indices
    oerrind = magind+1
    ;; Take weighed means in case there is overlap of some of the chip images
    gd = where(phot.(magind) lt 50,ngd)
    totalfluxwt[gd,expind[0]] += 2.5118864d^phot[gd].(magind) * (1.0d0/phot[gd].(errind)^2)
    totalwt[gd,expind[0]] += 1.0d0/phot[gd].(errind)^2
  endfor
  newflux = totalfluxwt/totalwt
  newmag = 2.50*alog10(newflux)
  newerr = sqrt(1.0/totalwt)
  bdmag = where(finite(newmag) eq 0,nbdmag)
  if nbdmag gt 0 then begin
    newmag[bdmag] = 99.99
    newerr[bdmag] = 9.99
  endif
  ;; Put in the structure
  for i=0,nexp-1 do begin
    inst_phot.(i*2) = newmag[*,i]   ; mags
    inst_phot.(i*2+1) = newerr[*,i]   ; errs
  endfor
endif


;; Calibrated magnitudes
;;----------------------
;; e.g. RMAG1, R1ERR, GMAG2, G2ERR
;; NOTE, if there's only ONE image/exposure for a band
;;  the it has the XMAG and XERR name WITHOUT a 1 at the end!!!
;; Create the schema, ALL exposures in chronological order
;;   even if they are not represented in this tile
;; EXPSTR is already in chronological order
;; Loop over filters
magnames = strarr(nexp)
errnames = strarr(nexp)
for i=0,nexpfilt-1 do begin
  filt = strupcase(expfiltindex.value[i])
  ind = expfiltindex.index[expfiltindex.lo[i]:expfiltindex.hi[i]]
  nind = expfiltindex.num[i]
  magnames[ind] = filt+'MAG'+strtrim(lindgen(nind)+1,2)
  errnames[ind] = filt+strtrim(lindgen(nind)+1,2)+'ERR'
endfor
for i=0,nexp-1 do begin
  if i eq 0 then calib_schema = create_struct(magnames[0],0.0) else $
    calib_schema = create_struct(calib_schema,magnames[i],0.0)
  calib_schema = create_struct(calib_schema,errnames[i],0.0)    
endfor
calib_phot = replicate(calib_schema,nphot)

;; All PHOT calibrated MAG/ERR column indices in order
inpmagnames = strarr(nfiles)
inperrnames = strarr(nfiles)
for i=0,nfilefilt-1 do begin
  filt = strupcase(filefiltindex.value[i])
  ind = filefiltindex.index[filefiltindex.lo[i]:filefiltindex.hi[i]]
  nind = filefiltindex.num[i]
  ;; If there is only ONE observation in a band then the name
  ;; is XMAG and XERR instead of XMAG1 and X1ERR, but it is in
  ;; the right order with the rest of the calibrated magnitude columns
  ;; In the output photometry structure, EVERY exposure gets a
  ;; calibrated magnitude column with a number.
  if nind eq 1 then begin
    inpmagnames[ind] = filt+'MAG'
    inperrnames[ind] = filt+'ERR'
  endif else begin
    inpmagnames[ind] = filt+'MAG'+strtrim(lindgen(nind)+1,2)
    inperrnames[ind] = filt+strtrim(lindgen(nind)+1,2)+'ERR'
  endelse
endfor
calibphmagind = lonarr(nfiles)
calibpherrind = lonarr(nfiles)
for i=0,nfiles-1 do begin
  magind = where(phtags eq inpmagnames[i],nmagind)
  calibphmagind[i] = magind[0]
  errind = where(phtags eq inperrnames[i],nerrind)
  calibpherrind[i] = errind[0]
endfor

;; Stuff in the calibrated photometry
totalfluxwt = dblarr(nphot,nexp)  ; temporary array for calculations
totalwt = dblarr(nphot,nexp)
for i=0,nfiles-1 do begin
  ;; Which column are we using
  expind = where(expstr.expnum eq filestr[i].expnum,nexpind)
  magind = calibphmagind[i]  ; input indices
  errind = calibpherrind[i]
  omagind = expind*2        ; output indices
  oerrind = magind+1
  ;; Take weighed means in case there is overlap of some of the chip images
  gd = where(phot.(magind) lt 50,ngd)
  totalfluxwt[gd,expind[0]] += 2.5118864d^phot[gd].(magind) * (1.0d0/phot[gd].(errind)^2)
  totalwt[gd,expind[0]] += 1.0d0/phot[gd].(errind)^2
endfor
newflux = totalfluxwt/totalwt
newmag = 2.50*alog10(newflux)
newerr = sqrt(1.0/totalwt)
bdmag = where(finite(newmag) eq 0,nbdmag)
if nbdmag gt 0 then begin
  newmag[bdmag] = 99.99
  newerr[bdmag] = 9.99
endif
;; Put in the structure
for i=0,nexp-1 do begin
  calib_phot.(i*2) = newmag[*,i]   ; mags
  calib_phot.(i*2+1) = newerr[*,i]   ; errs
endfor

;; Average calibrated magnitudes
;;------------------------------
;; e.g. ZMAG, ZERR, UMAG, UMAG, UERR
if avgmag eq 1 then begin
  ;; Don't need to do any thing, just use the average magnitudes that
  ;; are already there

  ;; All PHOT average calibrated MAG/ERR column indices in order
  avgphmagind = where(stregex(phtags,'MAG$',/boolean) eq 1,ncalibphmagind)
  avgpherrind = avgphmagind+1   ;; errors should be right after mags
  ;; Make the structure
  magnames = strupcase(expfiltindex.value)+'MAG'
  errnames = strupcase(expfiltindex.value)+'ERR'
  for i=0,nexpfilt-1 do begin
    if i eq 0 then avg_schema = create_struct(magnames[0],99.99) else $
      avg_schema = create_struct(avg_schema,magnames[i],99.99)
    avg_schema = create_struct(avg_schema,errnames[i],9.99)    
  endfor
  avg_phot = replicate(avg_schema,nphot)
  ;; Copy over the values
  STRUCT_ASSIGN,phot,avg_phot,/nozero
endif

;; Now combine everything in the final structure
;;----------------------------------------------

;; Make the schema. The order:
;; - pre-photometry columns, e.g. ID, X, Y
;; - instrumental photometry columns, e.g. I_R1, I_R1ERR
;; - calibrated photometry columns, e.g. RMAG1, RMAG1
;; - average calibrated photometry columns, e.g. UMAG, UERR, IMAG, IERR
;; - post-photometry collumns, e.g. chi, sharp, ra, dec
;;
;; Pre-photometry columns
if instmag eq 1 then hi=instphmagind[0]-1 else hi=calibphmagind[0]-1
fieldnames = phtags[0:hi]
fieldtypes = strarr(hi+1)
for i=0,hi do fieldtypes[i]=size(phot[0].(i),/type)
fschema = create_struct(fieldnames[0],fix(phot[0].(0),type=fieldtypes[0]))
for i=1,hi do fschema = create_struct(fschema,fieldnames[i],fix(phot[0].(i),type=fieldtypes[i]))
struct_assign,{dum:''},fschema  ; zero it out
;; Photometry columns
if instmag eq 1 then fschema = create_struct(fschema,inst_schema)  ;; Instrumental photometry columns
fschema = create_struct(fschema,calib_schema)                      ;; Calibrated photometry columns
fschema = create_struct(fschema,avg_schema)                        ;; Average calibrated photometry columns
;; Post-photometry columns
if avgmag eq 1 then lo=max(avgpherrind)+1 else lo=max(calibphmagind)+1
hi = n_elements(phtags)-1
fieldnames = phtags[lo:hi]
fieldtypes = strarr(hi-lo+1)
for i=lo,hi do fieldtypes[i-lo]=size(phot[0].(i),/type)
for i=lo,hi do fschema = create_struct(fschema,fieldnames[i-lo],fix(phot[0].(i),type=fieldtypes[i-lo]))
struct_assign,{dum:''},fschema  ; zero it out
ftags = tag_names(fschema)

;; Create the final structure
newphot = replicate(fschema,nphot)
;; Stuff in the information
if instmag eq 1 then hi=instphmagind[0]-1 else hi=calibphmagind[0]-1
for i=0,hi do newphot.(i)=phot.(i)    ;; pre-photometry columns
if instmag eq 1 then struct_assign,inst_phot,newphot,/nozero
struct_assign,calib_phot,newphot,/nozero
struct_assign,avg_phot,newphot,/nozero
if avgmag eq 1 then lo=max(avgpherrind)+1 else lo=max(calibphmagind)+1
hi = n_elements(phtags)-1
for i=lo,hi do begin
  ind = where(ftags eq phtags[i],nind)
  newphot.(ind) = phot.(i)
endfor

end
