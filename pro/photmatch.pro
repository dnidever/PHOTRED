;+
;
; This matches up two photometric catalogs using coordinates and
; photometry.
;
; INPUTS:
;  str1        First photometry structure
;  str2        Second photometry structure
;  =dcr        Critical matchup radius.  Stars closer than this will be
;                further evaluated for matching.  dcr=0.5 arcsec by default.
;  =magindarr  The structure field indices for the magnitudes.
;              If this is not input then the structures must be in this format:
;              ID, X, Y, MAG1, ERR1, MAG2, ERR2, ..., CHI, SHARP, other tags
;  /astlinear  Use a linear astrometric fit instead of constant
;  /astquad    Use a quadratic astrometric fit instead of constant
;  /stp        Stop at end of program
;  /silent     Don't print anything
; 
; OUTPUTS:
;  ind1        The matched indices for str1
;  ind2        The matched indices for str2
;  =magoffset  The photometric offsets.  magoff = str1 - str2
;  =magoffsig  The error in the photometric offset.
;  =astoffset  The astrometric constant offsets [RAoff, DECoff] in
;                degrees.  astoff = str1 - str2
;  =count      The number of matches
;  =error      The error message if there was one, else undefined
;
; USAGE:
;  IDL>photmatch,str1,str2,ind1,ind2,dcr=0.5,count=count
;
; By D.Nidever Jan 2007
;-


pro phot_overlap_dummy
FORWARD_FUNCTION trans_coord, trans_coord_dev, stardiff
end

;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function trans_coord,ra,dec,par

new = par[0] + par[1]*ra + par[2]*dec + par[3]*ra*dec + $
      par[4]*ra^2.0 + par[5]*dec^2.0 + par[6]*dec*ra^2.0  + par[7]*ra*dec^2.0
;new = par[0] + par[1]*ra + par[2]*ra^2.0 + par[3]*dec + par[4]*ra*dec + $
;      par[5]*dec*ra^2.0 + par[6]*dec^2.0 + par[7]*ra*dec^2.0
;;new = par[0] + par[1]*ra + par[2]*dec + par[3]*ra*dec

return,new

end

;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function trans_coord_dev,par,ra=ra,dec=dec,y=y

; Transform coordinates and return deviates

newy = trans_coord(ra,dec,par)
return,y-newy

end

;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function stardiff,star1,star2,rms=rms,magindarr1=magindarr1,magindarr2=magindarr2,errindarr1=errindarr1,$
                  errindarr2=errindarr2,magoff=magoff

; Find the difference between stars
; STAR1  need to be a single star
; STAR2  can be an array

; How similar are these stars
; always divide by the "expected" variation in the quantity

; RA/DEC distances
dist = sphdist(star1.ra,star1.dec,star2.ra,star2.dec,/deg)*3600.0
diff = dist/rms

; Compare all of the Photometric BANDS
nmagindarr = n_elements(magindarr)
for k=0,nmagindarr-1 do begin
  magdiff = star1.(magindarr1[k]) - star2.(magindarr1[k]) - magoff[k]
  err = sqrt( star1.(errindarr1[k])^2.0 + star2.(errindarr2[k]) )
  gg = where(abs(magdiff) lt 50. and err lt 8.0,ngg)
  if ngg gt 0 then diff[gg] = sqrt( diff[gg]^2.0 + (magdiff[gg]/err[gg])^2.0 ) 
endfor

; Compare PROB, if it exists
if TAG_EXIST(star1,'PROB') and TAG_EXIST(star2,'PROB') then begin
  probdiff = star1.prob - star2.prob
  diff = sqrt( diff^2.0 + (probdiff/0.2)^2.0 )
endif

return,diff

end

;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro photmatch,str1,str2,ind1,ind2,dcr=dcr,magindarr=magindarr,count=count,$
    stp=stp,silent=silent,error=error,magoffset=magoffset,magoffsig=magoffsig,$
    astlinear=astlinear,astquad=astquad,astoffset=astoffset

undefine,error,magoffset,magoffsig,astoffset
ind1 = -1
ind2 = -1
count = 0

nstr1 = n_elements(str1)
nstr2 = n_elements(str2)

; Not enough inputs
;--------------------
if nstr1 eq 0 or nstr2 eq 0 then begin
  print,'Syntax - photmatch,str1,str2,ind1,ind2,dcr=dcr,magindarr=magindarr,count=count,'
  print,'    stp=stp,silent=silent,error=error,astlinear=astlinear'
  error = 'Not enough inputs'
  return
endif

; Error Handling
;------------------
; Establish error handler. When errors occur, the index of the  
; error is returned in the variable Error_status:  
CATCH, Error_status 

;This statement begins the error handler:  
if (Error_status ne 0) then begin 
   print,'PHOTMATCH ERROR: ', !ERROR_STATE.MSG  
   CATCH, /CANCEL 
   error = !ERROR_STATE.MSG
   return
endif


;--------------------------------
; Checking the input structures
;--------------------------------

; Checking that they are structures
type1 = size(str1,/type)
type2 = size(str2,/type)
if type1 ne 8 then begin
  print,'str1 IS NOT A STRUCTURE'
  error = 'str1 IS NOT A STRUCTURE'
  return
endif
if type2 ne 8 then begin
  print,'str2 IS NOT A STRUCTURE'
  error = 'str2 IS NOT A STRUCTURE'
  return
endif

tags1 = tag_names(str1)
ntags1 = n_elements(tags1)
tags2 = tag_names(str2)
ntags2 = n_elements(tags2)

;; Checking that the structures are the same
;tags1 = tag_names(str1)
;ntags1 = n_elements(tags1)
;tags2 = tag_names(str2)
;ntags2 = n_elements(tags2)
;if ntags1 ne ntags2 then begin
;  print,'DIFFERENT DATA STRUCTURES'
;  error = 'DIFFERENT DATA STRUCTURES'
;  return
;endif
;
;alltags = [tags1,tags2]
;ui = uniq(alltags,sort(alltags))
;nui = n_elements(ui)
;if nui ne ntags1 then begin
;  print,'DIFFERENT DATA STRUCTURES'
;  error = 'DIFFERENT DATA STRUCTURES'
;  return
;endif
tags = tags1

; Make sure there are RA/DEC tags
gdra = where(tags eq 'RA',ngdra)
if ngdra eq 0 then begin
  print,'NO "RA" TAG'
  error = 'NO "RA" TAG'
  return
endif
gddec = where(tags eq 'DEC',ngddec)
if ngddec eq 0 then begin
  print,'NO "DEC" TAG'
  error = 'NO "DEC" TAG'
  return
endif

if not keyword_set(silent) then begin
  print,'-------------------------------------------------------'
  print,'MATCHING STR1 (',strtrim(nstr1,2),' sources) and STR2 (',strtrim(nstr2,2),' sources)'
  print,'-------------------------------------------------------'
endif


;#####################################
;# STARTING THE MATCHING PROCESS

; Find the overlap region
rar1 = minmax(str1.ra)
decr1 = minmax(str1.dec)
rar2 = minmax(str2.ra)
decr2 = minmax(str2.dec)

ind1 = where(str1.ra ge rar2[0] and str1.ra le rar2[1] and $
             str1.dec ge decr2[0] and str1.dec le decr2[1],nind1)
ind2 = where(str2.ra ge rar1[0] and str2.ra le rar1[1] and $
             str2.dec ge decr1[0] and str2.dec le decr1[1],nind2)

; No overlap
if (nind1 eq 0 or nind2 eq 0) then return


; Getting the Magnitude/Error field indices
;------------------------------------------
; NO magindarr input
if n_elements(magindarr) eq 0 then begin

  ; Assuming that the file is in the format:
  ; ID, X, Y, MAG1, ERR1, MAG2, ERR2, ..., CHI, SHARP, other tags
  gdy = where(tags eq 'Y',ngdy)
  if (ngdy eq 0) then begin
    print,'No "Y" TAG FOUND'
    error = 'No "Y" TAG FOUND'
    ind1 = -1
    ind2 = -1
    count = 0
    return
  endif
  gdchi = where(tags eq 'CHI',ngdchi)
  if (ngdchi eq 0) then begin
    print,'NO "CHI" TAG FOUND'
    error = 'NO "CHI" TAG FOUND'
    ind1 = -1
    ind2 = -1
    count = 0
    return
  endif
  lo = gdy[0] + 1         ; first magnitude field is AFTER Y
  hi = gdchi[0] - 1       ; last error field is BEFORE CHI

  ; Need an even number of magnitude/error fields
  if odd(hi-lo+1) eq 1 then begin
    print,'NOT THE RIGHT NUMBER OF MAGNITUDE/ERROR TAGS BETWEEN "Y" AND "CHI"'
    error = 'NOT THE RIGHT NUMBER OF MAGNITUDE/ERROR TAGS BETWEEN "Y" AND "CHI"'
    ind1 = -1
    ind2 = -1
    count = 0
    return
  endif

  nfilters = (hi-lo+1)/2
  magindarr = lindgen(nfilters)*2+lo      ; indices for the magnitudes fields

endif
nmagindarr = n_elements(magindarr)
nfilters = nmagindarr
errindarr = magindarr+1                  ; indices for the error fields
filters = tags[magindarr]

; Make sure the error tag names end in "ERR"
;---------------------------------------------
for i=0,nfilters-1 do begin
  errname = strtrim(tags[errindarr[i]],2)
  len = strlen(errname)
  ending = strmid(errname,len-3,3)
  if strupcase(ending) ne 'ERR' then begin
    print,errname,' DOES NOT END IN "ERR"'
    error = errname+' DOES NOT END IN "ERR"'
    ind1 = -1
    ind2 = -1
    count = 0
    return
  endif
end

; MAKING sure that STR2 has the same MAGNITUDES
; Getting the Magnitude and Error indices for STR2
; They might not be in the same place as in STR1
;-------------------------------------------------
magindarr2 = lonarr(nfilters)-1
errindarr2 = lonarr(nfilters)-1
for i=0,nfilters-1 do begin
  magind = where(tags2 eq tags1[magindarr[i]],nmagind)  
  magindarr2[i] = long(magind[0])
  errind = where(tags2 eq tags1[errindarr[i]],nerrind)  
  errindarr2[i] = long(errind[0])
end

; Some magnitudes missing
bdmag = where(magindarr2 eq -1,nbdmag)
; NO magnitudes in COMMON
if nbdmag gt 0 and nbdmag eq nmagindarr then begin
  print,'NO magnitudes in common.  Using just RA/DEC to match'
  undefine,magindarr,errindarr
  undefine,magindarr2,errindarr2
  nmagindarr = 0
endif
; Some missing, but Some in common
if nbdmag gt 0 then begin
  REMOVE,bdmag,magindarr,errindarr,magindarr2,errindarr2
  nmagindarr = n_elements(magindarr)
  print,'Only ',strtrim(nmagindarr,2),' magnitudes in common'
endif



;-------------------------------
; FIRST MATCHUP WITH COORDINATES
;-------------------------------
if n_elements(dcr) eq 0 then dcr=0.5
SRCMATCH,str1[ind1].ra,str1[ind1].dec,str2[ind2].ra,str2[ind2].dec,dcr,mind1,mind2,count=nmatch,/sph

; No matches
if (nmatch eq 0) then begin
  print,'No matches'
  ind1 = -1
  ind2 = -1
  count = 0
  return
end

; Matched indices
match1 = ind1[mind1]
match2 = ind2[mind2]

if not keyword_set(silent) then begin
  print,strtrim(nmatch,2),' initial matches within ',stringize(dcr,ndec=2),' arcsec'
  print,''
endif


;-------------------------------
; REMOVE ANY ASTROMETRIC OFFSETS
;-------------------------------
str1m = str1[match1]
str2m = str2[match2]

;-- Get only "good" stars --
if (nmagindarr gt 0) then begin

  g1 = where(str1m.(magindarr[0]) lt 50.,ng1)
  if nmagindarr ge 2 then begin
    g1 = where(str1m.(magindarr[0]) lt 50. and str1m.(magindarr[1]) lt 50.,ng1)
  endif
  ; Remove sources with bad CHI values
  if TAG_EXIST(str1m,'CHI') and ng1 gt 0 then begin
    bd = where(str1m[g1].chi gt 1.5,nbd)
    if nbd gt 0 and nbd lt ng1 then REMOVE,bd,g1
    if nbd eq ng1 then undefine,g1
    ng1 = n_elements(g1)
  endif
  ; Remove sources with bad SHARP values
  if TAG_EXIST(str1m,'SHARP') and ng1 gt 0 then begin
    bd = where(abs(str1m[g1].sharp) gt 1.0,nbd)
    if nbd gt 0 and nbd lt ng1 then REMOVE,bd,g1
    if nbd eq ng1 then undefine,g1
    ng1 = n_elements(g1)
 endif

; No magnitudes, use ALL stars
endif else begin
  g1 = findgen(nmatch)
endelse

; Some "good" stars found, REMOVE any astrometric offsets
;--------------------------------------------------------
if (ng1 gt 0) then begin

  ; These ASTROMETRIC fits are ONLY for MATCHING stars
  ; and shouldn't be applied to the NON-MATCHING stars.

  ; Mean offsets BEFORE Astronometric Fit
  if not keyword_set(silent) then $
    print,'--PRE-FIT Astrometric Comparison--'
  raresid1 = (str1m[g1].ra-str2m[g1].ra)*cos(str1m[g1].dec/!radeg)*3600.
  decresid1 = (str1m[g1].dec-str2m[g1].dec)*3600.
  rms1 = sqrt( mean( raresid1^2.0 + decresid1^2.0 ) )
  if ng1 gt 1 then mnradiff1=median(raresid1,/even) else mnradiff1=raresid1[0]
  if ng1 gt 1 then mndecdiff1=median(decresid1,/even) else mndecdiff1=decresid1[0]
  if not keyword_set(silent) then begin
    print,'RMS = ',rms1,'arcsec',format='(A6,F8.5,A7)'
    print,'Mean RA offset:  ',mnradiff1,'arcsec',format='(A17,F8.5,A7)'
    print,'Mean DEC offset: ',mndecdiff1,'arcsec',format='(A17,F8.5,A7)'
  endif

  if not keyword_set(silent) then $
    print,'Fitting astrometric offset'

  ; Constant, linear or quadratic fit?  Constant by default
  ; par = [0, x, y, xy, xx, yy, yxx, xyy]
  parinfo = replicate({fixed:1},8)
  parinfo[0].fixed = 0                                        ; constant, default
  if keyword_set(astlinear) then parinfo[[0,1,2,3]].fixed=0   ; linear
  if keyword_set(astquad) then parinfo.fixed=0                ; quadratic

  ;ra = str1m[g1].ra
  ;dec = str1m[g1].dec
  ra = str2m[g1].ra
  dec = str2m[g1].dec
  radiff = (str1m[g1].ra - str2m[g1].ra)
  decdiff = (str1m[g1].dec - str2m[g1].dec)

  ; Now do a robust fit for the transformation
  func = 'trans_coord_dev'
  fa = {ra:double(ra), dec:double(dec), y:double(radiff)}
  par = dblarr(8)
  fpar_ra = MPFIT(func,par,functargs=fa, perror=perror, niter=iter, status=status1,$
               bestnorm=chisq, dof=dof, autoderivative=1, /quiet, parinfo=parinfo)  ;ftol=1d-10
  fa = {ra:double(ra), dec:double(dec), y:double(decdiff)}
  par = dblarr(8)
  fpar_dec = MPFIT(func,par,functargs=fa, perror=perror, niter=iter, status=status2,$
               bestnorm=chisq, dof=dof, autoderivative=1, /quiet, parinfo=parinfo)  ;ftol=1d-10

  ; MPFIT okay
  if (status1 ge 1 and status2 ge 1) then begin

    ; Print the transformation
    if not keyword_set(silent) then begin
      if not keyword_set(astlinear) and not keyword_set(astquad) then begin
        print,'RA fit:  ',fpar_ra[0],format='(A9,F11.6)'
        print,'DEC fit: ',fpar_dec[0],format='(A9,F11.6)'
      endif
      if keyword_set(astlinear) and not keyword_set(astquad) then begin
        print,'RA fit:  ',fpar_ra[0:3],format='(A9,4F11.6)'
        print,'DEC fit: ',fpar_dec[0:3],format='(A9,4F11.6)'
      endif
      if keyword_set(astquad) then begin
        print,'RA fit:  ',fpar_ra,format='(A9,8F11.6)'
        print,'DEC fit: ',fpar_dec,format='(A9,8F11.6)'
      endif
    endif

    ; Get RMS
    newra = trans_coord(str2m[g1].ra,str2m[g1].dec,fpar_ra) + str2m[g1].ra
    newdec = trans_coord(str2m[g1].ra,str2m[g1].dec,fpar_dec) + str2m[g1].dec
    raresid = (str1m[g1].ra-newra)*cos(str1m[g1].dec/!radeg)*3600.
    decresid = (str1m[g1].dec-newdec)*3600.
    rms = sqrt( mean( raresid^2.0 + decresid^2.0 ) )

    ; Mean offsets AFTER offsetting
    if not keyword_set(silent) then $
      print,'--POST-FIT Astrometric Comparison--'
    if ng1 gt 1 then mnradiff=median(raresid,/even) else mnradiff=raresid[0]
    if ng1 gt 1 then mndecdiff=median(decresid,/even) else mndecdiff=decresid[0]
    if not keyword_set(silent) then begin
      print,'RMS = ',rms,'arcsec',format='(A6,F8.5,A7)'
      print,'Mean RA offset:  ',mnradiff,'arcsec',format='(A17,F8.5,A8)'
      print,'Mean DEC offset: ',mndecdiff,'arcsec',format='(A17,F8.5,A8)'
    endif

    ; Constant Astrometric Offset BEFORE offsetting
    if keyword_set(astlinear) or keyword_set(astquad) then begin
      astoffset = [mnradiff1, mndecdiff1]/3600.0    ; previously derived offset
    endif else begin
      astoffset = [fpar_ra[0], fpar_dec[0]]  ; Use the MPFIT value
    endelse


    ; Get NEW RA/DEC
    ra2 = trans_coord(str2.ra,str2.dec,fpar_ra) + str2.ra
    dec2 = trans_coord(str2.ra,str2.dec,fpar_dec) + str2.dec


  ; MPFIT problem, just find constant offsets
  endif else begin

    ; Constant offsets
    radiff = (str1m[g1].ra - str2m[g1].ra)
    decdiff = (str1m[g1].dec - str2m[g1].dec)
    if ng1 gt 1 then begin
      mnradiff = median(radiff,/even)
      mndecdiff = median(decdiff,/even)
    endif else begin
      mnradiff = radiff
      mndecdiff = decdiff
    endelse

    ; Astrometric constant offset
    astoffset = [mnradiff, mndecdiff]/3600.0

    ; Get RMS
    raresid = radiff*cos(str1m[g1].dec/!radeg)*3600.
    decresid = decdiff*3600. 
    rms = sqrt( mean( raresid^2.0 + decresid^2.0 ) )

    if not keyword_set(silent) then begin
      print,'MPFIT Problem.  Just finding constant offset'
      print,'--POST-FIT Astrometric Comparison--'
      print,'RA constant offset  = ',string(mnradiff,format='(F8.5)')
      print,'DEC constant offset = ',string(mndecdiff,format='(F8.5)')
      print,'RMS = ',rms,'arcsec',format='(A6,F8.5,A7)'
    endif

    ; Get NEW RA/DEC
    ra2 = str2.ra + mnradiff
    dec2 = str2.dec + mndecdiff

  endelse


  ; REDO matching
  ;dcr2 = 5.0*rms > dcr
  dcr2 = dcr
  SRCMATCH,str1[ind1].ra,str1[ind1].dec,ra2[ind2],dec2[ind2],dcr2,mind1,mind2,count=nmatch,/sph

  ; Matched indices
  match1 = ind1[mind1]
  match2 = ind2[mind2]

  if not keyword_set(silent) then $
    print,strtrim(nmatch,2),' matches within ',stringize(dcr2,ndec=2),' arcsec AFTER astrometric offsets'

; NO "good" stars
endif else begin

  if not keyword_set(silent) then print,'No "good" stars.  No astrometric fit or offset.'

  ; Get RMS
  radiff = (str1m.ra - str2m.ra)
  decdiff = (str1m.dec - str2m.dec)
  raresid = radiff*cos(str1m.dec/!radeg)*3600.
  decresid = decdiff*3600.
  rms = sqrt( mean( raresid^2.0 + decresid^2.0 ) )

  ; Mean offsets
  if nmatch gt 1 then mnradiff=median(raresid,/even) else mnradiff=raresid[0]
  if nmatch gt 1 then mndecdiff=median(decresid,/even) else mndecdiff=decresid[0]
  if not keyword_set(silent) then begin
    print,'RMS = ',rms,'arcsec',format='(A6,F8.5,A7)'
    print,'Mean RA offset:  ',mnradiff,'arcsec',format='(A17,F8.5,A8)'
    print,'Mean DEC offset: ',mndecdiff,'arcsec',format='(A17,F8.5,A8)'
  endif

  ; Astrometric offset
  astoffset = [mnradiff, mndecdiff]/3600.0

  ra2 = str2.ra
  dec2 = str2.dec
endelse

;if not keyword_set(silent) then $
;  print,'RMS = ',stringize(rms,ndec=2),' arcsec'


;-------------------------------
; REMOVE ANY PHOTOMETRIC OFFSETS
;-------------------------------
str1m = str1[match1]
str2m = str2[match2]

; Some magnitudes in common
undefine,magoffset
IF (nmagindarr gt 0) then begin

  ;-- Get only "good" stars --
  g1 = where(str1m.(magindarr[0]) lt 50.,ng1)
  if nmagindarr ge 2 then begin
    g1 = where(str1m.(magindarr[0]) lt 50. and str1m.(magindarr[1]) lt 50.,ng1)
  endif
  ; Remove sources with bad CHI values
  if TAG_EXIST(str1m,'CHI') and ng1 gt 0 then begin
    g1_orig = g1
    bd = where(str1m[g1].chi gt 1.5,nbd)
    if nbd gt 0 and nbd lt ng1 then REMOVE,bd,g1
    if nbd eq ng1 then undefine,g1
    ng1 = n_elements(g1)
    if ng1 eq 0 then g1 = g1_orig   ; ALL CHIs bad
    ng1 = n_elements(g1)
  endif
  ; Remove sources with bad SHARP values
  if TAG_EXIST(str1m,'SHARP') and ng1 gt 0 then begin
    bd = where(abs(str1m[g1].sharp) gt 1.0,nbd)
    if nbd gt 0 and nbd lt ng1 then REMOVE,bd,g1
    if nbd eq ng1 then undefine,g1
    ng1 = n_elements(g1)
  endif

  if (ng1 gt 0) then begin

    ; Getting the PHOTOMETRIC OFFSETS
    ;--------------------------------
    if not keyword_set(silent) then begin
      print,''
      print,'Photometric Offsets:'
    endif

    ; Looping through the bands
    ; The offset is STR1 - STR2
    nmagindarr = n_elements(magindarr)
    magoffset = fltarr(nmagindarr)
    magoffsig = fltarr(nmagindarr)
    for i=0,nmagindarr-1 do begin

      magdiff = str1m[g1].(magindarr[i])-str2m[g1].(magindarr2[i])
      magerr = sqrt( str1m[g1].(errindarr[i])^2.0 + str2m[g1].(errindarr2[i])^2.0 )
      gdmag = where(str1m[g1].(magindarr[i]) lt 50. and str2m[g1].(magindarr2[i]) lt 50. and $
                    magerr lt 0.07,ngdmag)
      if ngdmag lt 10. then $   ; not enough points, lower error threshold
        gdmag = where(str1m[g1].(magindarr[i]) lt 50. and str2m[g1].(magindarr2[i]) lt 50. and $
                      magerr lt 0.1,ngdmag)
      if ngdmag lt 10. then $   ; not enough points, lower error threshold
        gdmag = where(str1m[g1].(magindarr[i]) lt 50. and str2m[g1].(magindarr2[i]) lt 50. and $
                      magerr lt 0.2,ngdmag)
      if ngdmag lt 10. then $   ; not enough points, lower error threshold
        gdmag = where(str1m[g1].(magindarr[i]) lt 50. and str2m[g1].(magindarr2[i]) lt 50. and $
                      magerr lt 0.5,ngdmag)

      ; Some sources with decent photometry
      if ngdmag gt 0 then begin

        ;WMEANERR,magdiff[gdmag],magerr[gdmag],magoff1,magoff1err
        ;RESISTANT_MEAN,magdiff[gdmag],2.5,magoff2,magoff2sig,numrej
        ROBUST_MEAN,magdiff[gdmag],magoff1,magoff1err,sig=magerr[gdmag],numrej=numrej

        magoffset[i] = magoff1
        magoffsig[i] = magoff1err/sqrt(ngdmag-numrej)  ; stdev of mean
        com=''
      endif else begin
        ; No good sources, use 0.0
        com=' No good sources' 
      endelse

      ; Print out the offsets
      if not keyword_set(silent) then $
        print,tags1[magindarr[i]],' offset (1-2) = ',strtrim(string(magoffset[i],format='(F10.4)'),2),$
             ' +/- '+strtrim(string(magoffsig[i],format='(F10.4)'),2)+com
    end


    ; CHECK FOR PHOTOMETRIC DIFFERENCES
    ;---------------------------------
    ; Looping through the bands
    bdarr = fltarr(nmagindarr,nmatch) 
    for i=0,nmagindarr-1 do begin

      ; Look for stars with large differences in magnitude
      mag1 = str1m.(magindarr[i])
      mag2 = str2m.(magindarr2[i])
      magdiff = mag1 - mag2 - magoffset[i]
      bd = where(mag1 gt 50. or mag2 gt 50.,nbd)
      if nbd gt 0 then magdiff[bd] = 0.0

      err1 = str1m.(errindarr[i])
      err2 = str2m.(errindarr2[i])
      toterr = sqrt( err1^2.0 + err2^2.0 )

      ; Very different mags but need good photometry for both
      bdmask = long( abs(magdiff) gt 5.0*toterr )
      bdarr[i,*] = bdmask

    end

    ; Check all bands for possible bad matches
    bdall = MAX(bdarr,dim=1)
    bdmatch = where(bdall eq 1,nbdmatch)


    ; REMEMBER THESE COULD ALSO BE VARIABLES!!!

    ; Get "better" matches
    ;---------------------
    if (nbdmatch gt 0) then begin

      ; Loop through the ones with large differences to see if
      ; there is a better match
      for j=0,nbdmatch-1 do begin

        temp = str1m[bdmatch[j]]
        dist = sphdist(temp.ra,temp.dec,ra2,dec2,/deg)*3600.0
        closeind = where(dist lt 5.0*rms,ncloseind)

        ; More than one possibility
        ;--------------------------
        if (ncloseind gt 1) then begin

          ; Compare them
          COMPARE:
          temp2 = str2[closeind]
          temp2.ra = ra2[closeind]
          temp2.dec = dec2[closeind]
          diff = stardiff(temp,temp2,rms=rms,magindarr1=magindarr,errindarr1=errindarr,$
                          magindarr2=magindarr2,errindarr2=errindarr2,magoff=magoffset)

          bestdiff = min(diff)
          best = first_el(minloc(diff))
          bestmatch = closeind[best[0]]

          ; Different match than before
          ;----------------------------
          if (bestmatch ne match2[bdmatch[j]]) then begin
            ;print,'Different match'

            ; Another star from list1 has ALREADY been matched to this star
            ;--------------------------------------------------------------
            otherind = where(match2 eq bestmatch[0],ndoubles)
            otherind = otherind[0]       ; should only be ONE
            if ndoubles gt 0 then begin
              ;print,'DOUBLE MATCHES!!!'

              ; How well does this one compare
              tempother = str1m[otherind]
              tempbest = str2[bestmatch[0]]
              tempbest.ra = ra2[bestmatch[0]]
              tempbest.dec = dec2[bestmatch[0]]
              diff2 = stardiff(tempother,tempbest,rms=rms,magindarr1=magindarr,errindarr1=errindarr,$
                               magindarr2=magindarr2,errindarr2=errindarr2,magoff=magoffset)

              ; OTHER Star matches BETTER
              ;--------------------------
              if (diff2 lt bestdiff) then begin

                ; Remove this star and re-compare
                if n_elements(closeind) gt 1 then begin
                  REMOVE,best,closeind
                  goto,COMPARE
                endif else begin
                  ; Do nothing, leave previous match
                endelse

              ; OTHER Star matches WORSE
              ;-------------------------
              endif else begin

                ; Make this the new match for the star we are currently checking
                match2[bdmatch[j]] = bestmatch[0]

                ; Need to find a new match for the other star
                ;   "tempother", who's match we just "stole".
                ; only use stars that have NOT been matched already
                left = findgen(n_elements(str2))
                gmatch2 = where(match2 ne -1,ngmatch2)
                if ngmatch2 eq n_elements(str2) then goto,BOMB   ; all matched already
                if ngmatch2 gt 0 then $                          ; none matched yet
                  REMOVE,match2[gmatch2],left
                ;dcr2 = 5.0*rms > dcr
                dcr2 = dcr
                SRCMATCH,tempother.ra,tempother.dec,ra2[left],dec2[left],dcr2,indx1,indx2,count=nnewmatch,/sph

                ; New match
                if (nnewmatch gt 0) then begin
                  match2[otherind] = indx2[0]

                ; No good match
                endif else begin
                  match1[otherind] = -1
                  match2[otherind] = -1
                endelse

              endelse ; other star matches WORSE

            ; Not been matched yet
            endif else begin
              match2[bdmatch[j]] = bestmatch[0]    ; NEW match
            endelse

            ;print,'Different match'
            ;plot,str1.ra,str1.dec,ps=1,xr=temp.ra+[-1,1]*0.005,yr=temp.dec+[-1,1]*0.005,xs=1,ys=1
            ;oplot,ra2,dec2,ps=4,co=250
            ;oplot,[temp.ra],[temp.dec],ps=6,sym=5,co=200
            ;stop

          endif  ; different match

          ;print,'More than one possibility'
          ;plot,str1.ra,str1.dec,ps=1,xr=temp.ra+[-1,1]*0.005,yr=temp.dec+[-1,1]*0.005,xs=1,ys=1
          ;oplot,ra2,dec2,ps=4,co=250
          ;oplot,[temp.ra],[temp.dec],ps=6,sym=5,co=200
          ;stop

        ; Only ONE possibility
        endif else begin
          ; Should be the same one we already have matched
        endelse

        BOMB:

      endfor  ; loop

    endif  ; some possible bad matches

  endif  ; some "good" stars

endif  ; some magnitudes in commone


; Removing "bad" matches
bd = where(match1 eq -1,nbd)
if nbd gt 0 then begin
  if nbd eq n_elements(match1) then undefine,match1,match2
  if nbd lt n_elements(match1) then REMOVE,bd,match1,match2
endif


; Final Matches
ind1 = match1
ind2 = match2
count = n_elements(ind1)

if not keyword_set(silent) then $
  print,strtrim(count,2),' final matches'

if keyword_set(stp) then stop

end
