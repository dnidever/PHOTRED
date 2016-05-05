pro phot_overlap,inpstr1,inpstr2,outstr,dcr=dcr,stp=stp,silent=silent,$
    s1=s1,s2=s2,magindarr=magindarr,error=error,posonly=posonly

;+
;
; This combines overlapping photometry
;
; INPUTS:
;  str1        First photometry structure
;  str2        Second photometry structure
;  =dcr        Critical matchup radius.  Stars closer than are combined,
;                dcr=0.5 by default.
;  =magindarr  The structure field indices for the magnitudes.
;              If this is not input then the structures must be in this format:
;              ID, X, Y, MAG1, ERR1, MAG2, ERR2, ..., CHI, SHARP, other tags
;  /posonly    Only use astrometric positions for the matching.
;  /stp        Stop at end of program
;  /silent     Don't print anything
; 
; OUTPUTS:
;  outstr      The final combined photometry structure
;  =error      The error message if there was one, else undefined
;
; USAGE:
;  IDL>phot_overlap,str1,str2,outstr
;
; By D.Nidever Jan 2007
;-

undefine,outstr,error

nstr1 = n_elements(inpstr1)
nstr2 = n_elements(inpstr2)

; Not enough inputs
;--------------------
if nstr1 eq 0 or nstr2 eq 0 then begin
  print,'Syntax - phot_overlap,str1,str2,outstr,stp=stp'
  return
endif

; Error Handling
;------------------
; Establish error handler. When errors occur, the index of the  
; error is returned in the variable Error_status:  
CATCH, Error_status 

;This statement begins the error handler:  
if (Error_status ne 0) then begin 
   print,'PHOT_OVERLAP ERROR: ', !ERROR_STATE.MSG  
   undefine,outstr
   CATCH, /CANCEL 
   error = !ERROR_STATE.MSG
   return
endif


; Copy to temporary structures
str1 = inpstr1
str2 = inpstr2


;--------------------------------
; Checking the input structures
;--------------------------------

; Checking that they are structures
type1 = size(str1,/type)
type2 = size(str2,/type)
if type1 ne 8 then begin
  print,'str1 IS NOT A STRUCTURE'
  return
endif
if type2 ne 8 then begin
  print,'str2 IS NOT A STRUCTURE'
  return
endif

; Checking that the structures are the same
tags1 = tag_names(str1)
ntags1 = n_elements(tags1)
tags2 = tag_names(str2)
ntags2 = n_elements(tags2)
if ntags1 ne ntags2 then begin
  print,'DIFFERENT DATA STRUCTURES'
  return
endif

alltags = [tags1,tags2]
ui = uniq(alltags,sort(alltags))
nui = n_elements(ui)
if nui ne ntags1 then begin
  print,'DIFFERENT DATA STRUCTURES'
  return
endif
tags = tags1

; Check the tag types
type1 = lonarr(ntags1)
for i=0,ntags1-1 do type1[i]=size(str1[0].(i),/type)
type2 = lonarr(ntags2)
for i=0,ntags2-1 do type2[i]=size(str2[0].(i),/type)
; Use same types for both
if total(abs(type1-type2)) gt 0 then begin
  print,'TYPES NOT THE SAME. MAKING THEM THE SAME'

  for i=0,ntags1-1 do begin
    itype = max([type1[i],type2[i]])
    if itype eq 7 then zero='' else zero=fix(0,type=itype)
    if i eq 0 then dum=create_struct(tags1[i],zero) else $
      dum=create_struct(dum,tags1[i],zero)
  endfor
  new1 = REPLICATE(dum,nstr1)
  STRUCT_ASSIGN,str1,new1        ; this copies everything properly
  str1 = new1

  for i=0,ntags2-1 do begin
    itype = max([type1[i],type2[i]])
    if itype eq 7 then zero='' else zero=fix(0,type=itype)
    if i eq 0 then dum=create_struct(tags2[i],zero) else $
      dum=create_struct(dum,tags2[i],zero)
  endfor
  new2 = REPLICATE(dum,nstr2)
  STRUCT_ASSIGN,str2,new2        ; this copies everything properly
  str2 = new2
endif


; If the tags are not in the same order
; copy STR2 to STR1's structure type
tagsmatch = where(tags1 eq tags2,ntagsmatch)
if ntagsmatch ne ntags1 then begin
  print,'TAGS ARE NOT IN THE SAME ORDER.  COPYING STR2 to STR1 structure type'

  ; Sometimes the filters are in different orders
  ; which can mess things up when concatenating
  ; e.g. M gets copied into D.
  ; For some reason IDL does NOT produce an error.

  new = REPLICATE(str1[0],nstr2)
  STRUCT_ASSIGN,str2,new        ; this copies everything properly
  str2 = new

endif

; Make sure there are RA/DEC tags
gdra = where(tags eq 'RA',ngdra)
if ngdra eq 0 then begin
  print,'NO "RA" TAG'
  return
endif
gddec = where(tags eq 'DEC',ngddec)
if ngddec eq 0 then begin
  print,'NO "DEC" TAG'
  return
endif


; Getting the Magnitude/Error field indices
;------------------------------------------
; NO magindarr input
if n_elements(magindarr) eq 0 then begin

  ; Assuming that the file is in the format:
  ; ID, X, Y, MAG1, ERR1, MAG2, ERR2, ..., CHI, SHARP, other tags
  gdy = where(tags eq 'Y',ngdy)
  if (ngdy eq 0) then begin
    print,'No "Y" TAG FOUND'
    return
  endif
  gdchi = where(tags eq 'CHI',ngdchi)
  if (ngdchi eq 0) then begin
    print,'NO "CHI" TAG FOUND'
    return
  endif
  lo = gdy[0] + 1         ; first magnitude field is AFTER Y
  hi = gdchi[0] - 1       ; last error field is BEFORE CHI

  ; Need an even number of magnitude/error fields
  if odd(hi-lo+1) eq 1 then begin
    print,'NOT THE RIGHT NUMBER OF MAGNITUDE/ERROR TAGS BETWEEN "Y" AND "CHI"'
    return
  endif

  nfilters = (hi-lo+1)/2
  magindarr = lindgen(nfilters)*2+lo      ; indices for the magnitudes fields

endif else begin
  nfilters = n_elements(magindarr)
endelse
errindarr = magindarr+1                  ; indices for the error fields
filters = tags[magindarr]

; Make sure the error tag names end in "ERR"
;---------------------------------------------
for i=0,nfilters-1 do begin
  errname = strtrim(tags[errindarr[i]],2)
  len = strlen(errname)
  ending = strmid(errname,len-3,3)
  if strupcase(ending) ne 'ERR' then begin
    print,errname,' DOES NOT ENDING IN "ERR"'
    return
  endif
end


;#####################################
;# STARTING THE COMBINATION PROCESS
if n_elements(dcr) eq 0 then dcr=0.5
PHOTMATCH,str1,str2,match1,match2,dcr=dcr,count=nmatch,silent=silent,magindarr=magindarr,posonly=posonly


; No overlap, just concatenate
if (nmatch eq 0) then begin
  outstr = [str1,str2]
  noutstr = n_elements(outstr)
  if not keyword_set(silent) then print,'NO OVERLAP'
  if not keyword_set(silent) then $
    print,strtrim(noutstr,2),' stars'
  return
endif



;-------------------
; Combine the data
;-------------------
if not keyword_set(silent) then $
  print,'Combining ',strtrim(nmatch,2),' Stars'

; This is how the stars will be ORDERED
; [ str1 (plus matched stars), str2 (stars that did NOT match) ]

; Left-over indices for STR2
;----------------------------
leftind2 = lindgen(nstr2)
if (nmatch lt nstr2) then begin
  REMOVE,match2,leftind2
  leftstr2 = str2[leftind2]
endif else begin
  undefine,leftstr2    ; no stars left
endelse
nleft2 = n_elements(leftstr2)


; Structures of the two MATCHED observations
;---------------------------------------------
s1 = str1[match1]
s2 = str2[match2]

; Starting the output structure
;---------------------------------
; Use STR1 to start with
; The output will inherit from STR1: ID, CHI, SHARP, and other tags
out = s1

; Average the RA/DEC coordinates
;--------------------------------
out.ra = (s1.ra+s2.ra)*0.5
out.dec = (s1.dec+s2.dec)*0.5

; Average the X/Y coordinates
;-----------------------------
; Make sure there are X/Y tags
gdx = where(tags eq 'X',ngdx)
gdy = where(tags eq 'Y',ngdy)
if (ngdx gt 0 and ngdy gt 0) then begin
  out.x = (s1.x+s2.x)*0.5
  out.y = (s1.y+s2.y)*0.5
endif


;-------------------------------
; Looping through the filters
;-------------------------------
for j=0,nfilters-1 do begin

  ; Filter index
  magind = magindarr[j]
  errind = magind+1    ; assume the error is the next column

  ; Seeing which ones are "good" and which ones are "bad (i.e. 99.99)
  bothbad = where(s1.(magind) gt 50. and s2.(magind) gt 50.,nbothbad)
  onegood1 = where(s1.(magind) lt 50. and s2.(magind) gt 50.,nonegood1)
  onegood2 = where(s1.(magind) gt 50. and s2.(magind) lt 50.,nonegood2)
  bothgood = where(s1.(magind) lt 50. and s2.(magind) lt 50.,nbothgood)


  ; Both bad
  if (nbothbad gt 0) then begin
    out[bothbad].(magind) = 99.9999
    out[bothbad].(errind) = 9.9999
  endif

  ; Only one good (first one)
  if (nonegood1 gt 0) then begin
    out[onegood1].(magind) = s1[onegood1].(magind)
    out[onegood1].(errind) = s1[onegood1].(errind)
  endif

  ; Only one good (second one)
  if (nonegood2 gt 0) then begin
    out[onegood2].(magind) = s2[onegood2].(magind)
    out[onegood2].(errind) = s2[onegood2].(errind)
  endif

  ; Both good
  if (nbothgood gt 0) then begin
    flux1 = 2.511864^s1[bothgood].(magind)
    flux2 = 2.511864^s2[bothgood].(magind)
    err1 = s1[bothgood].(errind)
    err2 = s2[bothgood].(errind)
    wt1 = 1.0/(err1^2.0)
    wt2 = 1.0/(err2^2.0)
    totalwt = wt1 + wt2
    totalflux = (flux1*wt1 + flux2*wt2)
    totalerr = (err1^2.0)*wt1 + (err2^2.0)*wt2
    newflux = totalflux/totalwt
    newmag = 2.50*alog10(newflux)
    newerr = sqrt(1.0/totalwt)

    out[bothgood].(magind) = newmag
    out[bothgood].(errind) = newerr
  endif

  ; Should we combine the CHI/SHARP at all??

end  ; filter loop


; Making the final output structure
;-----------------------------------
outstr = str1
outstr[match1] = out             ; put matched stars in STR1 positions
if (nleft2 gt 0) then $
  outstr = [outstr,leftstr2]     ; add-on non-matched STR2 stars
noutstr = n_elements(outstr)

if not keyword_set(silent) then $
  print,strtrim(noutstr,2),' stars'

if keyword_set(stp) then stop

end
