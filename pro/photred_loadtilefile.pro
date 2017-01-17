;+
;
; PHOTRED_LOADTILEINFO
;
; This programs loads the information for a fields
; tiling scheme.
;
; INPUTS:
;  tilefile  The tiling file, e.g. "F2.tiling"
;  /silent   Don't print anything to the screen.
;  /stp      Stop at the end of the program.
;
; OUTPUTS:
;  tilestr   The structure with information on the tiling
;              scheme/projection and the individual tiles.
;  =error    The error message if one occurred.
;
; USAGE:
;  IDL>photred_loadtilefile,'F2.tiling',tilestr
;
; By D.Nidever  Jan 2017
;-

pro photred_loadtilefile,tilefile,tilestr,error=error,silent=silent,stp=stp

undefine,tilestr
  
; Not enough inputs
if n_elements(tilefile) eq 0 then begin
  print,'Syntax - photred_loadtilefile,tilefile,tilestr,error=error'
  error = 'Not enough inputs'
  return
endif

; File not found
if file_test(tilefile) eq 0 then begin
  error = tilefile+' NOT FOUND'
  if not keyword_set(silent) then print,error
  return
endif

; Load the lines
READLINE,tilefile,tilelines,count=ntilelines

; First part has the information on the tiling scheme/projection
; Second part has information on each tile
undefine,tilestr,tstr
projstr = {cenra:0.0d0,cendec:0.0d0,nx:0L,ny:0L,xstep:0.0d0,ystep:0.0d0,xref:0.0,yref:0.0,ntiles:0L}
tags = tag_names(projstr)
for i=0,ntilelines-1 do begin
   if strtrim(tilelines[i],2) eq '' then goto,BOMB
   
   ;; Tiling scheme/projection information
   if stregex(tilelines[i],'=',/boolean) eq 1 then begin
    dum = strsplit(tilelines[i],'=',/extract)
    name = strtrim(dum[0],2)
    value = strtrim(dum[1],2)
    if valid_num(value,/integer) eq 1 then value=long(value)  ; convert to integer
    if valid_num(value) eq 1 then value=double(value)         ; convert to float
    tind = where(tags eq strupcase(name),ntind)
    ;; Tag found
    if ntind gt 0 then begin
      projstr.(tind[0]) = value

    ;; Tag NOT found, add it
    endif else begin
      projstr = CREATE_STRUCT(projstr,name,value)
    endelse

  ;; Individual Tile information
  endif else begin
    ; ID, NAME, x0, x1, nx, y0, y1, ny, nimages
    dum = strtrim(strsplit(tilelines[i],' ',/extract),2)
    tstr1 = {num:0L,name:'',x0:0L,x1:0L,nx:0L,y0:0L,y1:0L,ny:0L,nimages:0L}
    types = [3,7,3,3,3,3,3,3,3]
    for j=0,n_elements(types)-1 do tstr1.(j)=fix(dum[j],type=types[j])
    PUSH,tstr,tstr1
  endelse
  BOMB:
endfor                          ; loop through tile lines

; XREF/YREF are in IRAF 1-indexed format

; Make header with the WCS
MKHDR,tilehead,fltarr(5,5)
SXADDPAR,tilehead,'NAXIS1',projstr.nx
SXADDPAR,tilehead,'CDELT1',projstr.xstep
SXADDPAR,tilehead,'CRPIX1',projstr.xref
SXADDPAR,tilehead,'CRVAL1',projstr.cenra
SXADDPAR,tilehead,'CTYPE1','RA---TAN'
SXADDPAR,tilehead,'NAXIS2',projstr.ny
SXADDPAR,tilehead,'CDELT2',projstr.ystep
SXADDPAR,tilehead,'CRPIX2',projstr.yref
SXADDPAR,tilehead,'CRVAL2',projstr.cendec
SXADDPAR,tilehead,'CTYPE2','DEC--TAN'

; Combine all of the information
tilestr = CREATE_STRUCT(projstr,'tiles',tstr,'head',tilehead)

if keyword_set(stp) then stp

end
