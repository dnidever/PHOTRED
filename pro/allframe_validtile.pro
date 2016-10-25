;+
;
; ALLFRAME_VALIDTILE
;
; This double-checks if the TILE structure is valid.
;
; INPUTS:
;  tile     The tile structure
;
; OUTPUTS:
;  check    1-if the tile is good and 0-if there is
;             a problem with it or it doesn't exist.
;  =error   The problem if there was one.
;
; USAGE:
;  IDL>check = allframe_validtile(tile)
;
; By D. Nidever  Oct 2016
;-

function allframe_validtile,tile,error=error,silent=silent

undefine,error

; It must exist
if n_elements(tile) eq 0 then begin
  error = 'TILE does not exist'
  return,0
endif

; It must be a structure
if size(tile,/type) ne 8 then begin
  error = 'TILE is not a structure'
  return,0
endif

; Must have TYPE column
if tag_exist(tile,'type') eq 0 then begin
  error = 'TILE does not have TYPE column'
  return,0
endif

; Do we have enough information for each type
CASE strupcase(strtrim(tile.type,2)) of
'ORIG': begin
  ; Don't need anything for ORIG
  return,1
end
'WCS': begin
  ; AST structure exists
  if tag_exist(tile,'AST')

  ; NO AST structure, check other needed values
  endif else begin
    ; Need CRVAL, CRPIX, CTYPE, CDELT
    if tag_exist(tile,'CRVAL') eq 0 then begin
      error = 'TYPE "WCS" requires AST or CRVAL'
       return,0
   endif
    if tag_exist(tile,'CRPIX') eq 0 then begin
      error = 'TYPE "WCS" requires AST or CRPIX'
      return,0
    endif
    if tag_exist(tile,'CTYPE') eq 0 then begin
      error = 'TYPE "WCS" requires AST or CTYPE'
      return,0
    endif
    if tag_exist(tile,'CDELT') eq 0 then begin
      error = 'TYPE "WCS" requires AST or CDELT'
      return,0
    endif
  endelse
end
else: stop,strtrim(tile.type,2)+' not supported yet'
ENDCASE

return,1

end