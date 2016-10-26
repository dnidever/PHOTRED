;+
;
; ADD_COORDS
;
; Add RA/DEC coordinates to a structure using
; the WCS in a given header
;
; INPUTS:
;  str      The IDL structure to add RA/DEC coordinates to.
;             The structure must have X/Y columns
;             (in FITS/IRAF 1-indexed format).
;  headinp  The header string array or the FITS filename
;             
;
; OUTPUTS:
;  RA/DEC columns are added to STR with the ra/dec
;  coordinates using the WCS in the header.
;
; USAGE:
;  IDL>add_coords,str,head
;
; By D. Nidever  Oct 2016
;-

pro add_coords,str,headinp

if n_elements(str) eq 0 or n_elements(headinp) eq 0 then begin
  print,'Syntax - add_coords,str,head/filename'
  return
endif
if size(str,/type) ne 8 then begin
  print,'STR must be a structure'
  return
endif
if tag_exist(str,'X') eq 0 or tag_exist(str,'Y') eq 0 then begin
  print,'STR must have X/Y columns'
  return
endif

; Headinp is likely a filename
if n_elements(headinp) eq 1 then begin
  if file_test(headinp) eq 0 then begin
    print,headinp,' NOT FOUND'
  endif
  head = headfits(headinp)
endif else head = headinp           ; header input

if tag_exist(str,'ra') eq 0 then add_tag,str,'RA',0.0d0,str
if tag_exist(str,'dec') eq 0 then add_tag,str,'DEC',0.0d0,str
HEAD_XYAD,head,str.x-1,str.y-1,ra,dec,/deg
str.ra = ra
str.dec = dec

end
