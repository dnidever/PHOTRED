pro head_adxy,head,a,d,x,y,degree=degree,stp=stp,error=error

;+
;
; HEAD_ADXY
;
; This program uses the WCS in a FITS header to convert
; RA/DEC values (in hours and degrees) to X/Y pixels values.
; Vectors can be input.  Unlike adxy.pro/xyad.pro this also
; works with the MOSAIC TNX WCS system.
;
; INPUTS:
;  head    A FITS header containing a WCS
;  a       RA values in HOURS, or degrees if /degree set.
;  d       DEC values in degrees.
;  /degree RA in degrees, otherwise RA is in HOURS.
;  /stp    Stop at the end of the program.
;
; OUTPUTS:
;  x       X pixel values.  In IDL format, starting at 0.
;  y       Y pixel values.  In IDL format, starting at 0.
;
; USAGE:
;  IDL>head_adxy,head,a,d,x,y,/degree
;
; By D. Nidever   January 2007
;-

nhead = n_elements(head)
na = n_elements(a)
nd = n_elements(d)

r2d = 180.0d0/!dpi

; Not enough inputs
if nhead eq 0 or na eq 0 or nd eq 0 then begin
  print,'Syntax - head_adxy,head,a,d,x,y,degree=degree'
  return
endif

ctype1 = sxpar(head,'CTYPE1',/silent)
; No WCS in header
if strtrim(ctype1,2) eq '0' then begin
  x = -1
  y = -1
  print,'NO WCS IN HEADER'
  error = 'NO WCS IN HEADER'
  return
endif
tnx = stregex(ctype1,'TNX',/boolean,/fold_case)
dum = strsplit(ctype1,'-',/extract)
wcstype = dum[1]

; WCS type
case wcstype of
'TNX': begin
         if keyword_set(degree) then a2=a/r2d else a2=a*15.0d0/r2d
         d2 = d/r2d
         wcs = hdr2wcstnx(head)
         WCSTNX_RD2XY, a2, d2, wcs, x, y    ; Needs ra/dec in RADIANS
       end
'TPV': begin
         if keyword_set(degree) then a2=a/r2d else a2=a*15.0d0/r2d
         d2 = d/r2d
         wcs = hdr2wcstpv(head)
         WCSTPV_RD2XY, a2, d2, wcs, x, y    ; Needs ra/dec in RADIANS
       end
else: ADXY,head,a,d,x,y    ; Needs ra/dec in degrees
endcase

; Return scalars if only one element
if n_elements(x) eq 1 then x=x[0]
if n_elements(y) eq 1 then y=y[0]

if keyword_set(stp) then stop

end
