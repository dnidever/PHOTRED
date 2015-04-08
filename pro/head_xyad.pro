pro head_xyad,head,x,y,a,d,degree=degree,stp=stp,error=error

;+
;
; HEAD_XYAD
;
; This program uses the WCS in a FITS header to convert
; X/Y pixel values to RA/DEC values (in hours and degrees).
; Vectors can be input.  Unlike adxy.pro/xyad.pro this also
; works with the MOSAIC TNX WCS system.
;
; INPUTS:
;  head    A FITS header containing a WCS
;  x       X pixel values.  In IDL convention, starting at 0.
;  y       Y pixel values.  In IDL convention, starting at 0.
;  /degree RA will be in degrees, otherwise RA is in HOURS.
;  /stp    Stop at the end of the program.
;
; OUTPUTS:
;  a       RA values in HOURS unless /degree set.
;  d       DEC values in degrees.
;
; USAGE:
;  IDL>head_xyad,head,x,y,a,d,/degree
;
; By D. Nidever   January 2007
;-


nhead = n_elements(head)
nx = n_elements(x)
ny = n_elements(y)

r2d = 180.0d0/!dpi

; Not enough inputs
if nhead eq 0 or nx eq 0 or ny eq 0 then begin
  print,'Syntax - head_xyad,head,x,y,a,d,degree=degree'
  return
endif

ctype1 = sxpar(head,'CTYPE1',/silent)
; No WCS in header
if strtrim(ctype1,2) eq '0' then begin
  a = -1
  d = -1
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
         wcs = hdr2wcstnx(head)
         WCSTNX_XY2RD, x, y, wcs, a2, d2     ; returns RA/DEC in radians!
         a = a2*r2d   ; convert to degrees
         d = d2*r2d   ; convert to degrees
       end
'TPV': begin
         wcs = hdr2wcstpv(head)
         WCSTPV_XY2RD, x, y, wcs, a2, d2     ; returns RA/DEC in radians!
         a = a2*r2d   ; convert to degrees
         d = d2*r2d   ; convert to degrees
       end
else: XYAD,head,x,y,a,d
endcase

;; WCS is not TNX
;if tnx eq 0 then begin
;  XYAD,head,x,y,a,d
;
;; WCS IS TNX
;endif else begin
;  wcs = hdr2wcstnx(head)
;  WCSTNX_XY2RD, x, y, wcs, a2, d2     ; returns RA/DEC in radians!
;  a = a2*r2d   ; convert to degrees
;  d = d2*r2d   ; convert to degrees
;endelse

if not keyword_set(degree) then a=a/15.0

; Return scalars if only one element
if n_elements(a) eq 1 then a=a[0]
if n_elements(d) eq 1 then d=d[0]

if keyword_set(stp) then stop

end
