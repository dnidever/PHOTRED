;+
; NAME:
;	interp
; PURPOSE: (one line)
;	One dimensional interpolation onto a new x grid (both irregular)
; DESCRIPTION:
; CATEGORY:
;	Mathematical
; CALLING SEQUENCE:
;	interp,x1,y1,x2,y2,e2,E1=e1
; INPUTS:
;	x1 - Input X array.
;	y1 - Matching input Y array.
;	x2 - New X array to interpolate Y to.
; OPTIONAL INPUT PARAMETERS:
; KEYWORD PARAMETERS:
;	E1 - Array of uncertainties for Y.  Default=equal uncertainties.
; OUTPUTS:
;	y2 - Interpolated values of y1 that match the x2 array.
;	e2 - Array of uncertainties for interpolated values.  For each
;          interpolation point, uncertainties are calculated using the
;          left and the right uncertainties.  The minimum of these
;          becomes the final uncertainty.  Returned only if E1 is specified.
; COMMON BLOCKS:
; SIDE EFFECTS:
; RESTRICTIONS:
; PROCEDURE:
; MODIFICATION HISTORY:
;	Written by Marc W. Buie, Lowell Observatory, 1992 Sep 28.
;    5/19/93, DWL, Modified the uncertainty computation.
;   97/11/13, MWB, Moved e1 to a keyword, uncertainties are now optional.
;   98/4/2, MWB, Fixed bug in returned e2 value.
;   2000/01/31, MWB, Fixed minor bug in loop index, now converted to long.
;-
pro interp,x1,y1,x2,y2,e2,E1=e1

   if badpar(x1,[1,2,3,4,5],  1,caller='INTERP: (x1) ',npts=n1) then return
   if badpar(y1,[1,2,3,4,5],  1,caller='INTERP: (y1) ',npts=n2) then return
   if badpar(x2,[1,2,3,4,5],[0,1],caller='INTERP: (x2) ', $
                                npts=nout,type=outtype) then return

   if badpar(e1,[0,1,2,3,4,5],  1,caller='INTERP: (e1) ', $
                                  npts=n3,type=e1type) then return

   calc_err = e1type ne 0

   IF calc_err THEN BEGIN

      if (min([n1,n2,n3]) ne max([n1,n2,n3])) then begin
         print,'ITERP: x1, y1, and e1 must all be the same length!'
         return
      endif

      if n1 eq 1 then begin
         print,'INTERP: x1, y1, and e1 must have more than one element!'
         return
      endif

   ENDIF ELSE BEGIN

      if n1 ne n2 then begin
         print,'ITERP: x1 and y1 must be the same length!'
         return
      endif

      if n1 eq 1 then begin
         print,'INTERP: x1 and y1 must have more than one element!'
         return
      endif

   ENDELSE

   IF outtype eq 5 THEN BEGIN
      y2 = dblarr(nout,/nozero)
      IF calc_err THEN e2 = dblarr(nout,/nozero)
   ENDIF ELSE BEGIN
      y2 = fltarr(nout,/nozero)
      IF calc_err THEN e2 = fltarr(nout,/nozero)
   ENDELSE

   l = lonarr(nout,/nozero)
   r = lonarr(nout,/nozero)
   for i=0L,nout-1 do begin
      l[i] = max(where(x1 le x2[i]))
      r[i] = min(where(x1 gt x2[i]))
   endfor

   z=where(l eq -1,count)
   IF count ne 0 THEN l[z]=0
   
   z=where(r eq -1,count)
   IF count ne 0 THEN r[z]=n1-1

   z=where(r eq l and r eq 0,count)
   IF count ne 0 THEN r[z]=1

   z=where(r eq l,count)
   IF count ne 0 THEN l[z]=r[z]-1

   xc=(x2-x1[l])/(x1[r]-x1[l])
   y2 = y1[l] + (y1[r]-y1[l])*xc

   IF calc_err THEN BEGIN
      e2a = sqrt(e1[l]^2 + (e1[r]^2+e1[l]^2)*xc^2)
      e2b = sqrt(e1[r]^2 + (e1[r]^2+e1[l]^2)*(1-xc)^2)
      e2  = e2a < e2b
   ENDIF

end
