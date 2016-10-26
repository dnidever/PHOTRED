pro curpdiff,sph=sph,arcsec=arcsec

; This program allows you to find the distance between two points.

print,'FIRST CLICK'
cursor,x1,y1
wait,0.2

print,'SECOND CLICK'
cursor,x2,y2

; Spherical ra/dec coordinates (or similar)
;  use cos(dec) correction 
if keyword_set(sph) then begin
  mndec = mean([y1,y2])
  cosdec = cos(mndec/!radeg)
  ; Convert to arc seconds
  if keyword_set(arcsec) then begin
    mfac = 3600.0d0
    units = 'arcsec'
  endif else units='deg'
endif else begin
  cosdec = 1.0d0
  mfac = 1.0d0
  units = 'pixels'
endelse

print,'Distance = ',stringize(sqrt(((x1-x2)*cosdec)^2. + (y1-y2)^2.)*mfac,ndec=4),' ',units
print,'Delta X = ',stringize((x2-x1)*mfac*cosdec,ndec=4),' ',units
print,'Delta Y = ',stringize((y2-y1)*mfac,ndec=4),' ',units
print,'Angle = ',stringize(atan(y2-y1,x2-x1)*!radeg,ndec=4),' (CCW from Right)'
lineq,x1,y1,x2,y2

;stop

end
