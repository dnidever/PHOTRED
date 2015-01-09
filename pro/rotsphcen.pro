pro rotsphcen,lon,lat,clon,clat,nlon,nlat,polar=polar,gnomic=gnomic,stp=stp,$
              reverse=reverse

;+
;
; ROTSPHCEN
;
; This is very similar to rotsph.pro except that the coordinates
; input are not for the north pole but for the new equator.
; Everything is in DEGREES.
;
; INPUTS:
;  lon       Array of longitudes to be rotated
;  lat       Array of latitudes to be rotated
;  clon      Longitude of the new EQUATOR in the old coordinate system
;  clat      Latitude of the new EQUATOR in the old coordinate system
;  /polar    Return polar coordinates (rad,phi) instead of LON/LAT.
;            phi starts at North.
;  /gnomic   Also do a gnomic (tangent plane) projection.
;  /reverse  The reverse operation.  In that case (nlon,nlat) should be input
;            as (lon,lat). E.g.
;
;            IDL>rotsphcen,ra,dec,cra,cdec,nlon,nlat
;            IDL>rotsphcen,nlon,nlat,cra,cdec,nra,ndec,/reverse
;            
;            (ra,dec) and (nra,ndec) should be identical to 1E-10.
;
;
; OUTPUTS:
;  nlon  Array of rotated longitudes.  If /polar then this is PHI
;        the polar angle (measured from N toward E).
;        
;  nlat  Array of rotated latitudes.  If /polar then this is RAD
;        the polar radial distance.
;
; By D.Nidever  Jan.2007
;-

; Not enough inputs
if n_elements(lon) eq 0 or n_elements(lat) eq 0 or n_elements(clon) eq 0 or $
   n_elements(clat) eq 0 then begin
  print,'Syntax - rotsphcen,lon,lat,clon,clat,nlon,nlat,polar=polar,gnomic=gnomic,stp=stp,'
  print,'                   reverse=reverse'
  return
endif

radeg = 180.0d0/!dpi

; NOT polar coordinates
IF not keyword_set(polar) and not keyword_set(gnomic) then begin

  ; Get coordinates for the north pole
  np_lon = double(clon[0])
  np_lat = double(clat[0])+90.0d0
  if np_lat gt 90.0d0 then begin
    np_lon = double(clon[0])+180.0d0
    np_lon = np_lon mod 360.0d0

    np_lat = 90.0d0-double(clat)
  endif

  ; Run rotsph.pro
  ; NORMAL
  if not keyword_set(reverse) then begin

    ROTSPH,lon,lat,np_lon,np_lat,nlon,nlat,/original

  ; REVERSE
  endif else begin

    ROTSPH,lon,lat,np_lon,np_lat,nlon,nlat,/reverse,/original

    ; need to flip them around by 180 deg b/c the zero-point
    ;  is set by the NP lon
    if (double(clat[0])+90.0d0 gt 90.0d0) then begin
      nlon = (nlon + 180.0d0) mod 360.0d0
      nlat = -nlat
    endif

  endelse

  ; Make the longitudes continuous
  nlon = (nlon+180.0d0) mod 360.0d0
  nlon = nlon-180.0d0


; POLAR or GNOMIC
ENDIF ELSE BEGIN

  ; Making polar coordinates

  ;-----------------
  ; NORMAL
  ;------------------
  if not keyword_set(reverse) then begin


    ;Run rotsph.pro and specify the center of the field (the origin) as the
    ;  the Npole
    ROTSPH,lon,lat,clon[0],clat[0],phi,theta,/original

    ; phi is now going clockwise and from South
    orig_phi = phi
    phi = -phi+180.0d0      ; counterclockwise
    phi = phi mod 360.0d0
    rad = 90.0d0-theta

    ; Making gnomic projection
    if keyword_set(gnomic) then begin

      ; Scale the radius
      rad = radeg * cos(theta/radeg)/sin(theta/radeg)

      ; Now convert from gnomic polar to X/Y
      ; phi is from N toward E
      ; x = R*sin(phi)
      ; y = R*cos(phi)
      nlon = rad*sin(phi/radeg)
      nlat = rad*cos(phi/radeg)

    endif

    ; Output polar coordinates
    if keyword_set(polar) then begin
      nlon = phi
      nlat = rad
    endif


  ;-----------------
  ; REVERSE
  ;-----------------
  endif else begin


    ; Polar
    if keyword_set(polar) then begin
      phi = lon
      rad = lat
      ;rad = 90.0d0-theta
      theta = 90.0d0-rad
    endif

    ; Removing gnomic projection
    if keyword_set(gnomic) then begin

      ; Now convert from X/Y to gnomic polar
      ; phi is from N toward E
      ; x = R*sin(phi)
      ; y = R*cos(phi)
      ;nlon = rad*sin(phi/radeg)
      ;nlat = rad*cos(phi/radeg)
      rad = sqrt(lon^2.0+lat^2.0)
      phi = radeg*atan(lon,lat)      ; in degrees

      ; Scale the radius
      ;rad = radeg * cos(theta/radeg)/sin(theta/radeg)
      theta = radeg*atan(radeg/rad)   ; in degrees

    endif

    ; phi is now going clockwise and from South
    ;phi = -phi+180.0d0      ; counterclockwise
    ;phi = phi mod 360.0d0
    ;rad = 90.0d0-theta
    phi = -phi+180.0d0       ; reverse phi
    

    ;Run rotsph.pro and specify the center of the field (the origin) as the
    ;  the Npole
    ROTSPH,phi,theta,clon[0],clat[0],nlon,nlat,/reverse,/original

  endelse

ENDELSE  ; /polar or /gnomic

if keyword_set(stp) then stop

end
