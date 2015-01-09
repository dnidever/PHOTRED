pro rotsph,lon,lat,clon,clat,nlon,nlat,anode=anode,stp=stp,reverse=reverse,$
           original=original

;+
;
; ROTSPH
;
; This rotates a spherical coordinate system to a new
; pole
;
; I got the equations for this from the paper
; Calabretta et al. 2002, A&A, 395, 1077
; Equation 5.
;
; Also, see rotate_lb.pro that uses a matrix method
; and for which you can specify the equator you'd like.
; rotsph.pro is faster than rotate_lb.pro
; By default, the origin is the point where the two equators
; cross (unless =ANODE is set).
; This should give you the same result (within ~1E-10")
; rotate_lb,lon,lat,[clon,clat],[clon+90,0.],nlon,nlat
;
; INPUTS:
;  lon       Array of longitudes to be rotated
;  lat       Array of latitudes to be rotated
;  clon      Longitude of the new NORTH POLE in the old coordinate system
;  clat      Latitude of the new NORTH POLE in the old coordinate system
;  =anode    The "Ascending Node" which is the longitude (in the new
;              system) of the first point where the old equator cross
;              the new one.  This sets the zero-point of the new
;              longitude.  By default the zero-point of the new
;              coordinate system is the point where the two equators
;              cross.
;  /original Set the new longitude zero-point to be clon (if clat>0)
;              and clon+180 (if clat<0).  This is the way it was
;              originally done.  DON'T USE WITH "ANODE"
;  /stp      Stop at the end of the program
;  /reverse  The reverse operation.  In that case (nlon,nlat) should be input
;            as (lon,lat). E.g.
;
;            IDL>rotsph,ra,dec,cra,cdec,nlon,nlat
;            IDL>rotsph,nlon,nlat,cra,cdec,nra,ndec,/reverse
;            
;            (ra,dec) and (nra,ndec) should be identical to 1E-10.
;
; OUTPUTS:
;  nlon  Array of rotated longitudes
;  nlat  Array of rotated latitudes
;
; By D.Nidever  Jan.2007
;-

; Not enough inputs
if n_elements(lon) eq 0 or n_elements(lat) eq 0 or n_elements(clon) eq 0 or $
   n_elements(clat) eq 0 then begin
  print,'Syntax - rotsph,lon,lat,clon,clat,nlon,nlat,reverse=reverse,stp=stp'
  return
endif

radeg = 180.0d0/!dpi

; I think the clon,clat in the equation is for the SOUTH pole
;clon2 = clon+180.0d0
;clon2 = clon2 mod 360.0d0 
;clat2 = clat

alphap = double(clon[0])/radeg
deltap = double(clat[0])/radeg
phip = 90.0d0/radeg   ;180.0d0/radeg
if keyword_set(original) then phip = 180.0d0/radeg   ; original way
thetap = 90.0d0/radeg

; By default the origin of the new coordinate system is the point
; where the two equators cross for the first time
;  Unless /original is set.  Then the zero-point is at clon
;   (if clat>0) and clon+180 (if clat<0)

; NORMAL
if not keyword_set(reverse) then begin

  alpha = double(lon)/radeg
  delta = double(lat)/radeg
  ;alphap = double(clon[0])/radeg
  ;deltap = double(clat[0])/radeg
  ;phip = 180.0d0/radeg

  ; arg(x,y) but atan(y,x)
  phi = phip + atan( -cos(delta)*sin(alpha-alphap), sin(delta)*cos(deltap)-$
                    cos(delta)*sin(deltap)*cos(alpha-alphap) )

  theta = asin( (-1)>(sin(delta)*sin(deltap)+cos(delta)*cos(deltap)*cos(alpha-alphap))<1 )

  ; Preparing the output
  nlon = phi*radeg
  nlat = theta*radeg

  ; Ascending Node
  ;  By default the origin of nlon is the point where the two equators
  ;  cross the first time
  if n_elements(anode) gt 0 then nlon = nlon + anode


; REVERSE
endif else begin

  phi = double(lon)/radeg  
  theta = double(lat)/radeg  
  ;alphap = double(clon[0])/radeg
  ;deltap = double(clat[0])/radeg
  ;phip = 180.0d0/radeg
  ;thetap = 90.0/radeg

  ; Ascending Node
  if n_elements(anode) gt 0 then phi = double(lon-anode)/radeg

  ; arg(x,y) but atan(y,x)
  alpha = alphap + atan( -cos(theta)*sin(phi-phip), sin(theta)*cos(deltap) - $
                         cos(theta)*sin(deltap)*cos(phi-phip))
  delta = asin( sin(theta)*sin(deltap) + cos(theta)*cos(deltap)*cos(phi-phip) )

  ; Preparing the output
  nlon = alpha*radeg
  nlat = delta*radeg

endelse

; Want everything less than 360.0
nlon = nlon mod 360.0

; Make negative points positive
bd = where(nlon lt 0.0,nbd)
if nbd gt 0 then nlon[bd]=nlon[bd]+360.0d0

if keyword_set(stp) then stop

end
