PRO VCONV,vhel,glon,glat,vlsr,vgsr,vcirc=vcirc,bovy=bovy,schoenrich2010=schoenrich2010,lsrk=lsrk

;+
;
; VCONV
;
; Given an input array called vhel with values of heliocentric
; radial velocities and input arrays of the same length with
; the Galactic coordinates of the object (gl,gb),
; this code calculates Vlsr and Vgsr.
;
;  code assumes gl & gb given in degrees.
;
; INPUTS:
;  vhel   Array of heliocentric velocities (in km/s).
;  glon   Galactic longitude (in degrees)
;  glat   Galactic latitude (in degrees)
;  =vcirc Circular rotation velocity at the sun (220 km/s by default)
;
; OUTPUTS:
;  vlsr   Array of Local Standard of Rest velocities
;  vgsr   Array of Galactic Standard of Rest velocities
;
;
; USAGE:
;  IDL>vconv,vhelio,glon,glat,vlsr,vgsr
;
; By ??
;  updated by D.NIdever
;-

; Are there enough inputs
if n_params() lt 3 then begin
  print,'Syntax - vconv,vhel,glon,glat,vlsr,vgsr,vcirc=vcirc'
  return
endif

if not keyword_set(vcirc) then vcirc=220.0d

gl = glon*(!dpi/180.0d0)
gb = glat*(!dpi/180.0d0)

cgl = cos(gl) & sgl=sin(gl)
cgb = cos(gb) & sgb=sin(gb)

;; Bovy et al. (2012,2014)
if keyword_set(bovy) then begin
   ; Bovy et al. (2012,2014)
   ; VR = -10.5 +/- 1.0 km/s
   ; Vphi-Vc = 23.1 +3.6-0.5 km/s 
   ; Vz not measured, use D&B98
   ; The sun's motion
   ; Vxsun = 10.5
   ; Vysun = 242.0
   ; Vzsun = 7.0
   ; Vcirc = 218 km/s
   vlsr = vhel + ((10.50d0*cgb*cgl)+(23.1d0*cgb*sgl)+(7.00d0*sgb))
endif

;; Schoenrich et al. (2010)
if keyword_set(schoenrich2010) then begin
  ;; 11.1, 12.24, 7.25
   vlsr = vhel + ((11.10d0*cgb*cgl)+(12.24d0*cgb*sgl)+(7.25d0*sgb))
endif

;; Kinematic Local Standard of Rest
;;  defined from an average of the velocities of stars in the Solar neighborhood (Delhaye 1965; Gordon 1976)
;;  Solar motion = 20.0 km/sec towards (18h +30°) at epoch 1900.0.
if keyword_set(lsrk) then begin

  ;; This is from GBTIDL pro/toolbox/chdoppler.pro
  ;;glactc,ra,dec,2000.0,glon,glat,2
  ;;rasource=ra*15.d*!dtor
  ;;decsource=dec*!dtor
  ;;
  ;;nin = n_elements(glon)
  ;;xxsource = dblarr( 3, nin)
  ;;xxsource[0, *] = cos(decsource) * cos(rasource)
  ;;xxsource[1, *] = cos(decsource) * sin(rasource)
  ;;xxsource[2, *] = sin(decsource)
  ;;pvlsrk = dblarr( nin)
  ;;
  ;;;-----------------------LSRK SECTION-------------------------
  ;;;THE STANDARD LSRK IS DEFINED AS FOLLOWS: THE SUN MOVES AT 20.0 KM/S
  ;;;TOWARD RA=18.0H, DEC=30.0 DEG IN 1900 EPOCH COORDS
  ;;;using PRECESS, this works out to ra=18.063955 dec=30.004661 in 2000
  ;;;coords.
  ;;ralsrk_rad= 2.d*!pi*18.d/24.d
  ;;declsrk_rad= !dtor*30.d
  ;;precess, ralsrk_rad, declsrk_rad, 1900.d, 2000.d,/radian
  ;;
  ;;;FIND THE COMPONENTS OF THE VELOCITY OF THE SUN WRT THE LSRK FRAME
  ;;xxlsrk = dblarr( 3, nin)
  ;;xxlsrk[ 0, *] = cos(declsrk_rad) * cos(ralsrk_rad)
  ;;xxlsrk[ 1, *] = cos(declsrk_rad) * sin(ralsrk_rad)
  ;;xxlsrk[ 2, *] = sin(declsrk_rad)
  ;;vvlsrk = 20.d*xxlsrk
  ;;
  ;;;PROJECTED VELOCITY OF THE SUN WRT LSRK TO THE SOURCE
  ;;for nr=0, nin-1 do pvlsrk[ nr]=total(vvlsrk*xxsource[ *, nr])

  ralsrk = 18.d * 15.d
  declsrk = 30.d
  precess, ralsrk, declsrk, 1900.d, 2000.d
  glactc,ralsrk,declsrk,2000.0,gllsrk,gblsrk,1,/deg
  ucomp = 20.d * cos(gblsrk/!radeg) * cos(gllsrk/!radeg)
  vcomp = 20.d * cos(gblsrk/!radeg) * sin(gllsrk/!radeg)
  wcomp = 20.d * sin(gblsrk/!radeg)
  ;;  10.270606       15.317361       16.611399
  vlsr = vhel + ((ucomp*cgb*cgl)+(vcomp*cgb*sgl)+(wcomp*sgb))

  ;; Just take the dot-product of the two vectors,
  ;;  i.e. cos(angular distance between them)
  ;;angdist = sphdist(gllsrk,gblsrk,glon,glat,/deg)  
  ;;voff = 20.d * cos(angdist/!radeg)
   
endif

;; Dehnen & Binney 1998 values
if not keyword_set(bovy) and not keyword_set(schoenrich2010) and not keyword_set(lsrk) then begin
  ;  This equation takes the solar motion w.r.t. the LSR as
  ;  (9,11,6) km/sec (Mihalas & Binney)
  ;   Using updated values from Dehnen & Binney 1998, MNRAS, 298, 387
  ;   U = 10.00 +/- 0.36 km/s
  ;   V = 5.25 +/- 0.62 km/s
  ;   W = 7.17 +/- 0.38 km/s
  ;   This is in a right-handed system

   vlsr = vhel + ((10.00d0*cgb*cgl)+(5.25d0*cgb*sgl)+(7.17d0*sgb))
endif
;vlsr=vhel+((9.0d0*cgb*cgl)+(11.0d0*cgb*sgl)+(6.0d0*sgb))

; Schoenrich+10
; U = 11.1
; V = 12.24
; W = 7.25

;  This equation takes the rotation velocity of the LSR to
;  be 220 km/sec

vgsr = vlsr + (vcirc*cgb*sgl)

return

end

