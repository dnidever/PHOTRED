;+
;
; NAME:
;  airmass
;
; PURPOSE:
;  Compute airmass for one or more times.
;
; DESCRIPTION:
;  This is should be a pretty good function for computing the air mass factor.
;  The default is to use the cosine based formula derived by David Tholen
;  but the older secant based formula from Hardie is still available.  The
;  zenith angle is corrected for refraction before using either formula (see
;  REFRAC).  The defaults on the atmospheric conditions are STP.  This function
;  should be quite good up to 5 airmasses.  This formula will work up to
;  a zenith angle of 80 degrees after which the computation id not done.
;
; CATEGORY:
;  Astronomy
;
; CALLING SEQUENCE:
;  am = airmass(jd,ra,dec,lat,lon)
;
; INPUTS:
;  jd  - Julian date (must be double precision to get nearest second).
;  ra  - Right ascension (of date) in degrees.
;  dec - Declination (of date) in degrees.
;  lat - Latitude of observatory in degrees.
;  lon - West longitude of observatory in degrees.
;
; OPTIONAL INPUT KEYWORD PARAMETERS:
;  /stp - Stop at the end of the program.
;
; KEYWORD OUTPUT PARAMETERS:
;  ALT - Optional return of the altitude for each airmass (in degrees)
;  AZI - Optional return of the azimuth (west from south, in degrees).
;  LHA - Optional return of the local hour angle (in degrees).
;  LST - Optional return of the local sidereal time (in degrees).
;
; OUTPUTS:
;  Return value is the airmass in single precision.
;
; RESTRICTIONS:
;  Any input may be a vector.  If more than one is a vector then the
;  lengths must match.  The return will have the same dimensions as
;  the input.
;
; PROCEDURS USED:
;  badpar.pro
;
; MODIFICATION HISTORY:
;  Written 1992 March, by Marc W. Buie, Lowell Observatory
;  94/05/05 - MWB, modified to split out the LST calculation, added LST and
;                LHA optional keyword outputs and UT keyword input.'
;  97/03/03 - MWB, added the Tholen airmass equation as the default.  Hardie
;                is still included as an option.  Also added refraction
;                to the calculation which added numerous inputs.
;  97/10/23 - MWB, added the AZI keyword.
;  2002/01/08, MWB, now allows 2-d inputs on jd,ra,dec
;  2005/10/0 - DLN, added /degree keyword
;  2008/05/06 - DLN, simplified
;
;
;-

;------------------------------------------------------------------

pro am_lsidtim,in_jd,lon,lst,UT=ut

   if badpar(in_jd,[5],[0,1,2],CALLER='lsidtim: (jd) ') then return
   if badpar(lon,[4,5],0,CALLER='lsidtim: (lon) ') then return
   if badpar(ut,[0,2,3,4,5],0,CALLER='lsidtim: (UT) ',DEFAULT=0) then return

   ; add time offset to what ever JD is provided (need not be at 0h.
   jd = in_jd + double(ut)/24.0d0

   ; Extract JD at 0h UT
   jd0 = long(jd+0.5d0) - 0.5d0

   ; compute auxillary quantity
   t = (jd0-2451545.0d0)/36525.0d0

   ; compute mean sidereal time at Greenwich at 0h UT
   gst = t*(8640184.812866d0 + $
         t*(    0.093104d0   - $
         t*     6.2d-6 ))
   gst = 1.753368559233d0 + gst/43200.0d0*!dpi

   hour = (jd-jd0)*2.0d0*!dpi * 1.00273790935d0

   st = gst + hour

   lst = st - lon

   lst = lst mod (2.0d0*!dpi)

end

;------------------------------------------------------------------

pro am_hangle,jd,ra,dec,lat,lon,lha,lst,UT=ut

   if badpar(jd,[5],[0,1,2],CALLER='hangle: (jd) ') then return
   if badpar(ra,[4,5],[0,1,2],CALLER='hangle: (ra) ') then return
   if badpar(dec,[4,5],[0,1,2],CALLER='hangle: (dec) ') then return
   if badpar(lat,[4,5],0,CALLER='hangle: (lat) ') then return
   if badpar(lon,[4,5],0,CALLER='hangle: (lon) ') then return
   if badpar(ut,[0,2,3,4,5],0,CALLER='hangle: (UT) ',DEFAULT=0) then return

   AM_LSIDTIM,jd,lon,lst,UT=ut

   lha = lst - ra
   z = where(lha gt !pi,count)
   if (count ne 0) then lha[z] = lha[z] - 2.0*!pi
   z = where(lha le -!pi,count)
   if (count ne 0) then lha[z] = lha[z] + 2.0*!pi

end

;------------------------------------------------------------------

function airmass,jd,ra,dec,lat,lon, $
            ALT=alt,LHA=lha,LST=lst,AZI=azi

   njd = n_elements(jd)
   nra = n_elements(ra)
   ndec = n_elements(dec)
   nlat = n_elements(lat)
   nlon = n_elements(lon)


   ; Not enough inputs
   if (njd eq 0 or nra eq 0 or ndec eq 0 or nlon eq 0 or nlat eq 0) then begin
     print,'Syntax - am = airmass,jd,ra,dec,lat,lon,ALT=alt,LHA=lha,LST=lst,AZI=azi)'
     return,-1
   endif

   ; Arrays must be of the same length
   if (njd ne nra) or (njd ne ndec) then begin
     print,'Arrays must be of the same length'
     return,-1
   endif

   ; Initialize airmass array
   am = replicate(-1.0,n_elements(jd))

   ; Checking data types
   if badpar(jd,[5],[0,1,2],CALLER='airmass: (jd) ') then return,am
   if badpar(ra,[4,5],[0,1,2],CALLER='airmass: (ra) ') then return,am
   if badpar(dec,[4,5],[0,1,2],CALLER='airmass: (dec) ') then return,am
   if badpar(lat,[4,5],0,CALLER='airmass: (lat) ') then return,am
   if badpar(lon,[4,5],0,CALLER='airmass: (lon) ') then return,am
   ;wave = 0.56
   ;pressure = 760.0
   ;temp = 0.0
   ;relhum = 0.0
   radeg = 180.0d0/!dpi

   ; Convert from degrees to radians
   ra2 = ra/radeg
   dec2 = dec/radeg
   lat2 = lat/radeg
   lon2 = lon/radeg

   ; Get hour angle
   AM_HANGLE,jd,ra2,dec2,lat2,lon2,lha,lst

   ; Get altitude and azimuth
   alt = asin( sin(lat2)*sin(dec2) + cos(lat2)*cos(dec2) *cos(lha) )

   azi = atan( sin(lha), cos(lha)*sin(lat2) - tan(dec2)*cos(lat2) )

   ; Zenith distance
   zenith = !dpi/2.0-alt
   z = where(zenith le 1.521,count)   ; <87.15 deg

   ; Some good ones
   if (count ne 0) then begin
      ;zenith[z] = REFRAC(zenith[z],wave,pressure,temp,relhum)      
      ;n = airindex(wave,pressure,temp,relhum)
      n = 1.00029   ; for default wave, pressure, temp, relhum
      zenith[z] = asin(sin(zenith[z])/n)

      cz = cos(zenith[z])
      am[z] = sqrt(235225.0*cz*cz + 970.0 + 1.0) - 485*cz

   endif

   ; Convert to degrees
   alt = alt*radeg
   azi = azi*radeg
   lha = lha*radeg
   lst = lst*radeg

   if keyword_set(stp) then stop

   return,float(am)

end
