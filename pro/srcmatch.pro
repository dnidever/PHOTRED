;+
;
; SRCMATCH
;
; This works like srcor except that it breaks the chunks up into different
; domains so it will go faster.
;
; INPUTS:
;  xarr1     The array of X (RA,LON) values for the 1st set
;  yarr1     The array of Y (DEC,LAT) values for the 1st set
;  xarr2     The array of X (RA,LON) values for the 2nd set
;  yarr2     The array of Y (DEC,LAT) values for the 2nd set
;  dcr       The matching radius.  In arcseconds if /sph set.
;  /sph      Spherical coordinates input.  Arrays need to be in  *DEGREES*
;             and dcr in *ARCSECONDS*.
;  option=   Different options:
;            OPTION=0  Closest match from list2 is found for each element
;                       of list1 (within DCR).
;            OPTION=1  One-to-one mapping.  The default.
;            OPTION=2  Same as OPTION=1 but DCR is ignored.
;  /stp      Stop at the end
;  domains=  Specify the number of domains desired.  Otherwise it
;              will return the number of domains used.
;  /usehist  Use the HISTOGRAM_ND programs MATCH_SPH/MATCH_2D to do
;              faster matching.  This is now the default.
;
; OUTPUTS:
;  ind1      Index of matched objects for 1st set
;  ind2      Index of matched objects for 2nd set
;  =count    The number of matches.  This is set to -1 if there was an error.
;
; USAGE:
;  IDL>srcmatch,xorig1,yorig1,xorig2,yorig2,dcr,ind1,ind2,count=count,sph=sph
;
; By D.Nidever   April 2007
;-


PRO srcor2,x1in,y1in,x2in,y2in,dcr,ind1,ind2,option=option,magnitude=magnitude,$
   spherical=spherical,silent=silent
;+
; NAME:
;       SRCOR
; PURPOSE:
;       Correlate the source positions found on two lists.
; CALLING SEQUENCE:
;       srcor,x1in,ylin,x2in,y2in,dcr,ind1,ind2
; INPUTS:
;       x1in,y1in - First set of x and y coordinates.  The program
;                   marches through this list element by element,
;                   looking in list 2 for the closest match.  So, the program
;                   will run faster if this is the shorter of the two lists.
;                   Unless you use the option or magnitude keyword, there is
;                   nothing to guarantee unique matches.  
;       x2in,y2in - Second set of x and y coordinates.  This list is
;                   searched in its entirety every time one element of list 1
;                   is processed.
;       dcr - Critical radius outside which correlations are rejected;
;             but see 'option' below.
; OPTIONAL KEYWORD INPUT:
;       option - Changes behavior of program and description of output
;                lists slightly, as follows: 
;       OPTION=0 or left out
;             Same as older versions of SRCOR.  The closest match from list2
;             is found for each element of list 1, but if the distance is
;             greater than DCR, the match is thrown out.  Thus the index
;             of that element within list 1 will not appear in the IND1 output
;             array.
;       OPTION=1
;             Forces the output mapping to be one-to-one.  OPTION=0 results,
;             in general, in a many-to-one mapping from list 1 to list 2.
;             Under OPTION=1, a further processing step is performed to
;             keep only the minimum-distance match, whenever an entry from
;             list 1 appears more than once in the initial mapping.
;       OPTION=2
;             Same as OPTION=1, except the critical distance parameter DCR
;             is ignored.  I.e., the closest object is retrieved from list 2
;             for each object in list 1 WITHOUT a critical-radius criterion,
;             then the clean-up of duplicates is done as under OPTION=1.
;       magnitude
;             An array of stellar magnitudes corresponding to x1in and y1in.  
;             If this is supplied, then the brightest star from list 1
;             within the selected distance of the star in list 2 is taken.
;             The option keyword is ignored in this case.
;       spherical
;             If SPHERICAL=1, it is assumed that the input arrays are in
;             celestial coordinates (RA and Dec), with x1in and x2in in
;             decimal hours and y1in and y2in in decimal degrees.  If
;             SPHERICAL=2 then it is assumed that the input arrays are in
;             longitude and latitude with x1in,x2in,y1in,y2in in decimal
;             degrees.  In both cases, the critial radius dcr is in
;             *arcseconds*.  Calculations of spherical distances are made
;             with the gcirc program.
; OUTPUTS:
;       ind1 - index of matched stars in first list
;       ind2 - index of matched stars in second list
; COMMON BLOCKS:
;       none
; SIDE EFFECTS:
;       none
; METHOD:
;       See under keyword LEVEL above.
; REVISON HISTORY:
;       Adapted from UIT procedure  J.Wm.Parker, SwRI 29 July 1997
;       Converted to IDL V5.0   W. Landsman   September 1997
;       
;-
;
 ON_Error,2   ; Return if error (incl. non-info message)

;;;
;   If not enough parameters, then print out the syntax.
;
IF N_params() lt 7 THEN BEGIN
  print,'SRCOR calling sequence: '
  print,'srcor,x1in,y1in,x2in,y2in,dcr,ind1,ind2 [,option={0, 1, or 2}] $'
  print,'      [,magnitude=mag_list_1, spherical={1 or 2}]'
  RETURN
ENDIF

;;;
;   Keywords.
;
IF not keyword_set(option) THEN option=0
;message,/info,'Option code = '+strtrim(option,2)
;IF (option lt 0) or (option gt 2) THEN MESSAGE,'Invalid option code.'

SphereFlag = keyword_set(Spherical)

;;;
;   Store the input variables into internal arrays that we can manipulate and
; modify.
;
x1=float(x1in)
y1=float(y1in)
x2=float(x2in)
y2=float(y2in)

;;;
;   If the Spherical keyword is set, then convert the input values (degrees
; and maybe hours) into radians, so GCIRC doesn't have to make this calculation
; each time it is called in the FOR loop.  Also convert the critical radius
; (which is in arcsec, so convert by 3600.) to radians
;
if SphereFlag then begin
   dcr2 = dcr
   if (Spherical eq 1) then XScale = 15.0 else XScale = 1.0 
   d2r  = !DPI/180.0d0
   x1 = x1 * (XScale * d2r)
   y1 = y1 * d2r
   x2 = x2 * (XScale * d2r)
   y2 = y2 * d2r
   dcr2 = dcr2 * (d2r / 3600.)
endif else dcr2=dcr^2


;;;
;   Set up some other variables.
;
n1=n_elements(x1) ;& message,/info,strtrim(n1,2)+' sources in list 1'
n2=n_elements(x2) ;& message,/info,strtrim(n2,2)+' sources in list 2'
nmch=0
ind1=-1L & ind2=-1L

;;;
;   The main loop.  Step through each index of list 1, look for matches in 2.
;
FOR i=0L,n1-1 DO BEGIN
   xx=x1[i] & yy=y1[i] 
   if SphereFlag then gcirc,0,xx,yy,x2,y2,d2 else d2=(xx-x2)^2+(yy-y2)^2
   dmch=min(d2,m)
   IF (option eq 2) or (dmch le dcr2) THEN BEGIN
      nmch=nmch+1
      IF nmch eq 1 THEN BEGIN 
         ind1=long(i) 
         ind2=long(m)
      ENDIF ELSE BEGIN
         ind1=[ind1,i]
         ind2=[ind2,m]
      ENDELSE
   ENDIF
ENDFOR
;message,/info,strtrim(nmch,2)+' matches found.'

;;;
;   Modify the matches depending on input options.
;
use_mag = (n_elements(magnitude) ge 1)
IF (option eq 0) and (not use_mag) THEN RETURN
IF use_mag THEN BEGIN
   ;message,/info,'Cleaning up output list using magnitudes.'
ENDIF ELSE BEGIN
   ;IF option eq 1 then message,/info,'Cleaning up output list (option = 1).'
   ;IF option eq 2 then message,/info,'Cleaning up output list (option = 2).'
ENDELSE

FOR i=0L,max(ind2) DO BEGIN
   csave = n_elements(ind2)
   ww = where(ind2 eq i,count) ; All but one of the list in WW must
                               ; eventually be removed.
   IF count gt 1 THEN BEGIN
      IF use_mag THEN BEGIN
         dummy = min(magnitude[ind1[ww]],m)
      ENDIF ELSE BEGIN
         xx=x2[i] & yy=y2[i]
         if SphereFlag then gcirc,0,xx,yy,x1[ind1[ww]],y1[ind1[ww]],d2 else $
                            d2=(xx-x1[ind1[ww]])^2+(yy-y1[ind1[ww]])^2
         ;IF n_elements(d2) ne count THEN MESSAGE,'Logic error 1'
         dummy = min(d2,m)
      ENDELSE
      remove,m,ww              ; Delete the minimum element
                               ; from the deletion list itself.

      remove,ww,ind1,ind2      ; Now delete the deletion list from
                               ; the original index arrays.
      ;IF n_elements(ind2) ne (csave-count+1) THEN MESSAGE,'Logic error 2'
      ;IF n_elements(ind1) ne (csave-count+1) THEN MESSAGE,'Logic error 3'
      ;IF n_elements(ind2) ne n_elements(ind1) THEN MESSAGE,'Logic error 4'
   ENDIF
ENDFOR

;message,/info,strtrim(n_elements(ind1),2)+' left in list 1'
;message,/info,strtrim(n_elements(ind2),2)+' left in list 2'

;
RETURN
end

;------------------------------------------------------------------------

pro factors,x,f

; This returns the array of factors of X
; This only works if X isn't too large

; Too large
if x gt 1d8 then begin
  print,'Number is too large'
  f=-1
  return
endif

; From 1 to x
;y = dindgen(x-2)+2
y = dindgen(x)+1
gd = where(x/y mod 1 eq 0,ngd)
f = y[gd]

end


;------------------------------------------------------------------------

pro srcmatch,xorig1,yorig1,xorig2,yorig2,dcr,ind1,ind2,sph=sph,stp=stp,$
             domains=domains,count=count,usehist=usehist

;+
;
; SRCMATCH
;
; This works like srcor except that it breaks the chunks up into different
; domains so it will go faster.
;
; INPUTS:
;  xarr1     The array of X (RA,LON) values for the 1st set
;  yarr1     The array of Y (DEC,LAT) values for the 1st set
;  xarr2     The array of X (RA,LON) values for the 2nd set
;  yarr2     The array of Y (DEC,LAT) values for the 2nd set
;  dcr       The matching radius.  In arcseconds if /sph set.
;  /sph      Spherical coordinates input.  Arrays need to be in  *DEGREES*
;             and dcr in *ARCSECONDS*.
;  /option   Different options:
;            OPTION=0  Closest match from list2 is found for each element
;                       of list1 (within DCR).
;            OPTION=1  One-to-one mapping.  The default.
;            OPTION=2  Same as OPTION=1 but DCR is ignored.
;  /stp      Stop at the end
;  domains=  Specify the number of domains desired.  Otherwise it
;              will return the number of domains used.
;  /usehist  Use the HISTOGRAM_ND programs MATCH_SPH/MATCH_2D to do
;              faster matching.  This is now the default.
;
; OUTPUTS:
;  ind1      Index of matched objects for 1st set
;  ind2      Index of matched objects for 2nd set
;  =count    The number of matches.  This is set to -1 if there was an error.
;
; USAGE:
;  IDL>srcmatch,xorig1,yorig1,xorig2,yorig2,dcr,ind1,ind2,count=count
;
; By D.Nidever   April 2007
;-

count=0

; Checking arguments
nx1 = n_elements(xorig1)
ny1 = n_elements(yorig1)
nx2 = n_elements(xorig2)
ny2 = n_elements(yorig2)

; Not enough inputs
if nx1 eq 0 or ny1 eq 0 or nx1 eq 0 or ny1 eq 0 then begin
  print,'Syntax - srcmatch,xarr1,yarr1,xarr2,yarr2,dcr,ind1,ind2,sph=sph,count=count,stp=stp'
  return
endif

; Arrays not of the same size
if (nx1 ne ny1) or (nx2 ne ny2) then begin
  print,'xarr1/2 and yarr1/2 not of the same size'
endif

; Error Handling
;------------------
; Establish error handler. When errors occur, the index of the  
; error is returned in the variable Error_status:  
CATCH, Error_status 

;This statement begins the error handler:  
if (Error_status ne 0) then begin 
   print,'SRCMATCH ERROR: ', !ERROR_STATE.MSG  
   ind1 = -1
   ind2 = -1
   count = -1
   CATCH, /CANCEL 
   return
endif


; Spherical, convert matching radius from arcseconds to degrees
dcr2 = dcr
if keyword_set(sph) then begin
  dcr2 = dcr2/3600.0d0

  xdcr = dcr2/cos(max([yorig1,yorig2])/!radeg)  ; in RA
  ydcr = dcr2
endif else begin
  xdcr = dcr2
  ydcr = dcr2
endelse


; --- Use the HISTOGRAM_ND matching programs instead ---
if n_elements(usehist) eq 0 then begin
  if running_gdl() eq 1 then usehist=0 else usehist=1  ; default
endif
if keyword_set(usehist) then begin
  ; Matching options
  ;  if option=1 or 2 then one_to_one=0, but default is to use one-to-one
  if n_elements(option) eq 1 then one_to_one=1-(option ge 1) else one_to_one=1
 
  ; Spherical
  if keyword_set(sph) then begin
    Result = MATCH_SPH([xorig1],[yorig1],[xorig2],[yorig2],dcr2,one_to_one=one_to_one)
  endif else begin  ; Cartesian
    p1 = dblarr(nx1,2) & p1[*,0]=xorig1 & p1[*,1]=yorig1
    p2 = dblarr(nx2,2) & p2[*,0]=xorig2 & p2[*,1]=yorig2
    Result = MATCH_ND(p1,p2,dcr2,one_to_one=one_to_one)
    ;Result = MATCH_2D([xorig1],[yorig1],[xorig2],[yorig2],dcr2) ; always one to one
  endelse

  ; Get the arrays and matches
  ind1 = where(result gt -1,count)
  if count gt 0 then ind2=result[ind1]

  ; Sometimes there is one index value that is out of range
  if count gt 0 then begin
    bd = where(ind1 gt (nx1-1) or ind2 gt (nx2-1),nbd,ncomp=ngd)
    ; Some bad ones to remove
    if nbd gt 0 then begin
      ; No good matches left
      if ngd eq 0 then begin
        ind1 = -1
        ind2 = -1
        count = 0
        return
      endif
      REMOVE,bd,ind1,ind2
    endif
  endif

  return
endif



  ; Ranges
  ; The overlapping range
mmx = [ min(xorig1) > min(xorig2), max(xorig1) < max(xorig2) ]
mmx = [ mmx[0]-1.5*xdcr, mmx[1]+1.5*xdcr ]
mmy = [ min(yorig1) > min(yorig2), max(yorig1) < max(yorig2) ]
mmy = [ mmy[0]-1.5*ydcr, mmy[1]+1.5*ydcr ]
;mmx = minmax([xarr1,xarr2])
;mmy = minmax([yarr1,yarr2])
minx = mmx[0]
miny = mmy[0]
xr = range(mmx)
yr = range(mmy)

; Making local copies
ind1orig = where(xorig1 ge mmx[0] and xorig1 le mmx[1] and yorig1 ge mmy[0] and yorig1 le mmy[1],nind1orig)
ind2orig = where(xorig2 ge mmx[0] and xorig2 le mmx[1] and yorig2 ge mmy[0] and yorig2 le mmy[1],nind2orig)
; Overlap
if (nind1orig gt 0 and nind2orig gt 0) then begin
  xarr1 = xorig1[ind1orig]
  yarr1 = yorig1[ind1orig]
  xarr2 = xorig2[ind2orig]
  yarr2 = yorig2[ind2orig]
endif else begin
  ind1 = -1
  ind2 = -1
  return
endelse


; How many domains do we want?
n = max([nx1,nx2])
;m = 10.0^(0.59875*alog10(n)-0.59322)
m = 10.0^(0.467*alog10(n)-0.471)
m = (round(m)>1.0)<1e4 ;100.0

; Round to nearest even integer, it looks like evens are faster
if m mod 2 ne 0 then m=m+1

if keyword_set(domains) then m=domains>1

;print,'Ndomains=',strtrim(m,2)

; Initializing the indices
undefine,ind1,ind2
;ind1 = fltarr(10.*n)  ; need extra in case of doubles
;ind2 = fltarr(10.*n)
;ngd = 0
mcount = 0

; Break up the domains into Mx, My
IF (m gt 1) then begin

  vol = xr*yr
  if keyword_set(sph) then vol=xr*yr*cos(mean(mmy)/!radeg)
  mvol = vol/m
  s = sqrt(mvol)

  ; Find the factors
  factors,m,f
  nf = n_elements(f)
  rf = reverse(f)

  ; Try all possibilities
  ; Which size is closest to the "ideal" size
  ddy = yr/f
  bestind = first_el(minloc(abs(s-float(ddy))))

  my = f[bestind]
  mx = rf[bestind]

  dx = xr/mx
  dy = yr/my

  ; Starting the domain looping

  ; Y domain loop
  for i=0,my-1. do begin

    ; Upper and lower Y limits for this domain
    y0 = miny+i*dy
    y1 = y0+dy + dcr2

    gdy1 = where(yarr1 ge y0 and yarr1 le y1,ngdy1)
    gdy2 = where(yarr2 ge y0 and yarr2 le y1,ngdy2)

    ; X domain loop
    for j=0,mx-1. do begin

      ; Spherical, correct for cos(dec)
      dcr3 = dcr2
      if keyword_set(spherical) then begin
        dcr3 = dcr2/cos(max([y0,y1])/!radeg)
      endif

      ; Upper and lower X limits for this domain
      x0 = minx+j*dx
      x1 = x0+dx + dcr3

      ; We've got some from arr1 in this domain
      gd1 = -1 & ngd1=0     ; bad by default
      if (ngdy1 gt 0) then begin

        gdx1 = where(xarr1[gdy1] ge x0 and xarr1[gdy1] le x1,ngdx1)

        if (ngdx1 gt 0) then begin
          gd1 = gdy1[gdx1]
          ngd1 = ngdx1
        endif
      endif

      ; We've got some from arr2 in this domain
      gd2 = -1 & ngd2=0   ; bad by default
      if (ngdy2 gt 0) then begin

        gdx2 = where(xarr2[gdy2] ge x0 and xarr2[gdy2] le x1,ngdx2)

        if (ngdx2 gt 0) then begin
          gd2 = gdy2[gdx2]
          ngd2 = ngdx2
        endif
      endif
      ;gd1 = where(xarr1 ge x0 and xarr1 le x1 and yarr1 ge y0 and yarr1 le y1,ngd1)
      ;gd2 = where(xarr2 ge x0 and xarr2 le x1 and yarr2 ge y0 and yarr2 le y1,ngd2)

      ;plot,xarr1,yarr1,ps=3
      ;oplot,xarr2,yarr2,ps=3,co=250
      ;if ngd1 gt 0 then oplot,xarr1[gd1],yarr1[gd1],ps=3,co=150
      ;if ngd2 gt 0 then oplot,xarr2[gd2],yarr2[gd2],ps=3,co=200

      ; There are some from both sets in this domain
      if (ngd1 gt 0) and (ngd2 gt 0) then begin
        ; Matching them
        if keyword_set(sph) then spherical=2   ; in degrees
        if not keyword_set(option) then opt=1 else opt=option
        SRCOR2,xarr1[gd1],yarr1[gd1],xarr2[gd2],yarr2[gd2],dcr,match1,match2,sph=spherical,opt=opt

        nmatch1 = n_elements(match1)
        nmatch2 = n_elements(match2)

        ; Do we have any matches
        if nmatch1 eq 1 then gmatch1 = where(match1 ne -1,ngmatch1) else ngmatch1=nmatch1
        if nmatch2 eq 1 then gmatch2 = where(match2 ne -1,ngmatch2) else ngmatch2=nmatch2

        ; Add to final indices
        if ngmatch1 gt 0 and ngmatch2 gt 0 then begin
          push,ind1,gd1[match1]
          push,ind2,gd2[match2]
          ;ind1[ngd] = gd1[match1]
          ;ind2[ngd] = gd2[match2]
          ;ngd = ngd+nmatch1
        end

        ;stop
      end

      ;stop
      mcount++
      ;if mcount mod 5 eq 0 then print,mcount

    end ; for j

  end ; for i

  ;; Remove elements at the end
  ;if ngd gt 0 then begin
  ;  ind1 = ind1[0:ngd-1]
  ;  ind2 = ind2[0:ngd-1]
  ;endif else begin
  ;  ind1 = -1
  ;  ind2 = -1
  ;endelse

  ; If no matches then return -1
  if n_elements(ind1) eq 0 then ind1=-1
  if n_elements(ind2) eq 0 then ind2=-1


  ; Get the unique ones
  ; Only want each star once

  ; First checking ind1
  ui1 = uniq(ind1,sort(ind1))
  ind1 = ind1[ui1]
  ind2 = ind2[ui1]

  ; Now check ind2
  ui2 = uniq(ind2,sort(ind2))
  ind1 = ind1[ui2]
  ind2 = ind2[ui2]


; Only 1 domain
ENDIF ELSE BEGIN

  ; Matching them
  SRCOR2,xarr1,yarr1,xarr2,yarr2,dcr,ind1,ind2,sph=sph,opt=1

ENDELSE

; Put the indices in terms of the original arrays
if (n_elements(ind1) gt 1 or ind1[0] ne -1) then begin
  ind1 = ind1orig[ind1]
  ind2 = ind2orig[ind2]
endif

; How many matches
dum = where(ind1 ne -1,count)

;stop

if keyword_set(stp) then stop

end
