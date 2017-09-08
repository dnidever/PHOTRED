;+
;
; WCSFIT
;
; This program automatically finds the WCS for an input FITS file.
; The USNO-B1 or 2MASS-PSC catalogs are used to do the fitting.
;
; INPUTS
;  input      Input of FITS filename(s).
;  =up        What is UP in the image (N, S, E, or W).  If this is
;               not input then up='N' is used by default.  MATCHSTARS.PRO
;               should be able to figure out the proper orientation.
;  =left      What is LEFT in the image (N, S, E, or W).  If this is
;               not input then left='E' is used by default.  MATCHSTARS.PRO
;               should be able to figure out the proper orientation.
;  =pixscale  The pixel scale in arc seconds/pixel.  This can sometimes
;               be obtained from the image header.
;  =cenra     An estimate of the central RA of the image (in DEGREES).
;               This can sometimes be obtained from the image header.
;  =cendec    An estimate of the central DEC of the image.  This can
;               sometimes be obtained from the image header.
;  =inpcat    The catalog (IDL structure) of detected sources in the image.
;               The structure must contain the tags X and Y.
;  =refname   The name of the reference catalog to use.  The choices are
;               'USNO-B1', '2MASS-PSC', 'UCAC4' or 'GAIA'.  USNO-B1 is the default.
;  =inprefcat  The reference catalog (IDL structure), normally USNO-B1 or 2MASS-PSC,
;               obtained from QUERYVIZIER.PRO.  The structure must contain the
;               tags RAJ2000/DEJ2000 or RA/DEC.
;  =refmaglim  The reference magnitude limit to use.  The default is 16.5 for 2MASS-PSC
;                and 21.0 for USNO-B1.
;  =catmaglim  The magnitude limit to use for the catalog of detected sources.  The default
;                is to use all "good" sources.
;  =caterrlim  The error limit to use for the catalog of detected sources.  The default
;                is to use all sources with err<1.0 mag.
;  =crpix1    The X value of the reference pixel.  The default is the center of the image.
;  =crpix2    The Y value of the reference pixel.  The default is the center of the image.
;  /noupdate  Do NOT Update the FITS header with the fitted WCS.  The
;               default is to update.
;  =projection The type of WCS projection to use.  "TAN" is used by default.
;  =searchdist The area (in arcmin) to be searched.  This is the area
;                for which reference stars (USNO-B1) will be acquired.
;                The default is the greater of Xsize, Ysize and 10 arcmin.
;  =rmslim    The maximum RMS to allow before the fit is "bad" and
;                rejected.  The default is rmslim=1.0 arcsec
;  =inpfwhm   Use this input FWHM for this image.  The default is to
;                derive it with IMFWHM.PRO.
;  /stp       Stop at the end of the program.
;
; OUTPUTS
;  =outcat       The catalog of detected sources in the image.
;  =outrefcat    The reference catalog of sources.
;  =matchrefcat  Matched reference catalog used for the final fitting
;  =matchcat     Matched image catalog used for the final fitting
;  =error        The error message if there was one, else undefined
;
; USAGE
;  IDL>wcsfit,'ccd1001.fits',up='N',left='E',pixscale=0.6995,cenra=90.55,cendec=-75.00
;
; By D.Nidever  Feb.2008
;-

pro wcsfit_dummy
FORWARD_FUNCTION trans_coo, trans_coo_dev, wcsfit_devcoo
end


;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function wcsfit_devcoo,par,x=x,y=y,ra=ra,dec=dec,astr=astr,wcs=wcs

; Make temporary astrometry structure
temp = astr
if n_elements(par) gt 4 then begin
  temp.cd[0,*] = par[0:1]
  temp.cd[1,*] = par[2:3]
  temp.crval = par[4:5]
endif else begin
  scale = par[0]
  theta = par[1]
  temp.cd[0,0] = scale*cos(theta)
  temp.cd[0,1] = scale*sin(theta)
  temp.cd[1,0] = -scale*sin(theta)
  temp.cd[1,1] = scale*cos(theta)
  temp.crval = par[2:3]
endelse

;temp.crpix = par[4:5]    ; Fixed now
;temp.crval = par[6:7]


; TNX/TPV projection
if n_elements(wcs) gt 0 then begin

  ; TNX
  if tag_exist(wcs,'TNX1') then begin
    ; Convert X/Y to RA/DEC
    tempwcs = wcs
    tempwcs.ast = temp
    WCSTNX_XY2RD,x,y,tempwcs,ra2,dec2,/degree

  ; TPV
  endif else begin
    ; Convert X/Y to RA/DEC
    tempwcs = wcs
    tempwcs.ast = temp
    WCSTPV_XY2RD,x,y,tempwcs,ra2,dec2,/degree
  endelse

; "Normal" projections
endif else begin

  ; Convert X/Y to RA/DEC
  xy2ad,x,y,temp,ra2,dec2

endelse

diff = sphdist(ra,dec,ra2,dec2,/deg)*3600.0

;print,par
;stop

return,diff

end


;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


pro wcsfit_find,filename,cat,exten=exten,fwhm=fwhm,inpfwhm=inpfwhm,inpim=inpim,silent=silent,error=error,stp=stp

; This runs IDL FIND/APER


undefine,error,im,fwhm,cat

; Not enough inputs
nfilename = n_elements(filename)
if (nfilename eq 0) then begin
  print,'Syntax - wcsfit_find,filename,cat,error=error'
  error = 'Not enough inputs'
  return
endif

; Does the file exist
test = FILE_TEST(filename)
if (test eq 0) and n_elements(inpim) eq 0 then begin
  print,filename,' NOT FOUND'
  error = filenamea+' NOT FOUND'
  return
endif


; Load the image
if n_elements(inpim) eq 0 then begin
  FITS_READ,filename,im,head,exten=exten,/no_abort,message=message
  if message ne '' then begin
    error = 'Problem loading '+filename
    print,error
    return    
  endif
  im = float(im)
endif else im=float(inpim)

; Use FIND, APER
;-----------------

; Get the FWHM
;IMFWHM,filename,fwhm,im=im
if n_elements(inpfwhm) eq 0 then begin
  IMFWHM,filename,fwhm,ellip,silent=silent,exten=exten,im=im,head=head
endif else begin
  fwhm = inpfwhm
  ellip = 0.0
endelse
;if fwhm gt 50 then fwhm=10
if fwhm gt 50 then begin
  error = 'Bad FWHM'
  print,error
  undefine,cat
  return
endif

; Find the background
skymode = -999999.
skysig = -999999.
SKY,im,skymode,skysig,/silent,highbad=max(im)-5000
if skysig lt 0.0 then skysig = mad(im)
if skymode lt 0.0 then skymode = median(im)

; MASK OUT BAD PIXELS AND COLUMNS!!!!
;skymask = float(abs(im-skymode) lt 3.0*skysig)
;gkernel = [ [1,1,1,1,1], [1,1,1,1,1], [1,1,1,1,1] ]
;checkmask = CONVOL(skymask,gkernel,/center,/edge_truncate) < 1
;sm = smooth(im,5,/edge)
;im2 = im-sm
;kernel = [ [0,-1,0], [-1,4,-1], [0,-1,0]]
;lorentz = CONVOL(checkmask*im2,kernel,/center,/edge_truncate,missing=0.0)
;mask = float(abs(lorentz) gt 1000 and im gt skymode+20*skysig)

; Bad pixels are causing problems, they are picked up as sources by FIND
; and APER returns 99.99 for any source with a "bad" pixel in it's aperture
; Temporarily mask out the bad pixels, set to SKYMODE. 5/19/2015.
tempim = im
bdpix = where(im gt max(im)-100,nbdpix)
if nbdpix gt 0 then tempim[bdpix]=skymode

; Find the sources
; X/Y start at 0, while daophot coordinates start at 1.
roundlim = [-1.0,1.0]* (ellip*2.5 > 1.0)  ; larger cutoff for high-ellipicity
;roundlim = [-1.0,1.0]
sharplim = [0.2,1.0]
FIND,tempim,x,y,flux,sharp,round,4.0*skysig,fwhm,roundlim,sharplim,/silent
;FIND,im,x,y,flux,sharp,round,4.0*skysig,fwhm,[-1.0,1.0],[0.2,1.0],/silent

; NO stars found
if (n_elements(x) eq 0) then begin
  print,'NO stars found'
  error = 'NO stars found'
  undefine,cat
  return
endif

; Get aperture photometry
lo = skymode-20.0*skysig
hi = max(im)
APER,tempim,x,y,mags,errap,sky,skyerr,1.0,3.0*fwhm,[40,50],[lo,hi],/silent,/meanback
;APER,im,x,y,mags,errap,sky,skyerr,1.0,3.0*fwhm,[40,50],[lo,max(im)],/silent,/meanback

; Convert coordinates from IDL to DAOPHOT/IRAF/FITS format (0 indexed to 1 indexed)
x = x + 1.0
y = y + 1.0

; Put it all in a structure
nstars = n_elements(x)
dum = {id:0L,x:0.0d0,y:0.0d0,mag:0.0,err:0.0,sky:0.0,skysig:0.0,flux:0.0,sharp:0.0,round:0.0,round2:0.0}
;dum = {id:0,x:0.0,y:0.0,flux:0.0,sharp:0.0,round:0.0,mag:0.0,err:0.0,sky:0.0,skyerr:0.0}
cat = REPLICATE(dum,nstars)
cat.id = lindgen(nstars)+1
cat.x = x
cat.y = y
;cat.flux = flux
cat.flux = flux
cat.sharp = sharp
cat.round = round
cat.round2 = round
cat.mag = reform(mags[0,*])
cat.err = reform(errap[0,*])
cat.sky = sky
cat.skysig = skyerr

; Get only good sources
; 1. Must have good magnitudes
; 2. sky level must not be too low
; 3. std.dev. in sky annulus must not be too high
cat_orig = cat
;stdskyerr = mad(cat.skyerr)    ; std.dev. in the sky error
;gd = where(cat.mag lt 50. and abs(cat.sky-skymode) lt (5.0*skysig) and $
;             cat.skyerr lt (skysig+5.0*stdskyerr),ngd)
gd = where(cat.mag lt 50.,ngd)
if (ngd eq 0) then begin
  print,'NO good stars found'
  error = 'NO good stars found'
  undefine,cat
  return
endif
cat = cat[gd]

undefine,tempim

if keyword_set(stp) then stop

end


;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


pro wcsfit_orient,ra,dec,info,x,y,stp=stp

; Getting information from INFO structure
cenra = info.cenra
cendec = info.cendec
pixscale = info.pixscale
xhalf = info.xhalf
yhalf = info.yhalf
up = info.up
left = info.left

; Use the gnomic projection
ROTSPHCEN,ra,dec,cenra,cendec,lon,lat,/gnomic

; Put on pixel scale
lon2 = lon * 3600.0d0 / pixscale
lat2 = lat * 3600.0d0 / pixscale

; What's the orientation
; UP
CASE up of
  'N': y = lat2   
  'S': y = -lat2
  'E': y = lon2
  'W': y = -lon2
  else: stop,'BAD UP'   
ENDCASE

; LEFT
CASE left of
  'N': x = -lat2
  'S': x = lat2
  'E': x = -lon2
  'W': x = lon2
  else: stop,'BAD LEFT'
ENDCASE   

; Have the coordinates start at (X/Y)=(0,0)
; Add xhalf/yhalf  
x = x + xhalf
y = y + yhalf

if keyword_set(stp) then stop

end


;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


pro wcsfit_findorient,refcat,cat,info,up,left,stp=stp

; This program finds the correct orientation
; Try the 8 different orientations

; Compile MATCHSTARS.PRO
RESOLVE_ROUTINE,'MATCHSTARS',/compile_full_file
;MATCHSTARS,/silent

; Temporary info file
info2 = info

; Ranges
minx = floor(min(cat.x))
miny = floor(min(cat.y))
maxx = ceil(max(cat.x))
maxy = ceil(max(cat.y))
xx2 = cat.x
yy2 = cat.y

n=8
uparr = ['N','N','S','S','E','E','W','W']
leftarr = ['E','W','E','W','N','S','N','S']
xcorrarr = fltarr(n)

; Looping through the possibilities
FOR i=0,n-1 do begin

  ; Orient the reference stars
  undefine,x,y
  info2.up = uparr[i]
  info2.left = leftarr[i]
  WCSFIT_ORIENT,double(refcat.raj2000),double(refcat.dej2000),info2,x,y

  ; Getting coordinates within the range
  gd1 = where(x ge minx and x le maxx and y ge miny and y le maxy,ngd1)
  xx1 = x[gd1] - minx
  yy1 = y[gd1] - miny

  ; Get the cross-correlation peak
  ; reusing the FFTs for coordinates2 and psf
  undefine,xshift,yshift,bestcorr,angle,xcorr
  if i eq 0 then undefine,fft2,fftp
  MATCHSTARS_XCORR,xx1,yy1,xx2,yy2,xshift,yshift,angle,bestcorr,xcorr,/smooth,xyscale=4,fwhm=5,$
                   fft2=fft2,fftp=fftp

  ; Need to compare bestcorr to std(xcorr)
  medxcorr = median(xcorr)
  stdxcorr = mad(xcorr)
  bestcorr2 = (bestcorr-medxcorr)/stdxcorr

  ; Put in array
  xcorrarr[i] = bestcorr2

  ;print,uparr[i],' ',leftarr[i],' ',bestcorr
  ;stop

END

bestind = first_el(maxloc(xcorrarr))
up = uparr[bestind]
left = leftarr[bestind]

info.up = up
info.left = left

;plot,xcorrarr,ps=-1
print,'BEST ORIENTATION: UP=',up,' LEFT=',left

if keyword_set(stp) then stop

end


;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro wcsfit_daomatch,refcat,cat,ind1,ind2,trans,count=count,rms=rms,stp=stp

; This does star matching with DAOMATCH which uses triangles
; on the bright stars

nrefcat = n_elements(refcat)
ncat = n_elements(cat)

undefine,ind1,ind2,trans
count=0
rms = 999999.

; Not enough inputs
if (nrefcat eq 0 or ncat eq 0) then begin
  print,'Syntax - wcsfit_daomatch,refcat,cat,ind1,ind2,trans,count=count,rms=rms,stp=stp'
  return
endif

refcat1 = refcat

; Only getting reference stars that overlap with the image
xmax = max(cat.x)
ymax = max(cat.y)
gdref1 = where(refcat1.x ge 0 and refcat1.y ge 0 and $
               refcat1.x le xmax and refcat1.y le ymax,ngdref1)
if (ngdref1 eq 0) then begin
  print,'NO OVERLAP'
  return
endif
refcat1 = refcat1[gdref1]
nrefcat1 = ngdref1

; Make the "fake" Reference ALS file
dum = {id:0L,x:0.0,y:0.0,mag:0.0,err:0.0,sky:10.0,niter:1,chi:1.0,sharp:0.0}
refals = replicate(dum,nrefcat1)
refals.id = lindgen(nrefcat1)+1
refals.x = refcat1.x
refals.y = refcat1.y

; 2MASS stars, or UCAC4 stars
if tag_exist(refcat1,'JMAG') then begin
  ;gdmag = where(finite(refcat.jmag) eq 1,ngdmag)
  ;gderr = where(finite(refcat.e_jmag) eq 1,ngderr)
  ;gdboth = where(finite(refcat.jmag) eq 1 and finite(refcat.e_jmag) eq 1,ngdboth)
  ;mag = refcat.jmag
  ;err = refcat.e_jmag
  ;
  ;; Fit the error distribution
  ;if (ngderr lt ngdmag) then begin
  ;  bd = where(finite(refcat.e_jmag) eq 0,nbd)
  ;  coef = poly_fit(10.^(mag[gdboth]/2.5),err[gdboth],1)
  ;  err[bd] = poly(10.^(mag[bd]/2.5),coef)
  ;endif
  ;refals.mag = mag
  ;refals.err = err

  refals.mag = refcat1.jmag
  refals.err = refcat1.e_jmag
endif else begin

  ; USNO-B1 stars
  if tag_exist(refcat1,'RMAG') then begin
    refals.mag = refcat1.rmag
    refals.err = refcat1.rerr
  endif else if tag_exist(refcat1,'R1MAG') then begin
    refals.mag = refcat1.r1mag
    refals.err = 0.05  ; just a guess
  endif

endelse

; Only getting reference stars with good photometry
gdref = where(finite(refals.mag) eq 1 and finite(refals.err) eq 1 and $
              refals.mag lt 50. and refals.err lt 50.,ngdref)
if (ngdref eq 0) then begin
  print,'NO REFERENCE STARS WITH GOOD PHOTOMETRY'
  return
endif
refals_orig = refals
refals = refals[gdref]
treffile = MKTEMP('dao')    ; absolute path
treffile = FILE_BASENAME(treffile[0])
FILE_DELETE,treffile,/allow
head=' NL    NX    NY  LOWBAD HIGHBAD  THRESH     AP1  PH/ADU  RNOISE    FRAD'
push,head,'  1  2048  4096   655.1 33000.0   35.52    3.00    2.30    3.09    4.94'
WRITEALS,treffile,refals,head


; Make the "fake" Image ALS file
dum = {id:0L,x:0.0,y:0.0,mag:0.0,err:0.0,sky:10.0,niter:1,chi:1.0,sharp:0.0}
catals = replicate(dum,ncat)
catals.id = cat.id
catals.x = cat.x
catals.y = cat.y
catals.mag = cat.mag
catals.err = cat.err
tcatfile = MKTEMP('dao')    ; absolute path
tcatfile = FILE_BASENAME(tcatfile[0])
FILE_DELETE,tcatfile,/allow
WRITEALS,tcatfile,catals,head


; Make a temporary script to run FIND
undefine,lines
push,lines,'#!/bin/csh'
push,lines,'daomatch << DONE'
push,lines,'${1}'
push,lines,'${1}.mch'
push,lines,'${2}/'     ; '/'=SAME SCALE
push,lines,'y'
push,lines,'DONE'
tempfile = MKTEMP('dao')    ; absolute path
WRITELINE,tempfile,lines
FILE_CHMOD,tempfile,'755'o

; Run the program
SPAWN,tempfile+' '+treffile+' '+tcatfile,out,errout
FILE_DELETE,tempfile,treffile,tcatfile,/allow    ; delete the temporary files


; Load the MCH file
mchfile = treffile+'.mch'
mchtest = FILE_TEST(mchfile)
if (mchtest eq 1) then begin
  LOADMCH,mchfile,files,daotrans
endif else begin
  print,'NO MCH FILE'
  return
endelse
FILE_DELETE,mchfile,/allow


; Apply the transformation
trans1 = reform(daotrans[1,*])
trans = trans1[0:5]       ; final transformation array
Apar = trans1[0]
Bpar = trans1[1]
Cpar = trans1[2]
Dpar = trans1[3]
Epar = trans1[4]
Fpar = trans1[5]

; from ccdpck.txt
;              x(1) = A + C*x(n) + E*y(n)
;              y(1) = B + D*x(n) + F*y(n)

; Apply transformation
xout = Apar + Cpar*cat.x + Epar*cat.y
yout = Bpar + Dpar*cat.x + Fpar*cat.y


; Get the matches
SRCMATCH,refcat.x,refcat.y,xout,yout,10.0,ind1a,ind2a,count=count1

if (count1 lt 2) then begin
  print,'DAOMATCH - NOT ENOUGH MATCHES'
  return
endif

; Calculate First RMS
xdiff1 = refcat[ind1a].x-xout[ind2a]
ydiff1 = refcat[ind1a].y-yout[ind2a]
xmed1 = median(xdiff1,/even)
ymed1 = median(ydiff1,/even)
diff1 = sqrt( (xdiff1-xmed1)^2.0 + (ydiff1-ymed1)^2.0 )
rms1 = sqrt( mad(xdiff1)^2.0 + mad(ydiff1)^2.0)

; Second match
dcr = ceil(rms1)
SRCMATCH,refcat.x,refcat.y,xout,yout,dcr,ind1,ind2,count=count

; Calculate final RMS
xdiff = refcat[ind1].x-xout[ind2]
ydiff = refcat[ind1].y-yout[ind2]
xmed = median(xdiff,/even)
ymed = median(ydiff,/even)
diff = sqrt( (xdiff-xmed)^2.0 + (ydiff-ymed)^2.0 )
rms = sqrt( mad(xdiff)^2.0 + mad(ydiff)^2.0)

; Print output
print,'DAOMATCH trans = ',string(trans1[0:5],format='(6F10.4)')
print,strtrim(count,2),' matches within ',strtrim(dcr,2),' pixels'
print,'RMS = ',strtrim(rms,2),' pixels'

if keyword_set(stp) then stop

end


;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


pro wcsfit_initwcs,cat3,refcat3,proj,head,crpix1=crpix1,crpix2=crpix2,error=error,stp=stp

; This puts an initial WCS in the header
;
; See http://iraf.noao.edu/projects/ccdmosaic/tnx.html
; for more info on the TNX WCS coordinate system

undefine,error

; Error Handling
;------------------
; Establish error handler. When errors occur, the index of the  
; error is returned in the variable Error_status:  
CATCH, Error_status 
;This statement begins the error handler:  
if (Error_status ne 0) then begin 
   print,'WCSFIT_INITWCS ERROR: ', !ERROR_STATE.MSG  
   error = !ERROR_STATE.MSG
   CATCH, /CANCEL 
   return
endif


; Check inputs
ncat3 = n_elements(cat3)
nrefcat3 = n_elements(refcat3)
nproj = n_elements(proj)
nhead = n_elements(head)

; Not enough inputs
if (ncat3 eq 0 or nrefcat3 eq 0 or nproj eq 0 or nhead eq 0) then begin
  print,'Syntax - wcsfit_initwcs,cat3,refcat3,proj,head,crpix1=crpix1,crpix2=crpix2,stp=stp'
  error = 'Not enough inputs'
  return
endif

; Need at least 3 stars
if (ncat3 lt 3) then begin
  print,'WCSFIT_INITWCS needs at least 3 stars'
  error = 'WCSFIT_INITWCS needs at least 3 stars'
  return
endif

; Randomly pick 3
RANDOMIZE,cat3.x,3,dum,indx=indx
ra3 = double(refcat3[indx].raj2000)
dec3 = double(refcat3[indx].dej2000)
x3 = cat3[indx].x
y3 = cat3[indx].y


;-----------
; WCS types
;-----------
case strupcase(proj) of

;----------------
; TNX Projection
;----------------
'TNX': begin


  ; This header already has a TNX WCS
  ;----------------------------------
  ctype1 = sxpar(head,'CTYPE1',/silent)
  if (strtrim(ctype1,2) ne '0') then begin


    ; The CRPIX *CANNOT* be changed if the header
    ; already has a WCS because the TNX distortion
    ; parameters were defined with specific CRPIX
    ; values.  Changing CRPIX would change the distortion!!!


    ; Get the TNX astrometry structure
    wcs = HDR2WCSTNX(head)

    ; Getting new CD matrix and CRVAL
    ;  Need to get iterate to get the best solution
    ;  sometimes the CD matrix that STARAST gives is
    ;  not that great.
    niter = 20 ;10
    rmsarr = fltarr(niter)
    cdarr = fltarr(niter,2,2)
    crvalarr = fltarr(niter,2)

    for i=0,niter-1 do begin

      ; Randomly pick 3
      RANDOMIZE,cat3.x,3,dum,indx=indx
      ra3 = double(refcat3[indx].raj2000)
      dec3 = double(refcat3[indx].dej2000)
      x3 = cat3[indx].x
      y3 = cat3[indx].y


      ; Remove the distortions

      ; Get TNX coordinates
      WCSTNX_RD2XY,ra3,dec3,wcs,tnx_x3,tnx_y3,/degree

      ; Create TAN astrometry structure
      ;  Make fake TAN header and extrast astrometry structure
      ;  This makes a difference in the PV1 array, PV1_2=90
      ;  but for TNX you get PV1_2=0 which is not what we want
      head2 = REPSTR(head,'TNX','TAN')
      EXTAST,head2,astr

      ; Get TAN coordinates
      AD2XY,ra3,dec3,astr,tan_x3,tan_y3

      ; The difference b/w TNX and TAN are the distortion terms
      ; Distortion terms
      x_distortion = tnx_x3 - tan_x3
      y_distortion = tnx_y3 - tan_y3
      ;x3b = x3 + x_distortion
      ;y3b = y3 + y_distortion
      x3b = x3 - x_distortion
      y3b = y3 - y_distortion

      ; Starting new TNX WCS structure
      wcs2 = wcs    

      ; Use STARAST to get the CD matrix
      STARAST,ra3,dec3,x3b,y3b,cd,projection='TAN'
      wcs2.ast.cd = cd

      ; Get NEW CRVAL, use all matched stars
      WCSTNX_XY2RD,cat3.x,cat3.y,wcs2,tra,tdec,/degree
      radiff = double(refcat3.raj2000) - tra
      decdiff = double(refcat3.dej2000) - tdec
      wcs2.ast.crval[0] = wcs2.ast.crval[0] + median(radiff)
      wcs2.ast.crval[1] = wcs2.ast.crval[1] + median(decdiff)

      ; Calculating RMS for this solution
      WCSTNX_XY2RD,cat3.x,cat3.y,wcs2,tra,tdec,/degree
      diff = SPHDIST(double(refcat3.raj2000),double(refcat3.dej2000),tra,tdec,/deg)
      rmsarr[i] = sqrt(mean(diff^2.0)) * 3600.0

      ; Saving the CD and CRVAL values
      cdarr[i,*,*] = cd
      crvalarr[i,*] = wcs2.ast.crval

    end

    ; Getting best CD solution
    bestind = first_el(minloc(rmsarr))
    bestcd = reform(cdarr[bestind,*,*])
    bestcrval = reform(crvalarr[bestind,*,*])

    print,'INIT WCS rms = ',rmsarr[bestind],' arcsec'

    wcs2 = wcs
    wcs2.ast.cd = bestcd
    wcs2.ast.crval = bestcrval
          
    ; Put new WCS into header
    WCSTNX2HDR, head, wcs2    


    ;HEAD_XYAD,head,cat3.x,cat3.y,aa,dd,/degree
    ;plot,refcat3.raj2000,(refcat3.dej2000-dd)*3600.,ps=1
    ;diff = SPHDIST(refcat3.raj2000,refcat3.dej2000,aa,dd,/deg)
    ;print,sqrt(mean(diff^2.))*3600.0

    ;stop


  ; NEW TNX WCS
  ;-------------
  endif else begin

    ; Use TAN for now
    temphead = head
    WCSFIT_INITWCS,cat3,refcat3,'TAN',temphead,crpix1=crpix1,crpix2=crpix2

    ; Get TAN astrometry structure
    EXTAST,temphead,astr

    ; Create TNX structure with zero distortion
    ;
    ; 1=chebyshev, 2=legendre, 3=polynomial.
    ; 
    ; This is how the xi/eta min/max values are used
    ; to normalize the xi/eta coordinates:
    ;   xin = (2 * xi - (ximax + ximin)) / (ximax - ximin)
    ;   etan = (2 * eta - (etamax + etamin)) / (etamax - etamin)
    ; max=+1 and min=-1 will leave xi/eta unchanged.
    tnx = {fun_type:3, xiorder:4, etaorder:4, cross_type:2, ximin:-1.0d0, ximax:1.0d0,$
            etamin:-1.0d0, etamax:1.0d0, m:[0,1,2,3,0,1,2,0,1,0], n:[0,0,0,0,1,1,1,2,2,3],$
            c:dblarr(10)}

    ; Create WCS structure
    wcs = {ast:astr, tnx1:tnx, tnx2:tnx, ccdsec:intarr(4), datasec:intarr(4),$
           crvaloffset_deg:dblarr(2), astcd:dblarr(2,2)}

    ; Change projection from TAN to TNX
    wcs.ast.ctype[0] = REPSTR(wcs.ast.ctype[0],'TAN','TNX')
    wcs.ast.ctype[1] = REPSTR(wcs.ast.ctype[1],'TAN','TNX')

    ; Put the WCS into the header
    WCSTNX2HDR, head, wcs   

  endelse  ; new TNX WCS

End ; TNX

'TPV': begin

  ; This header already has a TPV WCS
  ;----------------------------------
  ctype1 = sxpar(head,'CTYPE1',/silent)
  if (strtrim(ctype1,2) ne '0') then begin


    ; The CRPIX *CANNOT* be changed if the header
    ; already has a WCS because the TPV distortion
    ; parameters were defined with specific CRPIX
    ; values.  Changing CRPIX would change the distortion!!!


    ; Get the TPV astrometry structure
    wcs = HDR2WCSTPV(head)

    ; Getting new CD matrix and CRVAL
    ;  Need to get iterate to get the best solution
    ;  sometimes the CD matrix that STARAST gives is
    ;  not that great.
    niter = 20 ;10
    rmsarr = fltarr(niter)
    cdarr = fltarr(niter,2,2)
    crvalarr = fltarr(niter,2)

    for i=0,niter-1 do begin

      ; Randomly pick 3
      RANDOMIZE,cat3.x,3,dum,indx=indx
      ra3 = double(refcat3[indx].raj2000)
      dec3 = double(refcat3[indx].dej2000)
      x3 = cat3[indx].x
      y3 = cat3[indx].y


      ; Remove the distortions

      ; Get TPV coordinates
      WCSTPV_RD2XY,ra3,dec3,wcs,tpv_x3,tpv_y3,/degree

      ; Create TAN astrometry structure
      astr = wcs.ast
      astr.ctype[0] = REPSTR(astr.ctype[0],'TPV','TAN')
      astr.ctype[1] = REPSTR(astr.ctype[1],'TPV','TAN')

      ; Get TAN coordinates
      AD2XY,ra3,dec3,astr,tan_x3,tan_y3

      ; The difference b/w TPV and TAN are the distortion terms
      ; Distortion terms
      x_distortion = tpv_x3 - tan_x3
      y_distortion = tpv_y3 - tan_y3
      ;x3b = x3 + x_distortion
      ;y3b = y3 + y_distortion
      x3b = x3 - x_distortion
      y3b = y3 - y_distortion

      ; Starting new TPV WCS structure
      wcs2 = wcs    

      ; Use STARAST to get the CD matrix
      STARAST,ra3,dec3,x3b,y3b,cd,projection='TAN'
      wcs2.ast.cd = cd

      ; Get NEW CRVAL, use all matched stars
      WCSTPV_XY2RD,cat3.x,cat3.y,wcs2,tra,tdec,/degree
      radiff = double(refcat3.raj2000) - tra
      decdiff = double(refcat3.dej2000) - tdec
      wcs2.ast.crval[0] = wcs2.ast.crval[0] + median(radiff)
      wcs2.ast.crval[1] = wcs2.ast.crval[1] + median(decdiff)

      ; Calculating RMS for this solution
      WCSTPV_XY2RD,cat3.x,cat3.y,wcs2,tra,tdec,/degree
      diff = SPHDIST(double(refcat3.raj2000),double(refcat3.dej2000),tra,tdec,/deg)
      rmsarr[i] = sqrt(mean(diff^2.0)) * 3600.0

      ; Saving the CD and CRVAL values
      cdarr[i,*,*] = cd
      crvalarr[i,*] = wcs2.ast.crval

    end

    ; Getting best CD solution
    bestind = first_el(minloc(rmsarr))
    bestcd = reform(cdarr[bestind,*,*])
    bestcrval = reform(crvalarr[bestind,*,*])
    bestrms = rmsarr[bestind]

    wcs2 = wcs
    wcs2.ast.cd = bestcd
    wcs2.ast.crval = bestcrval

    ; RMS too high, this seems to have often for DECam images
    if rmsarr[bestind] gt 1.0 then begin
      ; Just fix the CRVAL values
      WCSTPV_XY2RD,cat3.x,cat3.y,wcs,tra,tdec,/degree
      radiff = double(refcat3.raj2000) - tra
      decdiff = double(refcat3.dej2000) - tdec
      wcs2 = wcs
      wcs2.ast.crval[0] = wcs2.ast.crval[0] + median(radiff)
      wcs2.ast.crval[1] = wcs2.ast.crval[1] + median(decdiff)

      ; Calculating RMS for this solution
      WCSTPV_XY2RD,cat3.x,cat3.y,wcs2,tra,tdec,/degree
      diff = SPHDIST(double(refcat3.raj2000),double(refcat3.dej2000),tra,tdec,/deg)
      bestrms = sqrt(mean(diff^2.0)) * 3600.0
    endif

    print,'INIT WCS rms = ',bestrms,' arcsec'

    ; Put new WCS into header
    WCSTPV2HDR, head, wcs2    


    ;HEAD_XYAD,head,cat3.x,cat3.y,aa,dd,/degree
    ;plot,refcat3.raj2000,(refcat3.dej2000-dd)*3600.,ps=1
    ;diff = SPHDIST(refcat3.raj2000,refcat3.dej2000,aa,dd,/deg)
    ;print,sqrt(mean(diff^2.))*3600.0


  ; NEW TPV WCS
  ;-------------
  endif else begin

    ; Use TAN for now
    temphead = head
    WCSFIT_INITWCS,cat3,refcat3,'TAN',temphead,crpix1=crpix1,crpix2=crpix2

    ; Get TAN astrometry structure
    EXTAST,temphead,astr

    ; Create TPV structure with zero distortion
    ; Create WCS structure
    pv1 = fltarr(40)
    pv1[1] = 1
    wcs = {ast:astr, pv1:pv1, pv2:pv1, ccdsec:intarr(4), datasec:intarr(4),$
           crvaloffset_deg:dblarr(2), astcd:dblarr(2,2)}

    ; Change projection from TAN to TPV
    wcs.ast.ctype[0] = REPSTR(wcs.ast.ctype[0],'TAN','TPV')
    wcs.ast.ctype[1] = REPSTR(wcs.ast.ctype[1],'TAN','TPV')

    ; Put the WCS into the header
    WCSTPV2HDR, head, wcs   

  endelse  ; new TPV WCS

End

;----------------------
; "Normal" projection
;----------------------
else: begin

  ; Get the intial astrometry structure
  extast,head,astr,noparams
  if noparams lt 1 then begin
    temphead = head
    STARAST,ra3,dec3,x3,y3,h=temphead,projection=proj
    EXTAST,temphead,astr
  endif

  ; Adding input CRPIX values
  if n_elements(crpix1) gt 0 and n_elements(crpix2) gt 0 then $
  astr.crpix = [crpix1,crpix2]

  ; Getting new CD matrix and CRVAL
  ;  Need to iterate to get the best solution
  ;  sometimes the CD matrix that STARAST returns is
  ;  not that great.
  niter = 20  ;10
  rmsarr = fltarr(niter)
  cdarr = fltarr(niter,2,2)
  crvalarr = fltarr(niter,2)

  for i=0,niter-1 do begin

    ; Randomly pick 3
    RANDOMIZE,cat3.x,3,dum,indx=indx
    ra3 = double(refcat3[indx].raj2000)
    dec3 = double(refcat3[indx].dej2000)
    x3 = cat3[indx].x
    y3 = cat3[indx].y

    ; Use STARAST to get the CD matrix
    STARAST,ra3,dec3,x3,y3,cd,projection=proj

    astr2 = astr
    astr.cd = cd

    ; Get NEW CRVAL, use all matched stars
    XY2AD,cat3.x,cat3.y,astr2,tra,tdec
    radiff = double(refcat3.raj2000) - tra
    decdiff = double(refcat3.dej2000) - tdec
    astr2.crval[0] = astr2.crval[0] + median(radiff)
    astr2.crval[1] = astr2.crval[1] + median(decdiff)

    ; Calculating RMS for this solution
    XY2AD,cat3.x,cat3.y,astr2,tra,tdec
    diff = SPHDIST(double(refcat3.raj2000),double(refcat3.dej2000),tra,tdec,/deg)
    rmsarr[i] = sqrt(mean(diff^2.0)) * 3600.0

    ; Saving the CD and CRVAL values
    cdarr[i,*,*] = astr2.cd
    crvalarr[i,*] = astr2.crval

  end


  ; Getting best CD solution
  bestind = first_el(minloc(rmsarr))
  bestcd = reform(cdarr[bestind,*,*])
  bestcrval = reform(crvalarr[bestind,*,*])

  print,'INIT WCS rms = ',rmsarr[bestind],' arcsec'

  astr2 = astr
  astr2.cd = bestcd
  astr2.crval = bestcrval
          
  ; Put new WCS into header
  PUTAST, head, astr2

End

Endcase


;stop

if keyword_set(stp) then stop

end


;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


pro wcsfit_refine,head,refcat,cat,iterate=iterate,maxiter=maxiter,error=error,rms=rms,nmatch=nmatch,stp=stp

; This refines the WCS solution
;
; The X/Y coordiantes in "cat" must be in the IDL convention (first
; pixel is 0).
;
; Copied from c40match.pro

undefine,error

; Error Handling
;------------------
; Establish error handler. When errors occur, the index of the  
; error is returned in the variable Error_status:  
CATCH, Error_status 

;This statement begins the error handler:  
if (Error_status ne 0) then begin 
   print,'WCSFIT_REFINE ERROR: ', !ERROR_STATE.MSG  
   error = !ERROR_STATE.MSG
   CATCH, /CANCEL 
   return
endif

; Check inputs
nhead = n_elements(head)
nrefcat = n_elements(refcat)
ncat = n_elements(cat)

if (nhead eq 0 or nrefcat eq 0 or ncat eq 0) then begin
  print,'Syntax - wcsfit_refine,head,refcat,cat,iterate=iterate,maxiter=maxiter,error=error'
  error = 'Not enough inputs'
  return
endif


; Initial matches
;refcatM = refcat1
;catM = cat1
refcatM = refcat
catM = cat
nmatch = n_elements(refcatM)

;------------------------------------------------
; Fitting the astrometry, matching stars, and 
; rejecting outliers
;------------------------------------------------


; Extract the astrometry from the header
EXTAST,head,astr,noparams


; Get TNX WCS structure
projtype = strmid(astr.ctype[0],5,3)
if projtype eq 'TNX' then wcs = HDR2WCSTNX(head)
if projtype eq 'TPV' then wcs = HDR2WCSTPV(head)


print,''
print,'FITTING WCS'
print,''

; Outlier rejection Loop
flag = 0
count = 0
nastr = astr
rms0 = 99.
chi0 = 99.
WHILE (flag ne 1) do begin

  ; Now fit with MPFIT

  ; Other inputs
  ; Get all stars for this chip
  ftol = 1d-10   ;1d-10

  ; Find close matches
  if (count gt 0) then begin

    ; TNX projection
    case projtype of
    'TNX': begin
      nwcs = wcs
      nwcs.ast = nastr
      WCSTNX_XY2RD,cat.x,cat.y,nwcs,ra2,dec2,/degree
    end
    'TPV': begin
      nwcs = wcs
      nwcs.ast = nastr
      WCSTPV_XY2RD,cat.x,cat.y,nwcs,ra2,dec2,/degree
    end
    else: XY2AD,cat.x,cat.y,nastr,ra2,dec2
    endcase

    dcr = (3.0*rms) < 3.0

    ;if count ge 1 then dcr=1.0
    ;if count ge 2 then dcr=0.7
    ;if count ge 3 then dcr=0.5

    ; Match them
    SRCMATCH,double(refcat.raj2000),double(refcat.dej2000),ra2,dec2,dcr,ind1,ind2,count=nind1,/sph
    if nind1 eq 0 then begin
      print,'NO MATCHES'
      return    
    endif
    refcatM = refcat[ind1]
    catM = cat[ind2]
    nmatch = n_elements(catM)
    ;print,strtrim(ncat3,2),' matches found'

    radiff = double(refcatM.raj2000)-ra2[ind2]
    rarms = mad(radiff)
    decdiff = double(refcatM.dej2000)-dec2[ind2]
    decrms = mad(decdiff)
    absdiff = sphdist(double(refcatM.raj2000),double(refcatM.dej2000),ra2[ind2],dec2[ind2],/deg)*3600.

    ; Removing outliers                                                                                                                                
    gd = where(abs(radiff-median(radiff)) le 3.0*rarms AND $
               abs(decdiff-median(decdiff)) le 3.0*decrms,ngd)
    ; need at least 6 stars/constraints                                                                                                                
    if ngd lt 6 then begin
      gd = where(abs(radiff-median(radiff)) le 5.0*rarms AND $
                 abs(decdiff-median(decdiff)) le 5.0*decrms,ngd)
    endif
    if ngd lt 6 then begin
      gd = (sort(absdiff))[0:n_elements(absdiff)-1]
      ngd = n_elements(gd)
    endif
  
    refcatM = refcatM[gd]
    catM = catM[gd]
  endif

  ; Inputs
  if projtype eq 'TNX' or projtype eq 'TPV' then begin
    fa = {x:double(catM.x), y:double(catM.y), ra:double(refcatM.raj2000), dec:double(refcatM.dej2000), astr:astr, wcs:wcs}
  endif else begin
    fa = {x:double(catM.x), y:double(catM.y), ra:double(refcatM.raj2000), dec:double(refcatM.dej2000), astr:astr}
  endelse

  ; Initial parameters
  ;par = [reform(nastr.cd[0,*]), reform(nastr.cd[1,*]), nastr.crpix, nastr.crval]
  if n_elements(catM) ge 6 then begin
    par = [reform(nastr.cd[0,*]), reform(nastr.cd[1,*]), nastr.crval]
  endif else begin
    ; not enough points to fit
    ;  just fix x/y offset and scale/rotation
    theta = atan(nastr.cd[0,1],nastr.cd[0,0])
    scale = nastr.cd[0,0] / cos(theta)
    par = [scale, theta, nastr.crval]
    ; cd[0,0] = scale*cos(theta)
    ; cd[0,1] = scale*sin(theta)
    ; cd[1,0] = -scale*sin(theta)
    ; cd[1,1] = scale*cos(theta)
  endelse

  ;parinfo = replicate({RELSTEP:0.0001},8)

  ; FITTING
  func = 'wcsfit_devcoo'
  maxit = 50L
  fpar = MPFIT(func, par, functargs=fa, perror=perror,niter=iter,status=status,$
               bestnorm=chisq, parinfo=parinfo, dof=dof, ftol=ftol, maxiter=maxit,/quiet)

  if (status lt 1) then begin
    print,'MPFIT error'
    error = 'MPFIT error'
    rms = 999999.
    return
  endif

  ; Getting new fitted values
  ;nastr = astr
  if n_elements(fpar) gt 4 then begin
    nastr.cd[0,*] = fpar[0:1]
    nastr.cd[1,*] = fpar[2:3]
    nastr.crval = fpar[4:5]
  endif else begin
    fscale = fpar[0]
    ftheta = fpar[1]
    nastr.cd[0,0] = fscale*cos(ftheta)
    nastr.cd[0,1] = fscale*sin(ftheta)
    nastr.cd[1,0] = -fscale*sin(ftheta)
    nastr.cd[1,1] = fscale*cos(ftheta)
    nastr.crval = fpar[2:3]
  endelse
  ;nastr.crpix = fpar[4:5]          ; These is fixed now
  ;nastr.crval = fpar[6:7]

  ; TNX projection
  case projtype of
  'TNX': begin
    nwcs = wcs
    nwcs.ast = nastr
    WCSTNX_XY2RD,catM.x,catM.y,nwcs,nra,ndec,/degree
  end
  'TPV': begin
    nwcs = wcs
    nwcs.ast = nastr
    WCSTPV_XY2RD,catM.x,catM.y,nwcs,nra,ndec,/degree
  end
  else:  XY2AD,catM.x,catM.y,nastr,nra,ndec
  endcase
  
  ; Error information
  diff = sphdist(nra,ndec,double(refcatM.raj2000),double(refcatM.dej2000),/deg)*3600.d0
  rms = sqrt(mean(diff^2.0))
  ;rms = sqrt(median(diff^2.0))
  sigpar = perror * sqrt(chisq/dof)
  chi = sqrt(chisq/dof)

  ; Print info
  print,'I=',strtrim(count+1,2),'  Nmatch=',strtrim(nmatch,2),'  RMS=',$
        string(rms,format='(F5.3)'),' arcseconds  Niter=',strtrim(iter,2)

  drms = rms0 - rms
  dchi = chi0 - chi

  if drms/rms lt 0.01 and dchi/chi lt 0.01 then flag=1
  ;if drms/rms lt 0.01 and dchi/chi lt 0.01 and count gt 0 then flag=1
  if (count+1) ge maxiter then flag=1
  if not keyword_set(iterate) and not keyword_set(maxiter) then flag=1

  rms0 = rms
  chi0 = chi

  count++

  ;stop

ENDWHILE       ; outlier rejection loop


; Final Statistics
; Calculating residuals
raresid = (double(refcatM.raj2000)-nra)*cos(double(refcatM.dej2000)/!radeg)*3600.
decresid = (double(refcatM.dej2000)-ndec)*3600.
resid = sqrt( raresid^2.0 + decresid^2.0 )
rms =  sqrt( mean( resid^2. ) )
;rms =  sqrt( median( resid^2. ) )

; Printing statistics 
if not keyword_set(silent) then begin
  form = '(A-15,A7,A8)'   
  print,''
  print,'** FIT STATISTICS **'
  print,'RA:'
  print,format=form,'Mean Resid',stringize(mean(raresid),ndec=3),' arcsec'
  print,format=form,'Stdev Resid',stringize(stdev(raresid),ndec=3),' arcsec'
  print,format=form,'Max Abs(Resid)',stringize(max(abs(raresid)),ndec=3),' arcsec'
  print,'DEC:'
  print,format=form,'Mean Resid',stringize(mean(decresid),ndec=3),' arcsec'
  print,format=form,'Stdev Resid',stringize(stdev(decresid),ndec=3),' arcsec'
  print,format=form,'Max Abs(Resid)',stringize(max(abs(decresid)),ndec=3),' arcsec'
  print,''
  print,format=form,'Total RMS',stringize(rms,ndec=3),' arcsec'
  print,''
endif ; not /silent


print,''
print,'Final WCS'
print,'-------------------------'
print,'CTYPE1 = ',strtrim(nastr.ctype[0],2)
print,'CTYPE2 = ',strtrim(nastr.ctype[1],2)
print,'CRVAL1 = ',strtrim(nastr.crval[0],2)
print,'CRVAL2 = ',strtrim(nastr.crval[1],2)
print,'CRPIX1 = ',strtrim(nastr.crpix[0],2)
print,'CRPIX2 = ',strtrim(nastr.crpix[1],2)
print,'CD1_1 =  ',strtrim(nastr.cd[0,0],2)
print,'CD1_2 =  ',strtrim(nastr.cd[0,1],2)
print,'CD2_1 =  ',strtrim(nastr.cd[1,0],2)
print,'CD2_2 =  ',strtrim(nastr.cd[1,1],2)
print,'-------------------------'

; Write the new WCS to the FITS file
head_orig = head

; TNX projection
case projtype of
'TNX': WCSTNX2HDR,head,nwcs
'TPV': WCSTPV2HDR,head,nwcs
else: PUTAST,head,nastr
endcase

;stop

if keyword_set(stp) then stop

end


;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


pro wcsfit,input,up=up0,left=left0,pixscale=pixscale0,cenra=cenra0,cendec=cendec0,outcat=cat,outrefcat=refcat,$
           inpcat=inpcat0,inprefcat=inprefcat0,stp=stp,error=error,projection=projection,noupdate=noupdate,$
           redo=redo,searchdist=searchdist,rmslim=rmslim,refname=refname,maxshift=maxshift,$
           refmaglim=refmaglim,catmaglim=catmaglim,caterrlim=caterrlim,inpfwhm=inpfwhm,head=head

t0 = systime(1)
undefine,error

; Not enough inputs
ninput = n_elements(input)
if ninput eq 0 then begin
  print,'Syntax - wcsfit,input,up=up,left=left,pixscale=pixscale,cenra=cenra,cendec=cendec,'
  print,'                inpcat=inpcat,inprefcat=inprefcat,stp=stp,error=error,projection=projection,'
  print,'                noupdate=noupdate,redo=redo,searchdist=searchdist,rmslim=rmslim,'
  print,'                refname=refname,refmaglim=refmaglim,catmaglim=catmaglim,caterrlim=caterrlim,'
  print,'                inpfwhm=inpfwhm,outcat=outcat,outrefcat=outrefcat'
  error = 'Not enough inputs'
  return
endif

; Error Handling
;------------------
; Establish error handler. When errors occur, the index of the  
; error is returned in the variable Error_status:  
CATCH, Error_status 

;This statement begins the error handler:  
if (Error_status ne 0) then begin 
   print,'WCSFIT ERROR: ', !ERROR_STATE.MSG  
   error = !ERROR_STATE.MSG
   CATCH, /CANCEL 
   return
endif


; Loading the input
LOADINPUT,input,files,count=nfiles

; Multiple files input
if nfiles gt 1 then begin
  for i=0,nfiles-1 do $
  WCSFIT,files[i],up=up0,left=left0,pixscale=pixscale0,cenra=cenra0,cendec=cendec0,$
           inpcat=inpcat0,inprefcat=inprefcat0,stp=stp,projection=projection,noupdate=noupdate,$
           redo=redo,searchdist=searchdist,rmslim=rmslim,refname=refname,maxshift=maxshift,$
           refmaglim=refmaglim,catmaglim=catmaglim,caterrlim=caterrlim,inpfwhm=inpfwhm
  return
endif

; Single file
filename = files[0]

; Does the file exist?
test = file_test(files)
if test eq 0 then begin
  print,'FILE ',filename,' NOT FOUND'
  error = 'FILE '+filename+' NOT FOUND'
  return
endif

; Fpack-compressed FITS file
if strmid(filename,6,7,/reverse_offset) eq 'fits.fz' then begin
  fpack = 1
  filebase = FILE_BASENAME(filename,'.fits.fz')
; Normal FITS file
endif else begin
  fpack = 0
  filebase = FILE_BASENAME(filename,'.fits')
endelse

; Defaults
if n_elements(noupdate) eq 0 then noupdate=0       ; updating the header?
if n_elements(rmslim) eq 0 then rmslim=1.0         ; maximum RMS to allow in arcsec


print,''
print,'========================================='
print,'RUNNING WCSFIT ON >>',filename,'<<'
print,'========================================='


; Getting image header
FITS_READ,filename,im,head,/no_abort,message=message
if message ne '' then begin
  error = 'Problem loading '+filename
  print,error
  return    
endif
if fpack eq 1 then begin
  head = headfits(filename,exten=1)  ; fits_read will modify the header improperly
  ; NAXIS1/NAXIS2 get modified by fpack, need to temporarily
  ;  put back the originals which are saved in ZNAXIS1/2
  ;  But save the fpack header so we can put things back
  ;  when we update the header at the end
  orig_head = head
  sxaddpar,head,'NAXIS1',sxpar(head,'ZNAXIS1')
  sxaddpar,head,'NAXIS2',sxpar(head,'ZNAXIS2')
endif
im = float(im)

; Information structure
info = {up:'', left:'', cenra:0.0d0, cendec:0.0d0, nx:0L, ny:0L,$
         pixscale:0.0, xhalf:0L, yhalf:0L}

; Does this image already have a WCS?
ctype1 = SXPAR(head,'CTYPE1',/silent)
ctype2 = SXPAR(head,'CTYPE2',/silent)
EXTAST,head,astr
nastr = n_elements(astr)
; We have a WCS
if nastr gt 0 then begin
  print,'This image has a WCS already'

  ; Getting orientation
  ; Left-Right
  HEAD_XYAD,head,[0,1],[0,0],a1,d1,/degree
  diff_a1 = a1[1]-a1[0]
  diff_a1 = diff_a1*cos(d1[0]/!radeg) * 3600.0
  diff_d1 = (d1[1]-d1[0])*3600.0
  ; X is east-west
  if abs(diff_a1) gt abs(diff_d1) then begin
    if diff_a1 gt 0 then left='W' else left='E'
  ; X is north-south
  endif else begin
    if diff_d1 gt 0 then left='S' else left='N'
  endelse
  if n_elements(left0) eq 0 then left0=left

  ; Up-Down
  HEAD_XYAD,head,[0,0],[0,1],a2,d2,/degree
  diff_a2 = a2[1]-a2[0]
  diff_a2 = diff_a2*cos(d2[0]/!radeg) * 3600.0
  diff_d2 = (d2[1]-d2[0])*3600.0
  ; Y is east-west
  if abs(diff_a2) gt abs(diff_d2) then begin
    if diff_a2 gt 0 then up='E' else up='W'
  ; Y is north-south
  endif else begin
    if diff_d2 gt 0 then up='N' else up='S'
  endelse
  if n_elements(up0) eq 0 then up0=up

  ; Projection type
  c1type = SXPAR(head,'CTYPE1',/silent)
  ;c1type = astr.ctype[0]
  projhead = first_el(strsplit(c1type,'-',/extract),/last)
  projhead = strtrim(projhead,2)

  ; Getting central RA/DEC
  if fpack eq 0 then begin
    nx = SXPAR(head,'NAXIS1',/silent)
    ny = SXPAR(head,'NAXIS2',/silent)
  endif else begin
    nx = SXPAR(head,'ZNAXIS1',/silent)
    ny = SXPAR(head,'ZNAXIS2',/silent)
  endelse
  HEAD_XYAD,head,0.5*nx,0.5*ny,cenra0,cendec0,/degree

endif

; WCS Projection Type
if n_elements(projection) eq 0 and n_elements(projhead) ne 0 then projection=projhead
if n_elements(projection) eq 0 then projection = 'TAN'


; Orientation Defaults
up='' & left=''
gdup=-1 & gdleft=-1
; UP INPUT
nup0 = n_elements(up0)
if nup0 gt 0 then dum = where(strtrim(up0,2) ne '',nup0)
if nup0 gt 0 then begin
  up = strmid(strtrim(strupcase(up0),2),0,1)

  ; Are the inputs okay
  gdup = first_el(where(strpos(['N','S','E','W'],up) ne -1,ngdup))
  ; Problem with UP
  if ngdup eq 0 then begin
    print,'UP MUST be one of N,S,E,W.  Using N by default'
    up = 'N'
    gdup = 0
  endif
endif else begin  ; UP input
   print,'No UP input.  Using N by default'
   up = 'N'
   gdup = 0
endelse
; LEFT INPUT
nleft0 = n_elements(left0)
if nleft0 gt 0 then dum = where(strtrim(left0,2) ne '',nleft0)
if nleft0 gt 0 then begin
  left = strmid(strtrim(strupcase(left0),2),0,1)

  ; Are the inputs okay
  gdleft = first_el(where(strpos(['N','S','E','W'],left) ne -1,ngdleft))
  ; Problem with LEFT
  if ngdleft eq 0 then begin
    print,'LEFT MUST be one of N,S,E,W.  Using E by default'
    left = 'E'
    gdleft = 2
  endif
endif else begin
  print,'No LEFT input.  Using E by default'
  left = 'E'
  gdleft = 2
endelse
; One of EACH
if ((gdup le 1 and gdleft le 1) or (gdup ge 2 and gdleft ge 2)) and $
   (gdup ne -1 and gdleft ne -1) then begin
  print,'IF UP=N/S then LEFT must be either E/W.  Fitting the orientation'
  up=''
  left=''
endif

info.up = up
info.left = left
if info.up eq '' then info.left=''
if info.left eq '' then info.up=''


; Getting pixel scale
if n_elements(pixscale0) eq 0 then begin
  GETPIXSCALE,filename,pixscale
  if pixscale lt 0 then begin
    print,'NO PIXEL SCALE'
    error = 'NO PIXEL SCALE'
    return
  endif
endif else pixscale = pixscale0
info.pixscale = pixscale

; Getting central RA
if n_elements(cenra0) eq 0 then begin
  sra = sxpar(head,'RA',/silent)

  ; RA not found, check CRVAL1
  if strtrim(sra,2) eq '0' then begin
    ctype1 = sxpar(head,'CTYPE1',/silent)
    crval1 = sxpar(head,'CRVAL1',/silent)

    ; Use CRVAL1
    if strupcase(strmid(strtrim(ctype1,2),0,2)) eq 'RA' and $
       strtrim(crval1,2) ne '0' then begin
       ra = double(crval1)       ; Assume it's in DEGREES

    ; CRVAL1 not there or not RA, Check CRVAL2
    endif else begin
      ctype2 = sxpar(head,'CTYPE2',/silent)
      crval2 = sxpar(head,'CRVAL2',/silent)

      ; It's there
      if strupcase(strmid(strtrim(ctype2,2),0,2)) eq 'RA' and $
         strtrim(crval2,2) ne '0' then begin
         ra = double(crval2)       ; Assume it's in DEGREES
      endif

    endelse

  ; RA FOUND
  endif else begin

    ra = double(sexig2ten(sra))
    ; If there is a colon then assume it is in HOURS
    colon = first_el(strpos(sra,':'))
    if colon ge 0 then ra=ra*15.0d0

  endelse

  ; NO RA
  if n_elements(ra) eq 0 then begin
    print,'NO RA'
    error = 'NO RA'
    return
  endif
  cenra = ra

endif else cenra=cenra0


; Getting central DEC
if n_elements(cendec0) eq 0 then begin
  sdec = sxpar(head,'DEC',/silent)
 
  ; DEC not found, check CRVAL2
  if strtrim(sdec,2) eq '0' then begin
    ctype2 = sxpar(head,'CTYPE2',/silent)
    crval2 = sxpar(head,'CRVAL2',/silent)

    ; Use CRVAL2
    if strupcase(strmid(strtrim(ctype2,2),0,3)) eq 'DEC' and $
       strtrim(crval2,2) ne '0' then begin
       dec = double(crval2)

    ; CRVAL2 not there or not DEC, Check CRVAL1
    endif else begin
      ctype1 = sxpar(head,'CTYPE1',/silent)
      crval1 = sxpar(head,'CRVAL1',/silent)

      ; It's there
      if strupcase(strmid(strtrim(ctype1,2),0,3)) eq 'DEC' and $
         strtrim(crval1,2) ne '0' then begin
         dec = double(crval1)
      endif  

    endelse

  ; DEC FOUND
  endif else begin
    dec = double(sexig2ten(sdec))
  endelse

  ; NO DEC
  if n_elements(dec) eq 0 then begin
    print,'NO DEC'
    error = 'NO DEC'
    return
  endif
  cendec = dec

endif else cendec=cendec0
info.cenra = cenra
info.cendec = cendec


; Get image size
nx = sxpar(head,'NAXIS1',/silent)
ny = sxpar(head,'NAXIS2',/silent)
; Get it from the image
if strtrim(nx,2) eq '0' or strtrim(ny,2) eq '0' then begin
  ;FITS_READ,filename,im,head
  sz = size(im)
  nx = sz[1]
  ny = sz[2]
endif
info.nx = nx
info.ny = ny

; Central pixel of image
xhalf = round(nx/2)
yhalf = round(ny/2)
info.xhalf = xhalf
info.yhalf = yhalf

; Print out image information
obj = SXPAR(head,'OBJECT',count=nobj,/silent)
if nobj eq 0 then obj=''
print,''
print,'OBJECT: ',obj
print,'ORIENTATION: UP-',up,', LEFT-',left
print,'PIXSCALE: ',strtrim(pixscale,2),' "/pixel'
print,'AREA: ',strtrim(long(info.nx*pixscale/60.),2),'x',strtrim(long(info.ny*pixscale/60.),2),' arcmin'
print,'Central RA=',strtrim(cenra,2),', DEC=',strtrim(cendec,2),' degrees'
print,'WCS Projection Type = ',projection
print,''


;########################################################
; STEP 1:  Get star lists
;########################################################
print,'------------------------------------------------------'
print,'STEP 1:  Get star lists'
print,'------------------------------------------------------'

; Check if USNO catalog input or already saved to file
; Does the image have a star list

;--------------------------------------------------------
; Image Catalog
;--------------------------------------------------------
if n_elements(inpcat0) gt 0 then begin

  type = size(inpcat0,/type)
  if type eq 8 then begin

    ; Check for X/Y tags
    tags = tag_names(inpcat0)
    xtag = where(stregex(tags,'^X',/boolean) eq 1,nxtag)
    ytag = where(stregex(tags,'^Y',/boolean) eq 1,nytag)

    ; We have X/Y tags
    if nxtag ge 0 and nytag ge 0 then begin
      cat = inpcat0
      ncat = n_elements(cat)

      print,'Using INPUT image star list, Nstars=',strtrim(ncat,2)
    endif

    if n_elements(cat) eq 0 then print,'INPUT catalog does NOT have X/Y tags'

  ; inpcat0 is NOT a structure
  endif else begin
    print,'INPUT catalog is NOT a structure'
  endelse

endif

; NO catalog input, check if there is a previously saved catalog
if n_elements(cat) eq 0 and not keyword_set(redo) then begin
  
  ; Check if there is already a DAOPHOT/FIND file for this file
  catfile = filebase+'_cat.dat'
  test = FILE_TEST(catfile)
  
  ; There IS a CAT file
  if (test eq 1) then begin
    RESTORE,catfile

    ; Check structure
    ncat = n_elements(cat)
    type = size(cat,/type)
    xtag = -1 & ytag = -1
    if ncat gt 0 and type eq 8 then begin
      tags = tag_names(cat)
      xtag = where(stregex(tags,'^X',/boolean) eq 1,nxtag)
      ytag = where(stregex(tags,'^Y',/boolean) eq 1,nytag)
    endif

    ; Catalog OKAY
    if ncat gt 0 and type eq 8 and nxtag ge 0 and nytag ge 0 then begin
      print,'Using PREVIOUSLY SAVED catalog file ',catfile,', Nstars=',strtrim(ncat,2)
    endif else begin
      print,'PROBLEMS with previously saved catalog ',catfile
      undefine,cat,ncat
    endelse

  endif ; NO previously saved catalog
endif  ; check for previously saved catalog file


; NO catalog input, get star list from DAOPHOT FIND or (SExtractor)
if n_elements(cat) eq 0 then begin

  ;print,'NO IMAGE star list.  Getting star list using DAOPHOT/FIND'
  print,'NO IMAGE star list.  Getting star list using FIND/APER'

  ; MASK OUT BAD PIXELS AND COLUMNS!!!!


 ; ; Run DAOPHOT FIND
 ; ;-----------------
 ;
 ; ; Make a .opt file
 ; base = file_basename(filename,'.fits')
 ; optfile = base+'.opt'
 ; test = FILE_TEST(optfile)
 ; if test eq 0 then begin
 ;   MKOPT,filename,fwhm=fwhm
 ; endif else begin
 ;   READLINE,optfile,optlines
 ;   fwhmind = where(stregex(optlines,'FW =',/boolean) eq 1,nfwhmind)
 ;   fwhmarr = strsplit(optlines[fwhmind[0]],'=',/extract)
 ;   fwhm = float(fwhmarr[1])
 ; endelse
 ;
 ; ; Copy the .opt file daophot.opt 
 ; if FILE_TEST('daophot.opt') eq 0 then $
 ;   FILE_COPY,base+'.opt','daophot.opt',/over
 ;
 ; ; Make temporary photo.opt file
 ; tphotofile = MKTEMP('photo')
 ; undefine,photlines
 ; push,photlines,'A1 = '+STRING(fwhm*3.0,format='(F7.4)')
 ; push,photlines,'IS = 45.0000'
 ; push,photlines,'OS = 50.0000'
 ; WRITELINE,tphotofile,photlines
 ;
 ; ; Make a temporary script to run FIND
 ; undefine,lines
 ; push,lines,'#!/bin/sh'
 ; push,lines,'daophot="/net/astro/bin/daophot"'
 ; push,lines,'export image=${1}'
 ; push,lines,'rm ${image}.temp.log      >& /dev/null'
 ; push,lines,'rm ${image}.temp.coo      >& /dev/null'
 ; push,lines,'rm ${image}.temp.ap       >& /dev/null'
 ; push,lines,'daophot << END_DAOPHOT >> ${image}.temp.log'
 ; push,lines,'OPTIONS'
 ; push,lines,'${image}.opt'
 ; push,lines,''
 ; push,lines,'ATTACH ${image}.fits'
 ; push,lines,'FIND'
 ; push,lines,'1,1'
 ; push,lines,'${image}.temp.coo'
 ; push,lines,'y'
 ; push,lines,'PHOTOMETRY'
 ; push,lines,FILE_BASENAME(tphotofile)
 ; push,lines,''
 ; push,lines,'${image}.temp.coo'
 ; push,lines,'${image}.temp.ap'
 ; push,lines,'EXIT'
 ; push,lines,'END_DAOPHOT'
 ; ;tempfile = maketemp('dao','.sh')
 ; tempfile = MKTEMP('dao')    ; absolute path
 ; WRITELINE,tempfile,lines
 ; FILE_CHMOD,tempfile,'755'o
 ;
 ; ; Run the program
 ; SPAWN,tempfile+' '+base,out,errout
 ; FILE_DELETE,tempfile    ; delete the temporary script
 ;
 ; ; Test the coo and ap file
 ; cootest = FILE_TEST(base+'.temp.coo')
 ; if cootest eq 1 then coolines=FILE_LINES(base+'.temp.coo') else coolines=0
 ; aptest = FILE_TEST(base+'.temp.ap')
 ; if aptest eq 1 then aplines=FILE_LINES(base+'.temp.ap') else aplines=0
 ;
 ;
 ; ; DAOPHOT ran properly
 ; if (coolines ge 4 and aplines ge 4) then begin
 ;
 ;   ; Load the coordinates file
 ;   LOADCOO,base+'.temp.coo',coo,coohead
 ;   ncoo = n_elements(coo)
 ;
 ;   ; Load the aperture photometry file
 ;   LOADAPER,base+'.temp.ap',aper,aperhead
 ;
 ;   ; Create the CAT structure
 ;   dum = {id:0L,x:0.0d0,y:0.0d0,mag:0.0,err:0.0,sky:0.0,skysig:0.0,sharp:0.0,round:0.0,round2:0.0}
 ;   cat = replicate(dum,ncoo)
 ;   cat.id = coo.id
 ;   cat.x = coo.x
 ;   cat.y = coo.y
 ;   cat.sharp = coo.sharp
 ;   cat.round = coo.round
 ;   cat.round2 = coo.round2
 ;   cat.mag = aper.mag[0]
 ;   cat.err = aper.err[0]
 ;   cat.sky = aper.sky
 ;   cat.skysig = aper.skysig
 ;
 ;   ; Remove the temporary files
 ;   FILE_DELETE,base+['.temp.log','.temp.coo','.temp.ap'],/allow
 ;   FILE_DELETE,'daophot.opt',/allow
 ;   FILE_DELETE,base+['.opt','.als.opt'],/allow
 ;   FILE_DELETE,tphotofile,/allow
 ;   junk = FILE_TEST(base+'jnk.fits')
 ;   if junk eq 1 then FILE_DELETE,base+'jnk.fits'
 ;
 ; ; DAOPHOT Problem, run IDL FIND/APER
 ; endif else begin
 ;
 ;   print,'DAOPHOT problem.  Using IDL FIND/APER.PRO'
 ;
 ;   ; Remove the temporary DAOPHOT files
 ;   FILE_DELETE,base+['.temp.log','.temp.coo','.temp.ap'],/allow
 ;   FILE_DELETE,'daophot.opt',/allow
 ;   FILE_DELETE,base+['.opt','.als.opt'],/allow
 ;   junk = FILE_TEST(base+'jnk.fits')
 ;   if junk eq 1 then FILE_DELETE,base+'jnk.fits'
 ;
 ;   ; Find the sources
 ;   WCSFIT_FIND,filename,cat,error=finderror
 ;
 ;   ; There was an error
 ;   if n_elements(finderror) gt 0 then begin
 ;     print,finderror
 ;     error = finderror
 ;     return
 ;   endif
 ;
 ; endelse

  ; Find the sources
  WCSFIT_FIND,filename,cat,error=finderror,inpfwhm=inpfwhm

  ; There was an error
  if n_elements(finderror) gt 0 then begin
    print,finderror
    error = finderror
    return
  endif

  ; Only keep "good" stars
  ;gd = where(cat.sharp ge 0.2 and cat.sharp le 1.0 and $
  ;           cat.round ge -1.0 and cat.round le 1.0 and $
  ; appropriate SHARP and ROUND cuts are already applied
  ; in WCSFIT_FIND
  gd = where(cat.mag lt 50.0 and cat.err lt 1.0,ngd)
  if (ngd eq 0) then begin
    print,'NO good stars'
    error = 'NO good stars'
    return
  endif
  cat_orig = cat
  cat = cat[gd]
  ncat = n_elements(cat)
  print,'Nstars = ',strtrim(ncat,2)


  ; The coordinates are in FITS/IRAF format (first pixel is 1)
  ; Changing to IDL convention (first pixel is 0)
  cat.x = cat.x - 1.0
  cat.y = cat.y - 1.0

  ; SAVING THE CATALOG FILE
  catfile = filebase+'_cat.dat'
  print,'Saving catalog file to FILE=',catfile
  SAVE,cat,file=catfile

endif  ; no catalog input
print,''


;--------------------------------------------------------
; Reference Catalog
;--------------------------------------------------------

; Reference catalog input
if n_elements(inprefcat0) gt 0 then begin

  type = size(inprefcat0,/type)
  if type eq 8 then begin

    ; Check for RAJ2000/DEJ2000 or RA/DEC tags
    tags = tag_names(inprefcat0)
    ra2000tag = first_el(where(strpos(tags,'RAJ2000') ne -1,nra2000tag))
    de2000tag = first_el(where(strpos(tags,'DEJ2000') ne -1,nde2000tag))

    ; We have RAJ2000/DEJ2000
    if nra2000tag ge 0 and nde2000tag ge 0 then begin
      refcat = inprefcat0
      nrefcat = n_elements(refcat)

      print,'Using INPUT reference star list, Nstars=',strtrim(nrefcat,2)

    ; NO RAJ2000/DEJ2000, check RA/DEC
    endif else begin

      ratag = first_el(where(strpos(tags,'RA') ne -1,nratag))
      dectag = first_el(where(strpos(tags,'DEC') ne -1,ndectag))

      ; The structure is okay, adding the RAJ2000/DEJ2000 tags
      if nratag ge 0 and ndectag ge 0 then begin
       
        refcat = inprefcat0
        nrefcat = n_elements(refcat)

        ADD_TAG,refcat,'RAJ2000',0.0d0,refcat
        ADD_TAG,refcat,'DEJ2000',0.0d0,refcat
        refcat.raj2000 = refcat.ra
        refcat.dej2000 = refcat.dec

        print,'Uing INPUT reference star list, Nstars=',strtrim(nrefcat,2)
        print,'Using RA/DEC tags in INPUT structure' 

      endif  ; structure okay
    endelse  ; no raj2000/dej2000

    if n_elements(refcat) eq 0 then print,'INPUT reference structure does NOT have RAJ2000/DEJ2000 or RA/DEC tags'

  ; inprefcat0 is NOT a structure
  endif else begin
    print,'INPUT reference catalog is NOT a structure'
  endelse

endif

; NO reference catalog input, check if there is a previously saved catalog
if n_elements(refcat) eq 0 and not keyword_set(redo) then begin
      
  ; Check if there is already a reference file for this file
  refcatfile = filebase+'_refcat.dat'
  test = file_test(refcatfile)

  ; There IS a reference file
  if (test eq 1) then begin
    RESTORE,refcatfile

    ; Check structure
    nrefcat = n_elements(refcat)
    type = size(refcat,/type)
    ratag = -1 & dectag = -1
    if nrefcat gt 0 and type eq 8 then begin
      reftags = tag_names(refcat)
      ra2000tag = first_el(where(strpos(reftags,'RAJ2000') ne -1,nra2000tag))
      de2000tag = first_el(where(strpos(reftags,'DEJ2000') ne -1,nde2000tag))
    endif
  
    ; Catalog OKAY
    if nrefcat gt 0 and type eq 8 and nra2000tag ge 0 and nde2000tag ge 0 then begin
      print,'Using PREVIOUSLY SAVED reference catalog file ',refcatfile,', Nstars=',strtrim(nrefcat,2)
    endif else begin 
      print,'PROBLEMS with previously saved reference catalog ',refcatfile
      undefine,refcat,nrefcat
    endelse
  
  endif ; NO previously saved reference catalog
endif  ; check for previously saved reference catalog file

; NO reference catalog input
if n_elements(refcat) eq 0 then begin

  ; Getting X/Y sizes in arcmin
  xdist = ceil( (pixscale*nx)/60. )
  ydist = ceil( (pixscale*ny)/60. )
  dist = xdist > ydist   ; use symmetric distance since we don't know the orientation yet
  ;dist = dist*2.0        ; get 100% more just to be safe
  ;dist = dist*4.0        ; get 100% more just to be safe
  if keyword_set(searchdist) then dist=float(searchdist)      ; input value
  ;dist = dist > 60.0      ; get at least 1 degree
  ;dist = dist > 30.0      ; get at least 1/2 degree
  ;dist = dist > 45.0      ; get at least 3/4 degree
  dist = dist > 10.0       ; get at least 10 arcmin

  ; Querying the catalog
  refcatname = 'USNO-B1'    ; the default  
  if keyword_set(refname) then refcatname=refname
  if refcatname ne 'USNO-B1' and refcatname ne '2MASS-PSC' and refcatname ne 'UCAC4' and refcatname ne 'GAIA/GAIA' then refcatname='USNO-B1'

  print,'NO Reference Catalog Input: QUERYING ',refcatname,' Catalog',$
       '  Area:',strtrim(long(dist),2),'x',strtrim(long(dist),2),' arcmin'
  cfa = 1  ;0 ; 1
  userefcatname = refcatname
  if refcatname eq '2MASS-PSC' and cfa eq 1 then userefcatname='II/246'   ; cfa issue
  refcat = QUERYVIZIER(userefcatname, [cenra,cendec], [dist,dist], cfa=cfa, /allcolumns)
  nrefcat = n_elements(refcat)
  type = size(refcat,/type)

  ; Problem
  if type[0] ne 8 then begin
    print,'QUERYVIZIER ERROR: Probably a problem with RA/DEC, or the WEB connection'
    error = 'QUERYVIZIER ERROR: Probably a problem with RA/DEC, or the WEB connection'
    return
  endif


  ; Add Rmag, Bmag, Rerr, Berr
  if (refcatname eq 'USNO-B1') then begin

    refcat_orig = refcat    

    ; Getting average Rmag, Bmag 
    rmag = median([[refcat.r1mag],[refcat.r2mag]],dim=2,/even)
    bmag = median([[refcat.b1mag],[refcat.b2mag]],dim=2,/even)
    rmagbd = where(finite(rmag) eq 0,nrmagbd)      ; replace NaNs with 99.9999
    if nrmagbd gt 0 then rmag[rmagbd] = 99.9999
    bmagbd = where(finite(bmag) eq 0,nbmagbd)
    if nbmagbd gt 0 then bmag[bmagbd] = 99.9999

    ; Getting Err, Berr
    rerr = abs(refcat.r1mag-refcat.r2mag)
    berr = abs(refcat.b1mag-refcat.b2mag)
    rerrbd = where(finite(rerr) eq 0,nrerrbd)      ; replace NaNs with 0.2
    if nrerrbd gt 0 then rerr[rerrbd] = 0.2
    if nrmagbd gt 0 then rerr[rmagbd] = 9.9999     ; bad mag -> bad error
    berrbd = where(finite(berr) eq 0,nberrbd)
    if nberrbd gt 0 then berr[berrbd] = 0.2
    if nbmagbd gt 0 then berr[bmagbd] = 9.9999     ; bad mag -> bad error

    ; Adding the BMAG, RMAG, BERR, RERR tags
    ADD_TAG,refcat,'BMAG',0.0,refcat
    ADD_TAG,refcat,'RMAG',0.0,refcat
    ADD_TAG,refcat,'BERR',0.0,refcat
    ADD_TAG,refcat,'RERR',0.0,refcat
    refcat.bmag = bmag
    refcat.rmag = rmag
    refcat.berr = berr
    refcat.rerr = rerr

    ; Remove stars with bad photometry
    gd = where(refcat.bmag lt 50. and refcat.rmag lt 50. and $
               refcat.berr lt 2. and refcat.rerr lt 2.,ngd)
    refcat = refcat[gd]

  endif


  ; 2MASS catalog
  if (refcatname eq '2MASS-PSC') then begin

    ; Maybe an error cut, but otherwise they are ALL good.
    ; gd = where(refcat.e_kmag lt 0.2,ngd)
    ; refcat = refcat[gd]

  end

  ; GAIA
  if (refcatname eq 'GAIA/GAIA') then begin
    ; they are probably all good
    ; Add RAJ2000 and DEJ2000 tags
    ADD_TAG,refcat,'RAJ2000',0.0d0,refcat
    ADD_TAG,refcat,'DEJ2000',0.0d0,refcat
    refcat.raj2000 = refcat.ra_icrs
    refcat.dej2000 = refcat.de_icrs
  endif

  ; UCAC4 ??

  nrefcat = n_elements(refcat)
  print,'Nstars = ',strtrim(nrefcat,2)

  ; SAVING THE REFERENCE CATALOG FILE
  refcatfile = filebase+'_refcat.dat'
  print,'Saving reference catalog file to FILE=',refcatfile
  SAVE,refcat,file=refcatfile

endif

; Reference catalog type
if n_elements(refname) eq 0 then begin
  ; 2MASS-PSC has _2MASS tag
  if tag_exist(refcat,'_2MASS') then refname='2MASS-PSC'
  ; USNO-B1 has USNO_B1_0 tag
  if tag_exist(refcat,'USNO_B1_0') then refname='USNO-B1'
  ; UCAC4 has UCAC4 tag
  if tag_exist(refcat,'UCAC4') then refname='UCAC4'
  if n_elements(refname) gt 0 then print,'Reference catalog type is ',refname else print,'Reference catalog type UNKNOWN'
endif


;########################################################
; STEP 2:  Convert references stars RA/DEC to pixel X/Y
; Need cenra, cendec, pixscale, orientation
;########################################################
print,'------------------------------------------------------'
print,'STEP 2:  Convert reference stars RA/DEC to pixel X/Y'
print,'------------------------------------------------------'

; IF THE IMAGE ALREADY HAS A WCS THEN USE THAT AS AN INITIAL GUESS
;stop

print,'UP=',info.up
print,'LEFT=',info.left

; Getting the properly oriented reference X/Y coordiantes
WCSFIT_ORIENT,double(refcat.raj2000),double(refcat.dej2000),info,x,y

; Add X/Y to the reference structure
reftags = tag_names(refcat)
xtag = first_el(where(strpos(reftags,'X') ne -1,nxtag))
if nxtag ge 0 then ADD_TAG,refcat,'X',0.0,refcat
ytag = first_el(where(strpos(reftags,'Y') ne -1,nytag))
if nytag ge 0 then ADD_TAG,refcat,'Y',0.0,refcat
refcat.x = x
refcat.y = y

; Use the WCS in the header to get initial X/Y coordinates
usedheadxy = 0
ctype1 = SXPAR(head,'CTYPE1',/silent)
if (strmid(ctype1,5,3) eq 'TNX' or strupcase(projection) eq 'TNX') or $
   (strmid(ctype1,5,3) eq 'TPV' or strupcase(projection) eq 'TPV') then begin
  print,'Using header WCS to get initial X/Y coordinates for Reference stars'

  ;HEAD_ADXY,head,double(refcat.raj2000),double(refcat.dej2000),x,y,/degree
  ADXY,head,double(refcat.raj2000),double(refcat.dej2000),x,y

  refcat.x = x
  refcat.y = y

  usedheadxy = 1
endif


;########################################################
; STEP 3:  Match star lists
; use matchstars.pro
;########################################################
print,'------------------------------------------------------'
print,'STEP 3:  Matching Star Lists'
print,'------------------------------------------------------'

; Keeping only REFERENCE stars with accurate coordinates
dum = where(stregex(reftags,'E_RAJ2000',/boolean) eq 1,ne_ra)
dum = where(stregex(reftags,'E_DEJ2000',/boolean) eq 1,ne_de)
; The structure has coordinate errors
if (ne_ra gt 0 and ne_de gt 0) then begin
  medera = median(refcat.e_raj2000)
  stdera = mad(refcat.e_raj2000)
  mededec = median(refcat.e_dej2000)
  stdedec = mad(refcat.e_dej2000)
  ralim = (medera+3.0*stdera) < 500
  declim = (mededec+3.0*stdedec) < 500
  gd = where(refcat.e_raj2000 le ralim and $
              refcat.e_dej2000 le declim,ngd)
  refcat1b = refcat[gd]

; No coordinate errors
endif else begin
  refcat1b = refcat
endelse

; Keeping only REFERENCE stars with low proper motions
dum = where(stregex(reftags,'PMRA',/boolean) eq 1,npmra)
dum = where(stregex(reftags,'PMDE',/boolean) eq 1,npmde)
; The structure has proper motion
;  GAIA-DR1 has NAN for most stars
if (npmra gt 0 and npmde gt 0) then begin
  gd = where( (finite(refcat1b.pmra) and abs(refcat1b.pmra) lt 200 and abs(refcat1b.pmde) lt 200) or $
              (not finite(refcat1b.pmra)),ngd)
  refcat1b = refcat1b[gd]
endif

; Keeping only REFERENCE stars within the magnitude limit
if tag_exist(refcat1b,'JMAG') then begin
  if n_elements(refmaglim) gt 0 then maglim=refmaglim else maglim=16.5
  gd = where(refcat1b.jmag lt maglim,ngd)
  print,'Keeping only REFERENCE stars with JMAG < ',strtrim(maglim,2),'  ',strtrim(ngd,2),' sources'
  refcat1b = refcat1b[gd]
endif else begin
  gd = lindgen(n_elements(refcat1b))
  if n_elements(refmaglim) gt 0 then maglim=refmaglim else maglim = 21.0
  if tag_exist(refcat1b,'RMAG') then begin
    gd = where(refcat1b.rmag lt maglim,ngd)
    print,'Keeping only REFERENCE stars with RMAG < ',strtrim(maglim,2),'  ',strtrim(ngd,2),' sources'
  endif else begin
    if tag_exist(refcat1b,'R1MAG') then begin
      gd = where(refcat1b.r1mag lt maglim,ngd)
      print,'Keeping only REFERENCE stars with RMAG1 < ',strtrim(maglim,2),'  ',strtrim(ngd,2),' sources'
    endif
  endelse
  refcat1b = refcat1b[gd]
endelse

; Keeping only DETECTED sources within the magnitude limit
if n_elements(catmaglim) gt 0 then begin
  gdcat = where(cat.mag lt catmaglim,ngdcat)
  print,'Keeping only DETECTED sources with MAG < ',strtrim(catmaglim,2),' ',strtrim(ngdcat,2),' sources'
  if ngdcat lt 5 then begin
    error = 'Not enough detected sources brighter than '+strtrim(catmaglim,2)+' mag to perform WCS fitting.'
    if not keyword_set(silent) then print,error
    return
  endif
  cat1 = cat[gdcat]
endif else cat1=cat

; Keeping only DETECTED sources within the error limit
if n_elements(caterrlim) gt 0 then begin
  gdcat2 = where(cat1.err lt caterrlim,ngdcat2)
  print,'Keeping only DETECTED sources with ERR < ',strtrim(caterrlim,2),' ',strtrim(ngdcat2,2),' sources'
  if ngdcat2 lt 5 then begin
    error = 'Not enough detected sources brighter than '+strtrim(catmerrlim,2)+' mag to perform WCS fitting.'
    if not keyword_set(silent) then print,error
    return
  endif
  cat1 = cat1[gdcat2]
endif

;-----------------------------
; FLOWCHART FOR MATCHING STARS
;
; -If a WCS already exists, then check if it is already good
;   . first check the RMS of stars matching (with SRCMATCH) with in 2.0 arcsec
;   . second use MATCHSTARS_XCORR
; -If no good match and the density of reference stars > than density
;  of image stars then try using just bright reference stars for the
;  matching
; -If no good match then use ALL reference stars for matching
; -If no good match and header WCS used to convert reference RA/DEC to
;  X/Y then try matching with the WCSFIT_ORIENT oriented coordinates
; -If no good match then try DAOMATCH for the matching.

RESOLVE_ROUTINE,'MATCHSTARS',/compile_full_file
nmatch = 0
matchrms = 999999.


; If a WCS already exists then check the rms
;-------------------------------------------
if (nastr gt 0) then begin

  print,'Checking WCS in header'

  ; Get X/Y coordinates for the reference stars
  HEAD_ADXY,head,refcat1b.raj2000,refcat1b.dej2000,xref,yref,/degree

  dcr = ceil(2.0/pixscale)
  SRCMATCH,xref,yref,cat1.x,cat1.y,dcr,ind1,ind2,count=nmatch
  ;if nmatch lt 10 then SRCMATCH,xref,yref,cat.x,cat.y,2*dcr,ind1,ind2,count=nmatch
  ;if nmatch lt 10 then SRCMATCH,xref,yref,cat.x,cat.y,4*dcr,ind1,ind2,count=nmatch

  ; Enough matches
  if (nmatch gt 5) then begin
    print,strtrim(nmatch,2),' MATCHES'

    xdiff = xref[ind1]-cat1[ind2].x
    ydiff = yref[ind1]-cat1[ind2].y
    xmed = median(xdiff,/even)
    ymed = median(ydiff,/even)

    ; Rematch
    if max(abs([xmed,ymed])) gt 1.0 then begin
      SRCMATCH,xref-xmed,yref-ymed,cat1.x,cat1.y,dcr,ind1,ind2,count=nmatch
      xdiff = xref[ind1]-cat1[ind2].x
      ydiff = yref[ind1]-cat1[ind2].y
      xmed = median(xdiff,/even)
      ymed = median(ydiff,/even)
    endif

    diff = sqrt( (xdiff-xmed)^2.0 + (ydiff-ymed)^2.0 )
    ;initrms = sqrt(mean(diff^2.0)) * pixscale
    initrms = sqrt( mad(xdiff)^2.0 + mad(ydiff)^2.0) * pixscale

    print,'Robust RMS = ',strtrim(initrms,2),' arcsec'

  endif else initrms=99999.

  ;; No good matches, try MATCHSTARS.PRO
  ;if (initrms ge 1.0) then begin
  ;
  ;  print,'Initial RMS bad.  Trying cross-correlation.'
  ;
  ;  MATCHSTARS,xref,yref,cat.x,cat.y,ind1,ind2,trans,maxshift=maxshift,count=nmatch,/norot
  ;
  ;  ; Enough matches
  ;  if (nmatch gt 10) then begin
  ;    xdiff = xref[ind1]-cat[ind2].x
  ;    ydiff = yref[ind1]-cat[ind2].y
  ;    xmed = median(xdiff,/even)
  ;    ymed = median(ydiff,/even)
  ;    diff = sqrt( (xdiff-xmed)^2.0 + (ydiff-ymed)^2.0 )
  ;    initrms = sqrt( mad(xdiff)^2.0 + mad(ydiff)^2.0)
  ;    print,'Robust RMS = ',strtrim(initrms,2),' arcsec'
  ;  endif else initrms=999999.
  ;endif ; no match


  ; No good matches, try MATCHSTARS_XCORR.PRO
  ;if (initrms ge 1.0) then begin
  ;if (nmatch lt 10) or (initrms ge 2.0) then begin
  ;if (nmatch lt 5) or (initrms ge 2.0) then begin
  if (nmatch lt 5) or (initrms ge 1.0) then begin

    print,'Initial RMS bad or not enough matches.  Trying cross-correlation.'

    MATCHSTARS_XCORR,xref,yref,cat1.x,cat1.y,xsh,ysh,ang,bestcorr,xcorr,xyscale=4,fwhm=4,$
                     smooth=4,nsig=nsig,matchnum=matchnum,maxshift=maxshift

    ; Good XCORR match
    if (nsig gt 7) then begin

      xref2 = xref-xsh
      yref2 = yref-ysh

      dcr = ceil(2.0/pixscale)
      SRCMATCH,xref2,yref2,cat1.x,cat1.y,dcr,ind1,ind2,count=nmatch

      if (nmatch gt 10) then begin
        xdiff = xref[ind1]-cat1[ind2].x
        ydiff = yref[ind1]-cat1[ind2].y
        xmed = median(xdiff,/even)
        ymed = median(ydiff,/even)
        diff = sqrt( (xdiff-xmed)^2.0 + (ydiff-ymed)^2.0 )
        initrms = sqrt( mad(xdiff)^2.0 + mad(ydiff)^2.0) * pixscale
        print,'Robust RMS = ',strtrim(initrms,2),' arcsec'
      endif else initrms=999999.
    endif ; good XCORR match
  endif ; no match


  ; Use the intial WCS
  ;if (initrms lt 2.0) then begin
  if (initrms lt 1.0) then begin

    print,'Using the intial WCS'
    refcat1b.x = xref
    refcat1b.y = yref

    ; Remove outliers, need at least 6 to continue
    ndiff = n_elements(diff)
    sidiff = sort(diff)
    difforder = lonarr(ndiff)
    difforder[sidiff] = lindgen(ndiff)  ; where are they in the ordered list
    bdind = where( abs(diff) gt 3.0*initrms and difforder gt 5,nbdind)
    if nbdind gt 0 then REMOVE,bdind,ind1,ind2

    refcat2 = refcat1b[ind1]
    cat2 = cat1[ind2]

    trans = [xmed, ymed, 1.0, 0.0, 0.0, 1.0]

    goto,FITWCS
    ;goto,REFINEWCS
 
  ; Initial WCS not good enough
  endif else begin
    print,'Initial WCS not good enough'
  endelse

endif else initrms=99999.

; Checking the density of stars
imarea = (info.nx*pixscale/60.0)*(info.ny*pixscale/60.0)
drefra = max(double(refcat1b.raj2000))-min(double(refcat1b.raj2000))
drefra = drefra*cos(median(double(refcat1b.dej2000))/!radeg)
drefdec = max(double(refcat1b.dej2000))-min(double(refcat1b.dej2000))
refarea = drefra*drefdec*3600.
imdensity = n_elements(cat)/imarea
refdensity = n_elements(refcat1b)/refarea
print,''
print,'Density of stars in image  = ',string(imdensity,format='(F7.2)'),' stars/arcmin^2'
print,'Density of reference stars = ',string(refdensity,format='(F7.2)'),' stars/arcmin^2'
print,''

; Using just the BRIGHT reference stars
;--------------------------------------
;if (nmatch lt 3 or matchrms*pixscale gt 1.5*rmslim) and $
if (nmatch lt 3 or initrms gt 1.5*rmslim) and (refdensity gt imdensity) then begin

  ; Sorting by magnitude
  if TAG_EXIST(refcat1b,'JMAG') then begin
    si = sort(refcat1b.jmag)
    refcat1b = refcat1b[si]
  endif else begin
    if tag_exist(refcat1b,'RMAG') then begin
      si = sort(refcat1b.rmag)
    endif else begin
      if tag_exist(refcat1b,'R1MAG') then si=sort(refcat1b.r1mag)
      if tag_exist(refcat1b,'_GMAG_') then si=sort(refcat1b._gmag_)
    endelse
    if n_elements(si) eq 0 then si=lindgen(n_elemetns(refcat1b))
    refcat1b = refcat1b[si]
  endelse

  ; Getting just the brightest reference stars
  nref = floor(imdensity*3*refarea) < n_elements(refcat1b)
  refcat1b = refcat1b
  refcat1b_bright = refcat1b[0:nref-1]

  print,'-- Matching Stars. Trying just the bright reference stars --'
  print,'Nref = ',strtrim(nref,2)
  print,''

  ; Get the matches
  MATCHSTARS,refcat1b_bright.x,refcat1b_bright.y,cat1.x,cat1.y,ind1,ind2,trans,count=nmatch,$
             rms=matchrms,maxshift=maxshift

  ; We successed
  if (nmatch ge 3 and matchrms*pixscale lt 1.5*rmslim) then begin

    ind1a = ind1
    ind2a = ind2
    nmatch1 = nmatch

    ; Transform the stars and match them to ALL the reference stars
    out = TRANS_COO(cat1.x,cat1.y,trans)
    xout = reform(out[0,*])
    yout = reform(out[1,*])
    dcr = (3.0*matchrms > 1.0) < 5.0
    SRCMATCH,refcat1b.x,refcat1b.y,xout,yout,dcr,ind1,ind2,count=nmatch

    print,strtrim(nmatch,2),' matches to ALL reference stars'
  endif

endif  ; use bright stars

; Match the stars with MATCHSTARS
;--------------------------------
if (nmatch lt 3 or matchrms*pixscale gt 1.5*rmslim) then begin

  print,''
  print,'-- Matching Stars.  Using ALL reference stars --'
  print,''

  MATCHSTARS,refcat1b.x,refcat1b.y,cat1.x,cat1.y,ind1,ind2,trans,count=nmatch,$
             rms=matchrms,maxshift=maxshift

endif


; NO MATCHES, If header WCS used, try WCSFIT_ORIENT X/Y
;------------------------------------------------------
if ((nmatch lt 3 or matchrms*pixscale gt 1.5*rmslim) and usedheadxy eq 1) then begin

  print,''
  print,'Header WCS initial X/Y coordinates gave BAD results.'
  print,'Trying basic orientation'
  print,''

  ; Getting the properly oriented reference X/Y coordinates
  WCSFIT_ORIENT,refcat1b.raj2000,refcat1b.dej2000,info,x,y
  refcat1b.x = x
  refcat1b.y = y

  ; Match the stars
  MATCHSTARS,refcat1b.x,refcat1b.y,cat1.x,cat1.y,ind1,ind2,trans,count=nmatch,$
             rms=matchrms,maxshift=maxshift

endif


; NO MATCHES, Try using DAOMATCH
;------------------------------------------------------
if (nmatch lt 3 or matchrms*pixscale gt 1.5*rmslim) then begin

  print,''
  print,'-- Trying DAOMATCH --'
  print,''

  ; Match the stars
  WCSFIT_DAOMATCH,refcat1b,cat1,ind1,ind2,trans,count=nmatch,rms=matchrms

endif


; NO matches
if nmatch lt 3 then begin
  error = 'TOO FEW MATCHES'
  print,'TOO FEW MATCHES'
  return
endif


; These are the matched catalogs
refcat2 = refcat1b[ind1]
cat2 = cat1[ind2]



;########################################################
; STEP 4:  Fit WCS
;########################################################
FITWCS:
print,'------------------------------------------------------'
print,'STEP 4:  Fitting WCS'
print,'------------------------------------------------------'


; Set CRPIX
;----------

; Input CRPIX1
if n_elements(crpix1) gt 0 then begin
  print,'Using input CRPIX1=',strtrim(crpix1,2)
  crpix1 = float(crpix1)
; CRPIX1 NOT input
endif else begin
  scrpix1 = SXPAR(head,'CRPIX1',/silent)
  ; Using value already in header
  if strtrim(scrpix1,2) ne '0' then begin
    crpix1 = float(scrpix1)
    print,'Using header CRPIX1=',strtrim(crpix1,2)
  ; Using center of image
  endif else begin
    crpix1 = nx*0.5
    print,'Using center of image, CRPIX1=',strtrim(crpix1,2)
  endelse
endelse

; Input CRPIX2
if n_elements(crpix2) gt 0 then begin
  print,'Using input CRPIX2=',strtrim(crpix2,2)
  crpix2 = float(crpix2)
; CRPIX2 NOT input
endif else begin
  scrpix2 = SXPAR(head,'CRPIX2',/silent)
  ; Using value already in header
  if strtrim(scrpix2,2) ne '0' then begin
    crpix2 = float(scrpix2)
    print,'Using header CRPIX2=',strtrim(crpix2,2)
  ; Using center of image
  endif else begin
    crpix2 = ny*0.5
    print,'Using center of image, CRPIX2=',strtrim(crpix2,2)
  endelse
endelse


;-------------------
; Get an initial WCS
;-------------------

; Get "good" stars
diff = trans_coo_dev(trans,x1=refcat2.x,y1=refcat2.y,x2=cat2.x,y2=cat2.y)
rms = sqrt(mean(diff^2.0))
gg = where(diff lt rms,ngg)
if ngg lt 3 then gg=(sort(diff))[0:2]  ; need at least 3
refcat3 = refcat2[gg]
cat3 = cat2[gg]

head_orig = head
if not keyword_set(projection) then proj='TAN' else proj=projection
WCSFIT_INITWCS,cat3,refcat3,proj,head,crpix1=crpix1,crpix2=crpix2,error=initerror

if (n_elements(initerror) gt 0) then begin
  print,'There was an error - RETURNING'
  print,initerror
  error = initerror
  return
endif


;-------------------------------
; FITTING THE WCS WITH ALL STARS
;-------------------------------
REFINEWCS:
WCSFIT_REFINE,head,refcat2,cat2,maxiter=5,error=referror,rms=rms
if n_elements(referror) gt 0 then begin
  print,'There was an error - RETURNING'
  print,referror
  error = referror
  return
endif

; Put the RMS in the header
SXADDHIST,'WCSFIT: RMS='+string(rms,format='(F5.3)')+' arcsec on '+systime(),head
SXADDHIST,'WCSFIT: NMATCH='+strtrim(nmatch,2),head
if n_elements(refname) gt 0 then $
  SXADDHIST,'WCSFIT: Reference catalog='+strtrim(refname,2),head


;-------------------------------
; UPDATING THE FITS HEADER
;-------------------------------
if not keyword_set(noupdate) then begin

  ; The final RMS is low enough
  if (rms le rmslim) then begin

    print,'Updating WCS in ',filename
    if fpack eq 0 then begin
      MWRFITS,im,filename,head,/create
      ;FITS_WRITE,filename,im,head    ; this sometimes puts in the 2nd extension
    endif else begin
      ; Put the original NAXIS1/2 values back
      sxaddpar,head,'NAXIS1',sxpar(orig_head,'NAXIS1')
      sxaddpar,head,'NAXIS2',sxpar(orig_head,'NAXIS2')
      ; Create temporary symbolic link to make modfits.pro think
      ; this is an ordinary FITS file
      tempfile = MAKETEMP('temp')
      FILE_LINK,filename,tempfile+'.fits'
      MODFITS,tempfile+'.fits',0,head,exten_no=1,errmsg=errmsg
      FILE_DELETE,[tempfile,tempfile+'.fits'],/allow  ; delete temporary files
    endelse

  ; RMS too high
  endif else begin
    print,'RMS is greater than limit of ',strtrim(rmslim,2)
    print,'File NOT updated'
    error = 'RMS is greater than limit of '+strtrim(rmslim,2)
  endelse
endif else begin
  print,'WCS *NOT* updated'
endelse

;MODFITS,filename,0,head,errmsg=errmsg  ; this can give problems with CHECKSUM


dt = systime(1)-t0
print,'Runtime = ',strtrim(dt,2)

if keyword_set(stp) then stop

end
