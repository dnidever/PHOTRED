pro wcsfit_imacs_dummy
FORWARD_FUNCTION trans_coo, trans_coo_dev
end

;------

pro wcsfit_imacs,filename,projection=projection,logfile=logfile,error=error,stp=stp,$
                 rmslim=rmslim,noupdate=noupdate,refname=refname,searchdist=searchdist,$
                 redo=redo


; The RA/DEC in the header are for the center
; of the ENTIRE frame!!


; Get coordinates stars in all eight chips
; (or just the central four) and then run MATCHSTARS.PRO


; For ROTANGLE=67.86  the chips are approximately oriented this way:
; 6  5  8  7  Up=South, Left=East  (these are rotated by 180 degrees)
; 1  2  3  4  Up=North, Left=West
; It's off by ~23 degrees, so I think ROTANGLE=90 gives the above setup.
; I had to rotate it ~23 degrees counter-clockwise (when the image was
; flipped in X to be sky-right).

; Basically once you rotated the top row by 180 deg then they will
; all have Up=North and Left=West.
; So then everything needs to flipped in the X direction and the
; coordinates will be sky-right.
; If ROTANGLE is not ~68, then the coordinates must be rotated
; further (not sure which way).

; QUESTION;  Should the image coordinates be transformed, or should
; the USNO-B1 coordinates be transformed to the image X/Y coordinates?
; If we do the latter then it might be easier to do each chip
; separately.

; Transform them both to an intermediate X/Y coordinate system.
; Basically properly rotated, sky-right, X/Y coordinates.


t0 = systime(1)
undefine,error


; Not enough inputs
nfilename = n_elements(filename)
if nfilename eq 0 then begin
  print,'Syntax - wcsfit_imacs,file,projection=projection,logfile=logfile,error=error,stp=stp'
  return
endif


; Error Handling
;------------------
; Establish error handler. When errors occur, the index of the  
; error is returned in the variable Error_status:  
CATCH, Error_status 

;This statement begins the error handler:  
if (Error_status ne 0) then begin 
   print,'WCSFIT_IMACS ERROR: ', !ERROR_STATE.MSG  
   error = !ERROR_STATE.MSG
   CATCH, /CANCEL 
   return
endif

if n_elements(logfile) gt 0 then logf=logfile else logf=-1

printlog,logf,'FITTING WCS for IMACS image ',filename


; Defaults
if n_elements(noupdate) eq 0 then noupdate=0       ; updating the header?
if n_elements(rmslim) eq 0 then rmslim=1.0         ; maximum RMS to allow

;; Looping through the chips
;FOR i=1,8 do begin
;
;  ;i=8
;
;  print,''
;  print,'-------------------'
;  print,' FITTING CHIP ',strtrim(i,2)
;  print,'-------------------'
;  print,''
;
;  ifile = filebase+'c'+strtrim(i,2)+'.fits'

printlog,logf,'==============================================='
printlog,logf,'Running WCSFIT_IMACS on >>',filename,'<< '
printlog,logf,'==============================================='

ifile = filename
filebase = FILE_BASENAME(ifile,'.fits')

  ; Make sure the file exists
  test = FILE_TEST(ifile)
  if (test eq 0) then begin
    print,ifile,' NOT FOUND'
    error = ifile+' NOT FOUND'
    goto,BOMB
    ;return
  endif

  ; Get the image header
  head = PHOTRED_READFILE(ifile,/header)

  ; Make sure this is an IMACS image
  chip = strtrim(sxpar(head,'CHIP',/silent),2)
  ra_d = strtrim(sxpar(head,'RA-D',/silent),2)
  dec_d = strtrim(sxpar(head,'DEC-D',/silent),2)
  choffx = strtrim(SXPAR(head,'CHOFFX',/silent),2)
  choffy = strtrim(SXPAR(head,'CHOFFY',/silent),2)
  rotangle = strtrim(SXPAR(head,'ROTANGLE',/silent),2)
  if (chip eq '0' or ra_d eq '0' or dec_d eq '0' or choffx eq '0' or $
      choffy eq '0' or rotangle eq '0') then begin
    printlog,logf,ifile,' IS NOT AN IMACS IMAGE'
    error = ifile+' IS NOT AN IMACS IMAGE'
    return
  endif


  ; Getting the chip number
  chip = long(chip)


  ; Make sure that |BITPIX| > 16
  ; Only works for single chip images
  bitpix = long(SXPAR(head,'BITPIX',/silent))
  if (bitpix eq 8 or bitpix eq 16) then begin
    printlog,logfile,'BIXPIX = ',strtrim(bitpix,2),'.  Making image FLOAT'

    ; Read in the image
    im = PHOTRED_READFILE(ifile,head,error=error)

    ; Write the FLOAT image
    if n_elements(error) eq 0 and size(im,/type) lt 4 then $
      FITS_WRITE_RESOURCE,ifile,float(im),head

    ; There was a problem reading the image
    if n_elements(error) gt 0 then begin
      printlog,logf,'PROBLEM READING IN ',file
      error = 'PROBLEM READING IN '+file
      return
    endif

  endif


  ;#########################################
  ;  IMAGE STARS
  ;#########################################

  ; FIND sources in the image
  ;--------------------------
  catfile = file_basename(ifile,'.fits')+'_cat.dat'
  test = FILE_TEST(catfile)

  if (test eq 0) then begin

    printlog,logf,'Getting star coordinates for ',ifile

    ; Load the file
    im = PHOTRED_READFILE(ifile,head)
    im = float(im)

    ; Get the IMACS mask
    MK_IMACSMASK,chip,mask

    ; Mask out the bad areas
    orig_im = im
    im = im * mask

    ; Use FIND, APER
    IMFWHM,ifile,fwhm,im=im
    ;im = float(im)

    ; Find the background
    SKY,im,skymode,skysig,/silent
    if skysig lt 0.0 then skysig = mad(im)

    ; Find the sources
    ; X/Y start at 0, while daophot coordinates start at 1.
    FIND,im,x,y,flux,sharp,round,4.0*skysig,fwhm,[-1.0,1.0],[0.2,1.0],/silent

    ; Get aperture photometry
    APER,im,x,y,mags,errap,sky,skyerr,1.0,3.0*fwhm,[40,50],[10,max(im)],/silent,/meanback

    ; Leave the coordinates in IDL convention!!!!
    ;; Convert coordinates from IDL to DAOPHOT/IRAF/FITS format (0 indexed to 1 indexed)
    ;x = x + 1.0
    ;y = y + 1.0

    ; Put it all in a structure
    nstars = n_elements(x)
    dum = {id:0,x:0.0,y:0.0,flux:0.0,sharp:0.0,round:0.0,mag:0.0,err:0.0,sky:0.0,skyerr:0.0}
    cat = REPLICATE(dum,nstars)
    cat.id = lindgen(nstars)+1
    cat.x = x
    cat.y = y
    cat.flux = flux
    cat.sharp = sharp
    cat.round = round
    cat.mag = reform(mags[0,*])
    cat.err = reform(errap[0,*])
    cat.sky = sky
    cat.skyerr = skyerr

    ; Get only good sources
    ; 1. Must have good magnitudes
    ; 2. sky level must not be too low
    ; 3. std.dev. in sky annulus must not be too high
    cat_orig = cat
    stdskyerr = mad(cat.skyerr)    ; std.dev. in the sky error
    gd = where(cat.mag lt 50. and abs(cat.sky-skymode) lt (5.0*skysig) and $
               cat.skyerr lt (skysig+5.0*stdskyerr),ngd)
    cat = cat[gd]

    ; SAVING THE CATALOG FILE
    catfile = filebase+'_cat.dat'
    print,'Saving catalog file to FILE=',catfile
    SAVE,cat,file=catfile

  ; Restore previously saved catalog
  endif else begin
    RESTORE,catfile
  endelse

  
  ; Get the FITS header
  head = PHOTRED_READFILE(ifile,/header)


  ; Transform to intermediate X/Y coordinate system
  ;------------------------------------------------
  ; Sky-right.
  ADD_TAG,cat,'X2',0.0,cat
  ADD_TAG,cat,'Y2',0.0,cat

  ; Get Chip offset from the header
  ; CHOFFX/Y are the chip offsets from the center of the chip
  ; to the center of the field (in arcsec).
  ; SKY-RIGHT sense.
  ; X(fieldcen) = X(chipcen) + XOFF
  ; Y(fieldcen) = Y(chipcen) + YOFF
  xoffcen = SXPAR(head,'CHOFFX',/silent)/0.20   ; convert to pix
  yoffcen = SXPAR(head,'CHOFFY',/silent)/0.20   ; convert to pix

  ; Star Coordinates
  x = cat.x
  y = cat.y

  ; Relative to chip center, chips are 2048x4096
  x = x - 1024.0
  y = y - 2048.0

  ; Lower row, flip coordinates in X
  ; to get sky-right
  if (chip lt 5) then begin
    x = -x

  ; Upper row, flip coordinates in Y
  ; to get sky-right
  endif else begin
    y = -y
  endelse

  ; Add the offset to frame center
  x = x + xoffcen
  y = y + yoffcen


  ; Add offset to bottom left of frame
  ; The gap b/w chips is about 41 pixels.
  x = x + 4096.0 + 41.0 + 20.50
  y = y + 4096.0 + 20.50

 
  ; Add to catalog
  cat.x2 = x
  cat.y2 = y

  ;stop


  ;#########################################
  ; REFERENCE STARS
  ;#########################################

  ;--------------------------------
  ; Now get the USNO-B1 catalog
  ;--------------------------------
  cenra = SXPAR(head,'RA-D',/silent)
  cendec = SXPAR(head,'DEC-D',/silent)

  usnocat = filebase+'_refcat.dat'
  test = FILE_TEST(usnocat)

  if (test eq 0) or keyword_set(redo) then begin

    if keyword_set(searchdist) then dist=float(searchdist) else dist=40.0      ; input value
    dist = dist > 40.0      ; get at least 40 arcmin

    ; Querying the catalog
    refcatname = 'USNO-B1'    ; the default  
    if keyword_set(refname) then refcatname=refname
    if refcatname ne 'USNO-B1' and refcatname ne '2MASS-PSC' then refcatname='USNO-B1'
    print,'NO Reference Catalog Input: QUERYING ',refcatname,' Catalog',$
         '  Area:',strtrim(long(dist),2),'x',strtrim(long(dist),2),' arcmin'
    cfa = 1
    userefcatname = refcatname
    if refcatname eq '2MASS-PSC' and cfa eq 1 then userefcatname='II/246'   ; cfa issue
    usno = QUERYVIZIER(userefcatname, [cenra,cendec], [dist,dist], cfa=cfa,/allcolumns)
    nusno = n_elements(usno)
    type = size(usno,/type)

    ; Problem
    if type ne 8 then begin
      print,'QUERYVIZIER ERROR: Probably a problem with RA/DEC, or the WEB connection'
      error = 'QUERYVIZIER ERROR: Probably a problem with RA/DEC, or the WEB connection'
      return
    endif

    print,'Nstars = ',strtrim(nusno,2)

    SAVE,usno,file=usnocat
  endif else begin
    RESTORE,usnocat
    if TAG_EXIST(usno,'USNO_B1_0') then refcatname='USNO-B1' else refcatname='2MASS-PSC'
  endelse


  ; Add Rmag, Bmag, Rerr, Berr
  if (refcatname eq 'USNO-B1') then begin

    ; Only get GOOD stars
    ; Getting average Rmag, Bmag 
    rmag = median([[usno.r1mag],[usno.r2mag]],dim=2,/even)
    bmag = median([[usno.b1mag],[usno.b2mag]],dim=2,/even)
    rmagbd = where(finite(rmag) eq 0,nrmagbd)      ; replace NaNs with 99.9999
    if nrmagbd gt 0 then rmag[rmagbd] = 99.9999
    bmagbd = where(finite(bmag) eq 0,nbmagbd)
    if nbmagbd gt 0 then bmag[bmagbd] = 99.9999

    ; Getting Err, Berr
    rerr = abs(usno.r1mag-usno.r2mag)
    berr = abs(usno.b1mag-usno.b2mag)
    rerrbd = where(finite(rerr) eq 0,nrerrbd)      ; replace NaNs with 0.2
    if nrerrbd gt 0 then rerr[rerrbd] = 0.2
    if nrmagbd gt 0 then rerr[rmagbd] = 9.9999     ; bad mag -> bad error
    berrbd = where(finite(berr) eq 0,nberrbd)
    if nberrbd gt 0 then berr[berrbd] = 0.2
    if nbmagbd gt 0 then berr[bmagbd] = 9.9999     ; bad mag -> bad error

    ; Adding the BMAG, RMAG, BERR, RERR tags
    reftags = TAG_NAMES(usno)
    dum = where(reftags eq 'BMAG',nbmag)
    if (nbmag eq 0) then ADD_TAG,usno,'BMAG',0.0,usno
    dum = where(reftags eq 'RMAG',nrmag)
    if (nrmag eq 0) then ADD_TAG,usno,'RMAG',0.0,usno
    dum = where(reftags eq 'BERR',nberr)
    if (nberr eq 0) then ADD_TAG,usno,'BERR',0.0,usno
    dum = where(reftags eq 'RERR',nrerr)
    if (nrerr eq 0) then ADD_TAG,usno,'RERR',0.0,usno
    usno.bmag = bmag
    usno.rmag = rmag
    usno.berr = berr
    usno.rerr = rerr

    ; Get good USNO stars
    blim = 30.0
    errlim = 400.
    pmlim = 100.
    ugd = where(usno.e_raj2000 lt errlim and usno.e_dej2000 lt errlim and abs(usno.pmra) lt pmlim and $
                abs(usno.pmde) lt pmlim and usno.bmag lt blim,nugd)
    usno_orig = usno
    usno = usno[ugd]

  endif


  ; Add X2/Y2 tags
  tags = TAG_NAMES(usno)
  dum = where(tags eq 'X2',nxtag)
  if nxtag eq 0 then ADD_TAG,usno,'X2',0.0,usno
  dum = where(tags eq 'Y2',nytag)
  if nytag eq 0 then ADD_TAG,usno,'Y2',0.0,usno


  ; Now transform to the intermediate X/Y coordinate system

  ; Convert to X/Y
  undefine,x,y

  ; Convert RA/DEC to GNOMIC
  ROTSPHCEN,usno.raj2000,usno.dej2000,cenra,cendec,lon,lat,/gnomic

  x = -lon  ; Make it sky-right
  y = lat

  ; The plate scale is 0.200 "/pix
  x = x * 3600.0d0 / 0.200d0
  y = y * 3600.0d0 / 0.200d0

  ; Check the rotation
  rotangle = SXPAR(head,'ROTANGLE',/silent)
  diffangle = 90.0 - rotangle
  ;diffangle = rotangle - 90.0

  x2 =  x*cos(diffangle/!radeg) + y*sin(diffangle/!radeg)
  y2 = -x*sin(diffangle/!radeg) + y*cos(diffangle/!radeg)


  ; Want origin at bottom-left of frame
  ; The IMACS images are ~8192x8192, the chip gaps are ~41 pixels
  x2 = x2 + 4096.0 + 41.0 + 20.5
  y2 = y2 + 4096.0 + 20.5

  ; Add to catalog
  usno.x2 = x2
  usno.y2 = y2


  ;#########################################
  ; Now MATCH the stars
  ;#########################################

  ; Only get stars around the chip
  xmin = min(cat.x2)-1000
  xmax = max(cat.x2)+1000
  ymin = min(cat.y2)-1000
  ymax = max(cat.y2)+1000
  gdusno = where(usno.x2 ge xmin and usno.x2 le xmax and $
                 usno.y2 ge ymin and usno.y2 le ymax,ngdusno)
  usno2 = usno[gdusno]

  ;plot,usno2.x2,usno2.y2,ps=1,/ysty
  ;oplot,cat.x2,cat.y2,ps=1,co=250
  ;stop

  ; Have the coordinates start at (0,0)
  ; makes the cross-correlation easier
  xx1 = usno2.x2 - xmin
  yy1 = usno2.y2 - ymin
  xx2 = cat.x2 - xmin
  yy2 = cat.y2 - ymin

  RESOLVE_ROUTINE,'MATCHSTARS',/compile_full_file

  ;; Now match them
  ;MATCHSTARS_XCORR,xx1,yy1,xx2,yy2,xshift,yshift,angle,bestcorr,xcorr,$
  ;                 smooth=5,xyscale=4,fwhm=5,matchnum=matchnum,nsig=nsig
  ;stop

  MATCHSTARS,xx1,yy1,xx2,yy2,ind1,ind2,trans
  if n_elements(ind1) eq 0 then ind1=-1
  dum = where(ind1 ne -1,nmatch)

  ; NO matches
  if nmatch eq 0 then begin
    error = 'NO MATCHES'
    print,'NO MATCHES'
    return
  endif

  newout = trans_coo(xx2,yy2,trans)
  newx = reform(newout[0,*])
  newy = reform(newout[1,*])
  diff = sqrt( (xx1[ind1]-newx[ind2])^2.0 + (yy1[ind1]-newy[ind2])^2.0 )

  ; Matched catalogs
  refcat2b = usno2[ind1]
  cat2b = cat[ind2]


  ;stop


  ;########################################################
  ; STEP 4:  Fit WCS
  ;########################################################
  print,'------------------------------------------------------'
  print,'STEP 4:  Fitting WCS'
  print,'------------------------------------------------------'

  ;---------------------------------
  ; Get an initial WCS, with 3 stars
  ;---------------------------------

  ; Get 3 good stars, randomly picked from "good" stars
  ;diff = trans_coo_dev(trans,x1=refcat2b.x,y1=refcat2b.y,x2=cat2b.x,y2=cat2b.y)

  ; Get "good" stars
  rms = sqrt(mean(diff^2.0))
  gg = where(diff lt 2.0*rms,ngg)
  refcat3 = refcat2b[gg]
  cat3 = cat2b[gg]

  ; Randomly pick 3
  RANDOMIZE,cat3.x,3,dum,indx=indx
  ra3 = refcat3[indx].raj2000
  dec3 = refcat3[indx].dej2000
  x3 = cat3[indx].x
  y3 = cat3[indx].y

  ; You can set the projection type with =PROJECTION
  head_orig = head
  if not keyword_set(projection) then proj='TAN' else proj=projection
  ; IF THERE IS ALREADY A WCS IN THE HEADER THEN USE THAT PROJECTION TYPE
  ; unless one is input
  STARAST,ra3,dec3,x3,y3,hdr=head,projection=proj


  ; Set CRPIX and update CRVAL appropriately
  ;--------------------------------------------
  nx = 2048.
  ny = 4096.

  ; Input CRPIX1
  if n_elements(crpix1) gt 0 then begin
    print,'Using input CRPIX1=',strtrim(crpix1,2)
    crpix1 = float(crpix1)  
  ; CRPIX1 NOT input
  endif else begin
    ;scrpix1 = SXPAR(head,'CRPIX1')
    scrpix1 = SXPAR(head_orig,'CRPIX1',/silent)
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
    ;scrpix2 = SXPAR(head,'CRPIX2')
    scrpix2 = SXPAR(head_orig,'CRPIX2',/silent)
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

  ; Updating CRPIX in astrometry structure
  EXTAST,head,astr
  astr.crpix = [crpix1, crpix2]
  XY2AD,cat2b.x,cat2b.y,astr,tra,tdec
  raoff = median(refcat2b.raj2000 - tra)
  decoff = median(refcat2b.dej2000 - tdec)
  astr.crval = astr.crval + [raoff,decoff]
  PUTAST,head,astr



  ;-------------------------------
  ; FITTING THE WCS WITH ALL STARS
  ;-------------------------------
  RESOLVE_ROUTINE,'WCSFIT',/compile_full_file
  WCSFIT_REFINE,head,refcat2b,cat2b,maxiter=5,error=error,rms=rms
  if n_elements(error) gt 0 then begin
    print,'There was an error - RETURNING'
    return
  endif



  ;-------------------------------
  ; UPDATING THE FITS HEADER
  ;-------------------------------
  if not keyword_set(noupdate) then begin

    ; The final RMS is low enough
    if (rms le rmslim) then begin
      print,'Updating WCS in ',filename
      FITS_WRITE,ifile,orig_im,head

    ; RMS too high
    endif else begin
      print,'RMS is greater than limit of ',strtrim(rmslim,2)
      print,'File NOT updated'
    endelse
  endif else begin
    print,'WCS *NOT* updated'
  endelse



  ;;-------------------------------
  ;; UPDATING THE FITS HEADER
  ;;-------------------------------
  ;MODFITS,ifile,0,head,errmsg=errmsg
  ;print,'UPDATING WCS in ',ifile


  ;; Plot the residuals
  ;RESOLVE_ROUTINE,'GET_ASTROM',/compile_full_file
  ;xyad,head,cat2b.x,cat2b.y,a,d
  ;plot_astrom_resid,refcat2b.raj2000,refcat2b.dej2000,a,d


  ; Get final good matches and save them
  XYAD,head,cat.x,cat.y,a,d
  SRCMATCH,usno.raj2000,usno.dej2000,a,d,0.5,ind1,ind2,/sph
  usnom = usno[ind1]
  catm = cat[ind2]
  savefile = filebase+'c'+strtrim(i,2)+'_astrom.dat'
  SAVE,usnom,catm,file=savefile


  BOMB:

  ;stop

;END


;--------------------------------
; NOW DO A GLOBAL FIT WITH TNX
;--------------------------------
; GET_ASTROM_GLOBAL_IMACS,filebase



dt = systime(1)-t0
print,'Runtime = ',strtrim(dt,2)

if keyword_set(stp) then stop

end
