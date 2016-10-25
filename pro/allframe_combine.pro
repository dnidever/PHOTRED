;+
;
; ALLFRAME_COMBINE
;
; This combines/stacks images for ALLFRAME.
;
; You need to have run daophot, allstar, daomatch and daomaster
; already.  There needs to be a fits, opt, als.opt, ap and als
; file for each file in the MCH file.  There also needs to be an
; associated RAW file for the MCH file.
;
; INPUTS:
;  file           The MCH filename
;  =nocmbimscale  Don't scale the images when combining them.  Not
;                   recommended, but the old way of doing it.  Bright
;                   stars can be missed this way.
;  =scriptsdir    The directory that contains all of the necessary scripts.
;  =irafdir       The IRAF home directory.
;  =logfile       A logfile to print to output to.
;  /fake          Run for artificial star tests.
;  /usecmn        Use the individual cmn.lst files to construct a
;                   cmn.lst file for the combined image.
;  /stp           Stop at the end of the program
;
; OUTPUTS:
;  The combined image and mask:
;    FILEBASE_comb.fits       combined image
;    FILEBASE_comb.bpm.fits   bad pixel mask
;    FILEBASE_comb.mask.fits  SExtractor weight map (very similar to bpm)
;    FILEBASE_comb.mch        The transformations of the individual frames
;                               to the combined frame.
;  =maskdatalevel  The "bad" data level above the highest "good" value
;                    in the combined image. 
;  =filestr        Information on all of the files used for the
;                    stacked iamge.
;  =error  The error message, if there was one, else undefined
;
; USAGE:
;  IDL>allframe_combine,'ccd1001.mch',scriptsdir='/net/home/dln5q/daophot/'
;
;
; By D.Nidever   February 2008 
; Automation of steps and scripts by J.Ostheimer and Rachael Beaton
; Major upgrade of the resampling and combination code, Oct 2016
;-

pro allframe_combine,file,tile=tileinp,stp=stp,scriptsdir=scriptsdir,error=error,$
             logfile=logfile,irafdir=irafdir,satlevel=satlevel,nocmbimscale=nocmbimscale,$
             fake=fake,maskdatalevel=maskdatalevel,usecmn=usecmn,filestr=filestr

COMMON photred,setup

FORWARD_FUNCTION trans_coo, trans_coo_dev
RESOLVE_ROUTINE,'MATCHSTARS',/compile_full_file

undefine,error

; Not enough inputs
nfile = n_elements(file)
if (nfile eq 0) then begin
  print,'Syntax - allframe_combine,file,stp=stp,scriptsdir=scriptsdir,satlevel=satlevel,'
  print,'                  nocmbimscale=nocmbimscale,error=error,logfile=logfile,'
  print,'                  irafdir=irafdir,fake=fake,usecmn=usecmn'
  return
endif

; Logfile
if keyword_set(logfile) then logf=logfile else logf=-1

; Error Handling
;------------------
; Establish error handler. When errors occur, the index of the  
; error is returned in the variable Error_status:  
CATCH, Error_status 

;This statement begins the error handler:  
if (Error_status ne 0) then begin 
   print,'ALLFRAME_COMBINE ERROR: ', !ERROR_STATE.MSG  
   error = !ERROR_STATE.MSG
   CATCH, /CANCEL 
   return
endif

; How many FIND iterations
if n_elements(finditer0) eq 0 then finditer=2 else finditer=finditer0
finditer = finditer < 10  ; maximum 10.

; Saturation level
if n_elements(satlevel) eq 0 then satlevel=6e4

; Scaling of the images to be combined
if n_elements(nocmbimscale) eq 0 then nocmbimscale=0

; Getting scripts directory and iraf directory
nsetup = n_elements(setup)
if nsetup gt 0 then begin
  scriptsdir = READPAR(setup,'SCRIPTSDIR')
  irafdir = READPAR(setup,'IRAFDIR')
endif


; No irafdir
if n_elements(scriptsdir) eq 0 then begin
  printlog,logf,'SCRIPTSDIR NOT INPUT'
  error = 'SCRIPTSDIR NOT INPUT'
  return
endif

; Run CHECK_IRAF.PRO to make sure that you can run IRAF from IDL
;---------------------------------------------------------------
CHECK_IRAF,iraftest,irafdir=irafdir
if iraftest eq 0 then begin
  print,'IRAF TEST FAILED.  EXITING'
  return
endif


; No scriptsdir
if n_elements(scriptsdir) eq 0 then begin
  printlog,logf,'SCRIPTSDIR NOT INPUT'
  error = 'SCRIPTSDIR NOT INPUT'
  return
endif
; Check if the scripts exist in the current directory
scripts = ['getpsf.sh','photo.opt','apcor.opt','lstfilter','goodpsf.pro','allframe.opt',$
           'default.sex','default.param','default.nnw','default.conv']
nscripts = n_elements(scripts)
; Loop through the scripts
for i=0,nscripts-1 do begin
  info = FILE_INFO(scriptsdir+'/'+scripts[i])
  curinfo = FILE_INFO(scripts[i])

  ; No file
  if info.exists eq 0 or info.size eq 0 then begin
    printlog,logf,scriptsdir+'/'+scripts[i],' NOT FOUND or EMPTY'
    error = scriptsdir+'/'+scripts[i]+' NOT FOUND or EMPTY'
    return
  endif

  ; Check if the two files are the same size, if not copy it
  if info.size ne curinfo.size then begin
    FILE_COPY,info.name,curinfo.name,/overwrite
  endif
endfor ; scripts loop

; FILENAME
mchfile = file_basename(file)
mchdir = file_dirname(file)
mchbase = file_basename(file,'.mch')

; CD to the directory
cd,current=curdir
cd,mchdir


; Check that the mch, als, and opt files exist
mchtest = file_test(mchfile)
if mchtest eq 0 then begin
  printlog,logf,mchfile,' NOT FOUND'
  return
endif

; Checking RAW file
rawtest = file_test(mchbase+'.raw')
if rawtest eq 0 then begin
  printlog,logf,mchbase+'.raw NOT FOUND'
  return
endif


;###########################################
; CHECK NECESSARY FILES

; Load the MCH file
LOADMCH,mchfile,files,trans,magoff
nfiles = n_elements(files)

; FAKE, check that we have all the files that we need
if keyword_set(fake) then begin
  ; weights, scale, zero, comb_psf, _shift.mch
  chkfiles = mchbase+['.weights','.scale','.zero','_comb.psf']
  bdfiles = where(file_test(chkfiles) eq 0,nbdfiles)
  if nbdfiles gt 0 then begin
    error = 'FAKE.  Some necessary files not found. '+strjoin(chkfiles[bdfiles],' ')
    printlog,logf,error
    return
  endif
endif


; Final files
base = file_basename(files,'.als')
fitsfiles = base+'.fits'
outfiles = base+'.shft.fits'
outmaskfiles = base+'.mask.shft.fits'


; Gather information on all of the files
printlog,logf,'Gathering file information'
ntrans = n_elements(trans[0,*])
filestr = replicate({fitsfile:'',catfile:'',nx:0L,ny:0L,trans:dblarr(ntrans),magoff:fltarr(2),head:ptr_new(),$
                      vertices_ra:dblarr(4),vertices_dec:dblarr(4),pixscale:0.0,saturate:0.0,$
                      background:0.0,comb_zero:0.0,comb_scale:0.0,comb_weights:0.0,$
                      resampfile:'',resampmask:'',resamptrans:dblarr(ntrans),resamptransrms:0.0},nfiles)
filestr.fitsfile = fitsfiles
filestr.catfile = files
filestr.trans = transpose(trans)
filestr.magoff = magoff
filestr.resampfile = outfiles
filestr.resampmask = outmaskfiles
for i=0,nfiles-1 do begin
  FITS_READ,filestr[i].fitsfile,im1,head1,/no_abort
  filestr[i].head = ptr_new(head1)
  filestr[i].nx = sxpar(head1,'NAXIS1')
  filestr[i].ny = sxpar(head1,'NAXIS2')
  HEAD_XYAD,head1,[0,filestr[i].nx-1,filestr[i].nx-1,0],[0,0,filestr[i].ny-1,filestr[i].ny-1],vra,vdec,/degree
  filestr[i].vertices_ra = vra
  filestr[i].vertices_dec = vdec
  GETPIXSCALE,'',pixscale,head=head1
  filestr[i].pixscale = pixscale
  saturate = sxpar(head1,'SATURATE',count=nsaturate,/silent)
  if nsaturate eq 0 then saturate=50000L
  filestr[i].saturate = saturate
  gdpix = where(im1 lt saturate,ngdpix,ncomp=nbdpix)
  background = median(im1[gdpix])
  filestr[i].background = background
endfor


;#################################################
; Create default reference frame if TILE not input
if n_elements(tileinp) gt 0 then tile=tileinp else tile={type:'WCS'}
if tile.type eq 'WCS' and n_tags(tile) eq 1 then begin
  printlog,logf,'Creating TILE projection'
  ; The default projection is a tangent plane centered
  ; at halfway between the ra/dec min/max of all of
  ; the images.  The mean pixel scale is used.
  rar = minmax(filestr.vertices_ra)
  decr = minmax(filestr.vertices_dec)
  cenra = mean(rar)
  cendec = mean(decr)
  pixscale = mean(filestr.pixscale)
  ; Set up the tangent plane projection
  step = pixscale/3600.0d0
  delta_dec = range(decr)
  delta_ra = range(rar)*cos(cendec/!radeg)
  nx = ceil(delta_ra*1.01/step)
  ny = ceil(delta_dec*1.01/step)
  xref = nx/2
  yref = ny/2
  MKHDR,tilehead,fltarr(5,5)
  SXADDPAR,tilehead,'NAXIS1',nx
  SXADDPAR,tilehead,'CDELT1',step
  SXADDPAR,tilehead,'CRPIX1',xref+1L
  SXADDPAR,tilehead,'CRVAL1',cenra
  SXADDPAR,tilehead,'CTYPE1','RA---TAN'
  SXADDPAR,tilehead,'NAXIS2',ny
  SXADDPAR,tilehead,'CDELT2',step
  SXADDPAR,tilehead,'CRPIX2',yref+1L
  SXADDPAR,tilehead,'CRVAL2',cendec
  SXADDPAR,tilehead,'CTYPE2','DEC--TAN'
  EXTAST,tilehead,tileast
  tileast.equinox = 2000

  printlog,logf,'RA range = [',strtrim(rar[0],2),',',strtrim(rar[1],2),'] deg'
  printlog,logf,'DEC range = [',strtrim(decr[0],2),',',strtrim(decr[1],2),'] deg'
  printlog,logf,'Central RA = ',strtrim(cenra,2)
  printlog,logf,'Central DEC = ',strtrim(cendec,2)
  printlog,logf,'NX = ',strtrim(nx,2)
  printlog,logf,'NY = ',strtrim(ny,2)

  ; Create the TILE structure
  tile = {type:'WCS',naxis:long([nx,ny]),cdelt:double([step,step]),crpix:double([xref+1L,yref+1L]),$
          crval:double([cenra,cendec]),ctype:['RA--TAN','DEC--TAN'],$
          head:tilehead,ast:tileast,xrange:[0,nx-1],yrange:[0,ny-1],nx:nx,ny:ny}
endif

; Check that the TILE is valid
if allframe_validtile(tile,error=tilerror) eq 0 then begin
  error = tilerror
  if not keyword_set(silent) then printlog,logf,error
  return
endif

; Add HEAD/AST to TILE if needed
if tag_exist(tile,'HEAD') eq 0 then begin
  MKHDR,tilehead,fltarr(5,5)
  SXADDPAR,tilehead,'NAXIS1',tile.naxis[0]
  SXADDPAR,tilehead,'CDELT1',tile.cdelt[0]
  SXADDPAR,tilehead,'CRPIX1',tile.crpix[0]
  SXADDPAR,tilehead,'CRVAL1',tile.crval[0]
  SXADDPAR,tilehead,'CTYPE1',tile.ctype[0]
  SXADDPAR,tilehead,'NAXIS2',tile.naxis[1]
  SXADDPAR,tilehead,'CDELT2',tile.cdelt[1]
  SXADDPAR,tilehead,'CRPIX2',tile.crpix[1]
  SXADDPAR,tilehead,'CRVAL2',tile.crval[1]
  SXADDPAR,tilehead,'CTYPE2',tile.crval[1]
  if tag_exist(tile,'CD') then begin
    SXADDPAR,tilehead,'CD1_1',tile.cd[0,0]
    SXADDPAR,tilehead,'CD1_2',tile.cd[0,1]
    SXADDPAR,tilehead,'CD2_1',tile.cd[1,0]
    SXADDPAR,tilehead,'CD2_2',tile.cd[1,1]
  endif
  EXTAST,tile.head,tileast
  tileast.equinox = 2000
  tile = CREATE_STRUCT(tile,'HEAD',tilehead,'AST',tileast)
endif
if tag_exist(tile,'AST') eq 0 then begin
  EXTAST,tile.head,tileast
  tileast.equinox = 2000
  tile = CREATE_STRUCT(tile,'AST',ast)
endif
; Add XRANGE/YRANGE
if tag_exist(tile,'XRANGE') eq 0 then begin
  tile = CREATE_STRUCT(tile,'XRANGE',long([0,tile.naxis[0]-1]),'NX',tile.naxis[0])
endif
if tag_exist(tile,'YRANGE') eq 0 then begin
  tile = CREATE_STRUCT(tile,'YRANGE',long([0,tile.naxis[1]-1]),'NY',tile.naxis[1])
endif


;###########################################
; STEP 1: IMALIGN PREP

printlog,logf,'------------------------'
printlog,logf,'STEP 1: GETTING WEIGHTS'
printlog,logf,'------------------------'

;-----------------------------------
; Computs Weights
if not keyword_set(fake) then begin
  ALLFRAME_GETWEIGHTS,mchfile,weights,scales,sky ;,raw2=raw2
  invscales = 1.0/scales
  bdscale = where(scales lt 1e-5 or invscales gt 900,nbdscale)
  if nbdscale gt 0 then begin
    scales[bdscale] = 1.0
    invscales[bdscale] = 1.0
    weights[bdscale] = 0.0
  endif
  weightfile = mchbase+'.weights'
  WRITECOL,weightfile,weights,fmt='(F10.6)'
  scalefile = mchbase+'.scale'
  WRITECOL,scalefile,invscales,fmt='(F10.6)'  ; want to scale it UP
  zerofile = mchbase+'.zero'
  WRITECOL,zerofile,-sky,fmt='(F12.4)'  ; want to remove the background, set to 1st frame

; FAKE, use existing ones
endif else begin
  weightfile = mchbase+'.weights'
  scalefile = mchbase+'.scale'
  zerofile = mchbase+'.zero'
  READCOL,weightfile,weights,format='F',/silent
  READCOL,scalefile,invscales,format='F',/silent
  scales = 1.0/invscales
  READCOL,zerofile,sky,format='F',/silent
  sky = -sky
endelse
; Stuff the information into the FILEINFO structure
filestr.comb_weights = weights
filestr.comb_scale = invscales
filestr.comb_zero = -sky


;###########################################
; STEP 2: Resample the images

printlog,logf,'---------------------------'
printlog,logf,'STEP 2: Resample the images'
printlog,logf,'---------------------------'

CASE tile.type of

; --- WCS projection on the sky ---
'WCS': begin

  ; Convert x/y and ra/dec grid for entire tile image
  ;  that we'll be using, ~15sec
  xb = (lindgen(tile.nx)+tile.xrange[0])#replicate(1,tile.ny)
  yb = replicate(1,tile.nx)#(lindgen(tile.ny)+tile.yrange[0])
  HEAD_XYAD,tile.head,xb,yb,rab,decb,/deg

  ; Loop through the files
  For i=0,nfiles-1 do begin
    FITS_READ,filestr[i].fitsfile,im1,head1

    ; Make the mask
    mask = bytarr(filestr[i].nx,filestr[i].ny)+1
    gdpix = where(im1 lt filestr[i].saturate,ngdpix,comp=bdpix,ncomp=nbdpix)
    ; 0-bad, 1-good  
    if nbdpix gt 0 then begin
      mask[bdpix] = 0
      im1[bdpix] = filestr[i].background
    endif

    ; Get X/Y range for this image in the final coordinate system
    HEAD_ADXY,tile.head,filestr[i].vertices_ra,filestr[i].vertices_dec,vx,vy,/deg
    xout = [floor(min(vx))-2 > tile.xrange[0], ceil(max(vx))+2 < tile.xrange[1]]
    xoutrel = xout-tile.xrange[0] ; relative to xrange[0]
    nxout = xout[1]-xout[0]+1
    yout = [floor(min(vy))-2 > tile.yrange[0], ceil(max(vy))+2 < tile.yrange[1]]
    youtrel = yout-tile.yrange[0] ; relative to yrange[0]
    nyout = yout[1]-yout[0]+1
    rr = rab[xoutrel[0]:xoutrel[1],youtrel[0]:youtrel[1]]
    dd = decb[xoutrel[0]:xoutrel[1],youtrel[0]:youtrel[1]]
    ALLFRAME_ADXYINTERP,head1,rr,dd,xx,yy,nstep=10

    ; The x/y position to bilinear need to be in the original system, ~1sec
    rim = BILINEAR(im1,xx,yy,missing=filestr[i].background)
    rmask = BILINEAR(mask,xx,yy,missing=0)

    ; Contruct final image
    fim = fltarr(tile.nx,tile.ny)+filestr[i].saturate
    fim[xoutrel[0]:xoutrel[1],youtrel[0]:youtrel[1]] = rim
    fmask = bytarr(tile.nx,tile.ny)
    fmask[xoutrel[0]:xoutrel[1],youtrel[0]:youtrel[1]] = rmask

    ; Contruct the final header
    fhead = head1
    PUTAST,fhead,tile.ast
    sxaddpar,fhead,'NAXIS1',tile.nx
    sxaddpar,fhead,'NAXIS2',tile.ny
    sxaddpar,fhead,'BPM',filestr[i].resampmask
    printlog,logf,filestr[i].resampfile,' ['+strtrim(xout[0],2)+':'+strtrim(xout[1],2)+','+$
                  strtrim(yout[0],2)+':'+strtrim(yout[1],2)+']'
    MWRFITS,fim,filestr[i].resampfile,fhead,/create
    mhead = fhead
    sxaddpar,mhead,'BITPIX',8
    sxdelpar,mhead,'BPM'
    MWRFITS,fmask,filestr[i].resampmask,mhead,/create
   ; this takes about ~37-50 sec for a 2kx4k image.
   ;  now it takes ~5 sec for a 2kx4k image
  Endfor
end

; --- Pixel based ---
'PIXEL': begin

  ; Expand the images to sizes that will allow all of the shifts
  hd1 = headfits(fitsfiles[0])
  nx = sxpar(hd1,'NAXIS1')
  ny = sxpar(hd1,'NAXIS2')
  ;  left, down, right, up
  pix_expand = [ abs(floor(min(xshift)) < 0), abs(floor(min(yshift)) < 0), $
                 ceil(max(xshift)) > 0 , ceil(max(yshift)) > 0 ]
  outmaskfiles = FILE_DIRNAME(tempfits)+'/'+FILE_BASENAME(tempfits,'.fits')+'.mask.shft.fits'
  printlog,logf,'Expanding images by [',strjoin(strtrim(pix_expand,2),','),'] pixels'
  LOADMCH,shiftmch+'.mch',files2,trans2
  nxf = nx+pix_expand[0]+pix_expand[2]
  nyf = ny+pix_expand[1]+pix_expand[3]
  xx1 = lindgen(nx)#replicate(1,ny)
  yy1 = replicate(1,nx)#lindgen(ny)
  for i=0,nfiles-1 do begin
    ; Image
    FITS_READ,tempfits[i],tim,thead
    background = median(tim)
    out = trans_coo(xx1[*],yy1[*],reform(trans[i,*]))
    xx2 = xx1*0.
    yy2 = yy1*0.
    xx2[*] = reform(out[0,*]) + pix_expand[0]  ; shift to expanded grid
    yy2[*] = reform(out[1,*]) + pix_expand[1]  ; shift to expanded grid
    triangulate,xx2,yy2,tr,b  ; triangulate
    xout = lindgen(nxf)
    yout = lindgen(nyf)
    tim2 = TRIGRID(xx2,yy2,tim, tr, XOUT = xout, YOUT = yout, missing=background)
    ;tim2 = fltarr(nxf,nyf)+background  ; set out of bounds pixels to background
    ;tim2[pix_expand[0]:pix_expand[0]+nx-1,pix_expand[1]:pix_expand[1]+ny-1]=tim
    thead2 = thead
    sxaddpar,thead2,'NAXIS1',nxf
    sxaddpar,thead2,'NAXIS2',nyf
    sxaddpar,thead2,'CRPIX1',sxpar(thead2,'CRPIX1')+pix_expand[0]
    sxaddpar,thead2,'CRPIX2',sxpar(thead2,'CRPIX2')+pix_expand[1]
    MWRFITS,tim2,outfiles[i],thead2,/create
    ; Mask
    mfile = file_basename(tempfits[i],'.fits')+'.mask.fits'
    FITS_READ,mfile,mim,mhead
    mim2 = TRIGRID(xx2,yy2,mim, tr, XOUT = xout, YOUT = yout, missing=background)
    ;mim2 = fltarr(nxf,nyf)   ; out of bounds pixels set to 0=bad
    ;mim2[pix_expand[0]:pix_expand[0]+nx-1,pix_expand[1]:pix_expand[1]+ny-1]=mim
    mhead2 = mhead
    sxaddpar,mhead2,'NAXIS1',nxf
    sxaddpar,mhead2,'NAXIS2',nyf
    sxaddpar,mhead2,'CRPIX1',sxpar(mhead2,'CRPIX1')+pix_expand[0]
    sxaddpar,mhead2,'CRPIX2',sxpar(mhead2,'CRPIX2')+pix_expand[1]
    ;outmaskfile = FILE_DIRNAME(outfiles[i])+'/'+FILE_BASENAME(outfiles[i],'.fits')+'.mask.shft.fits'
    MWRFITS,mim2,outmaskfiles[i],mhead2,/create
  endfor
end
else: stop,tile.type+' not implemented yet'
ENDCASE


; Creating new MCH file for the combined file
print,'Deriving new transformation equations for the resampled coordinate system'
for i=0,nfiles-1 do begin

  ; Convert X/Y of this system into the combined reference frame
  ngridbin = 50
  nxgrid = filestr[i].nx / ngridbin
  nygrid = filestr[i].ny / ngridbin
  xgrid = (lindgen(nxgrid)*ngridbin)#replicate(1,nygrid)
  ygrid = replicate(1,nxgrid)#(lindgen(nygrid)*ngridbin)
  HEAD_XYAD,(*filestr[i].head),xgrid,ygrid,ragrid,decgrid,/deg
  HEAD_ADXY,tile.head,ragrid,decgrid,refxgrid,refygrid,/deg

  ; Now fit the transformation
  xdiff = refxgrid-xgrid
  ydiff = refygrid-ygrid
  xmed = median([xdiff],/even)
  ymed = median([ydiff],/even)
  ; Fit rotation with linear fits if enough points
  coef1 = robust_poly_fitq(ygrid,xdiff,1)  ; fit rotation term
  coef1b = dln_poly_fit(ygrid,xdiff,1,measure_errors=xdiff*0+0.1,sigma=coef1err,/bootstrap)
  coef2 = robust_poly_fitq(xgrid,ydiff,1)  ; fit rotation term
  coef2b = dln_poly_fit(xgrid,ydiff,1,measure_errors=ydiff*0+0.1,sigma=coef2err,/bootstrap)
  ;theta = mean([-coef1[1],coef2[1]])
  WMEANERR,[-coef1[1],coef2[1]],[coef1err[1],coef2err[1]],theta,thetaerr

  ; [xoff, yoff, cos(th), sin(th), -sin(th), cos(th)]
  ;trans = [xmed, ymed, 1.0, 0.0, 0.0, 1.0]
  trans = [xmed, ymed, 1.0-theta^2, theta, -theta, 1.0-theta^2]
  ; Adjust Xoff, Yoff with this transformation
  xyout = trans_coo(xgrid,ygrid,trans)
  trans[0] += median([refxgrid - xyout[0,*]],/even) 
  trans[1] += median([refygrid - xyout[1,*]],/even)

  ; Fit full six parameters if there are enough stars
  fa = {x1:(refxgrid)(*),y1:(refygrid)(*),x2:(xgrid)(*),y2:(ygrid)(*)}
  initpar = trans
  fpar = MPFIT('trans_coo_dev',initpar,functargs=fa, perror=perror, niter=iter, status=status,$
               bestnorm=chisq, dof=dof, autoderivative=1, /quiet) 
  trans = fpar

  diff = trans_coo_dev(fpar,x1=refxgrid,y1=refygrid,x2=xgrid,y2=ygrid)
  rms = sqrt(mean(diff^2.))
  filestr[i].resamptrans = trans
  filestr[i].resamptransrms = rms

  ; The output is:
  ; filename, xshift, yshift, 4 trans, FRAD (from als file), 0.0
  format = '(A2,A-30,A1,2F10.2,4F10.5,2F10.3)'
  newline = STRING("'",filestr[i].catfile,"'",trans, filestr[i].magoff[0], rms, format=format)
  PUSH,mchfinal,newline

  ; Printing the transformation
  printlog,logf,format='(A-20,2F10.4,4F12.8,2F10.3)',filestr[i].catfile,trans,filestr[i].magoff[0],rms
endfor
; Write to the new MCH file
combmch = mchbase+'_comb.mch'
WRITELINE,combmch,mchfinal


;###########################################
; STEP 5: COMBINE IMAGES
printlog,logf,'-------------------'
printlog,logf,'STEP 5: IMCOMBINE'
printlog,logf,'-------------------'

; The imcombine input file
resampfile = mchbase+'.resamp'
WRITELINE,resampfile,filestr.resampfile

; SCALE the images for combining
;-------------------------------
if not keyword_set(nocmbimscale) then begin

  ;; Put BPM mask names in file headers
  ;;  these will be used by IMCOMBINE
  ;for i=0,nfiles-1 do begin
  ;  head = headfits(outfiles[i])
  ;  sxaddpar,head,'BPM',outmaskfiles[i]
  ;  modfits,outfiles[i],0,head
  ;endfor

  ; Combine the frames WITH scaling/offset/masking, for the bright stars
  ;printlog,logf,'Creating SCALED image'
  combfile = mchbase+'_comb.fits'
  FILE_DELETE,combfile,/allow
  FILE_DELETE,mchbase+'_comb.bpm.pl',/allow
  IRAF_IMCOMBINE,'@'+resampfile,combfile,combine='average',reject='avsigclip',$
                 weight='@'+weightfile,rdnoise='!rdnoise',gain='!gain',$
                 irafdir=irafdir,error=imcombineerror2,scale='@'+scalefile,zero='@'+zerofile,$
                 masktype='badvalue',maskvalue=0,bpmasks=mchbase+'_comb.bpm'

  if n_elements(imcombineerror2) ne 0 then begin
    printlog,logf,'ERROR in IRAF_IMCOMBINE'
    printlog,logf,imcombineerror2
    error = imcombineerror2
    return
  endif

  ; Convert BPM mask from PL to FITS
  FILE_DELETE,mchbase+'_comb.bpm.fits',/allow
  undefine,lines
  cd,current=curdir
  push,lines,'cd '+curdir
  push,lines,'imcopy '+mchbase+'_comb.bpm.pl '+mchbase+'_comb.bpm.fits'
  push,lines,'logout'
  tempfile = mktemp('tiraf')
  WRITELINE,tempfile,lines
  IRAF_RUN,tempfile,irafdir,silent=silent,out=out,error=error

  ; Delete temporary scripts and PL file
  FILE_DELETE,[tempfile,mchbase+'_comb.bpm.pl'],/allow


  ; Fix the rdnoise and background/sky level and saturate
  ;  the bad pixels for DAOPHOT
  ;------------------------------------------------------

  ; 10/02/12
  ; THE IMAGES ARE (1) ZERO-SUBTRACTED, (2) SCALED, AND (3) WEIGHT AVERAGED
  ; The algorithm is:
  ; 1.) add zero-level correction.  im = im+zero
  ; 2.) scale the images.  im = im*scale
  ; 3.) take weighted average.  combim=total(weight*im)
  ;      there is also clipping that takes place during the averaging
  ; The final RDNOISE is essentially: comb_rdnoise = sqrt(total((weights*rdnoise*scale)^2))
  ; A gain that changes from frame to frame could be problematic,
  ; but this shouldn't happen since it's the same chip from the same night.

  ; IMCOMBINE wants rdnoise in electrons and gain in electrons/DN.
  ; DAOPHOT expects rdnoise in DN.  That's why mkopt converts
  ;  it with the gain.  So we are fine.  The header should have
  ;  rdnoise in ELECTRONS.

  ; page 63-65 of daophot2.pdf shows how you need to modify rdnoise/gain
  ; when averaging/summing frames. in observing/mosaic/.

  ; Load the IMCOMBINE output combined file and BPM
  FITS_READ,combfile,combim,combhead
  FITS_READ,mchbase+'_comb.bpm.fits',badmask,maskhead  ; 0-good, 1-bad


  ; Fix the gain
  ; For N averaged frames gain(N)=N*gain(1)
  ; Leave the gain as is!  We are scaling everything to the reference
  ; and using its gain.  It's nearly impossible to figure out the real
  ; gain since we are scaling the images and then taking a weighted
  ; average with outlier rejection.  Find a gain that properly
  ; describes/follows Poisson errors for the final combined image is
  ; difficult/impossible.  But that's okay.  This is just for source
  ; detection and DAOPHOT FIND just cares about the noise in the
  ; background.  We just need to ensure that the sky and rdnoise
  ; are correct.

  ; Fix the rdnoise
  ; The final RDNOISE is essentially: comb_rdnoise = sqrt(total((weights*rdnoise*scale)^2))
  rdnoisearr = fltarr(nfiles)
  for i=0,nfiles-1 do rdnoisearr[i] = PHOTRED_GETRDNOISE(base[i]+'.fits')
  ;  the "scales" array here is actually 1/scales used by IMCOMBINE.
  rdnoise = sqrt(total((weights*rdnoisearr/scales)^2))
  rdnoise = rdnoise > 0.01   ; must be >=0.01 or it will be 0.00 in the opt file and daophot will crash
  dummy = PHOTRED_GETRDNOISE(combfile,keyword=rdnoisekey) ; get keyword
  sxaddpar,combhead,rdnoisekey,rdnoise

  ; Fix the sky
  ; DAOPHOT FIND computes the random error per pixel in ADU as
  ; noise = sqrt( sky level/gain + rdnoise^2)
  ; So it assumes that the noise in the background is sqrt(sky/gain)
  ; in ADU.  We need to set the sky level so this is correct.
  ; The final noise should be 
  ; final sky noise = sqrt(total((weights*scale*sqrt(sky/gain))^2)) 
  ; So the final sky level should be
  ; final sky = total((weights*scale*sqrt(sky/gain))^2)*gain
  gain = PHOTRED_GETGAIN(combfile,keyword=gainkey)
  comb_sky = total((weights*sqrt((sky>0)/gain)/scales)^2)*gain
  ; the "scales" array here is actually 1/scale
  combim += float(comb_sky)  ; keep it float


  ; set the maximum to a "reasonable" level
  ; Rescale the image and increase the gain
  if max(combim) gt 50000 then begin
    rescale = 50000./max(combim)
    combim = combim*rescale
    sxaddpar,combhead,gainkey,gain/rescale
    ; rdnoise does NOT get modified since it's in electrons
    ; we just need to modify the gain which takes you from ADUs to electrons
  endif

  maskdatalevel = max(combim) + 10000       ; set "bad" data level above the highest "good" value
  combim2 = combim*(1-badmask) + maskdatalevel*badmask    ; set bad pixels to maskdatalevel
  MWRFITS,combim2,combfile,combhead,/create  ; fits_write can create an empty PDU

  ; Create the weight map for Sextractor using the BPM output by IMCOMBINE
  ;  bad only if bad in ALL images
  weightmap = -2.0*float(badmask eq 1) + 1.0
  combweightfile = mchbase+'_comb.mask.fits'
  MWRFITS,weightmap,combweightfile,whead,/create

; NO SCALING of the images for combining
;---------------------------------------
Endif else begin

  combfile = mchbase+'_comb.fits'
  FILE_DELETE,combfile,/allow
  IRAF_IMCOMBINE,'@'+resampfile,combfile,combine='average',reject='avsigclip',$
                 weight='@'+weightfile,rdnoise='!rdnoise',gain='!gain',$
                 irafdir=irafdir,error=imcombineerror2

  if n_elements(imcombineerror2) ne 0 then begin
    printlog,logf,'ERROR in IRAF_IMCOMBINE'
    printlog,logf,imcombineerror2
    error = imcombineerror2
    return
  endif

  ; Fix the rdnoise and background/sky level and saturate
  ;  the bad pixels for DAOPHOT
  ;------------------------------------------------------

  ; See the explanations for all these steps above!!

  ; Load the IMCOMBINE output combined file and BPM
  FITS_READ,combfile,combim,combhead
  FITS_READ,mchbase+'_comb.bpm.fits',badmask,maskhead  ; 0-good, 1-bad

  ; Fix the rdnoise
  ; The final RDNOISE is essentially: comb_rdnoise = sqrt(total((weights*rdnoise)^2))
  rdnoisearr = fltarr(nfiles)
  for i=0,nfiles-1 do rdnoisearr[i] = PHOTRED_GETRDNOISE(base[i]+'.fits')
  rdnoise = sqrt(total((weights*rdnoisearr)^2))
  dummy = PHOTRED_GETRDNOISE(combfile,keyword=rdnoisekey) ; get keyword
  sxaddpar,combhead,rdnoisekey,rdnoise

  ; Fix the sky
  ; So the final sky level should be
  ; final sky = total((weights*scale*sqrt(sky/gain))^2)*gain
  gain = PHOTRED_GETGAIN(combfile,keyword=gainkey)
  comb_sky = total((weights*sqrt(sky/gain))^2)*gain
  ; the "scales" array here is actually 1/scale
  combim += comb_sky


  ; set the maximum to a "reasonable" level
  ; Rescale the image and increase the gain
  if max(combim) gt 50000 then begin
    rescale = 50000./max(combim)
    combim = combim*rescale
    sxaddpar,combhead,gainkey,gain/rescale
    ; rdnoise does NOT get modified since it's in electrons
    ; we just need to modify the gain which takes you from ADUs to electrons
  endif


  ; Making Sextractor "weight" map file
  ;------------------------------------
  ; masks have 0-bad, 1-good.
  ; anything with less than 1.0 is considered bad
  ; weight map, -1 is bad, +1 is good
  ; "bpm" is the SUM of the bad pixel masks
  ; consider a pixel bad that is bad in ANY image
  weightmap = -2.0*float(bpm lt nfiles) + 1.
  combweightfile = mchbase+'_comb.mask.fits'
  FITS_WRITE,combweightfile,weightmap,whead

  ;---------------------------------------------
  ; SATURATE BAD pixels in the COMBINED IMAGE
  ; DAOPHOT needs to have the bad pixels "saturated",
  ; SExtractor will know which pixels are bad from the "weight" map.
  ;
  ; We could skip the fiximage.pro step but we still need the
  ; individual bpm masks and setting the bad pixels to the background
  ; probably helps in the IMALIGN/IMCOMBINE steps.
  printlog,logf,''
  printlog,logf,'"Saturating" bad pixels in the COMBINED image'
  printlog,logf,''


  badmask = float(weightmap lt 0.5)
  maskdatalevel = max(combim) + 10000       ; set "bad" data level above the highest "good" value
  combim2 = combim*(1.0-badmask) + maskdatalevel*badmask    ; set bad pixels to 100,000
  FITS_WRITE,combfile,combim2,combhead

Endelse ; no scaling of images for combining

; Add TILETYPE to the combined image
combhead = headfits(combfile)
sxaddpar,combhead,'AFTILTYP',tile.type
MODFITS,combfile,0,combhead

; Delete the resampled images
FILE_DELETE,filestr.resampfile,/allow,/quiet
FILE_DELETE,filestr.resampmask,/allow,/quiet

; Make the common source file for the combined image
;---------------------------------------------------
if keyword_set(usecmn) then begin
  print,'Combining COMMON SOURCE files for the combined image.'
  ; Loop through the files and convert to coordinates to the comined file
  undefine,allcmn
  for i=0,nfiles-1 do begin
    cmnfile1 = file_basename(filestr[i].fitsfile,'.fits')+'.cmn.lst'
    if file_test(cmnfile1) eq 1 then begin
      cmn1 = IMPORTASCII(cmnfile1,fieldnames=['id','x','y','mag','err','sky','skysig','sharp','round','round2'],$
                         skipline=3,/silent)
      READLINE,cmnfile1,cmnlines1
      coohead1 = cmnlines1[0:1]
      ncmn1 = n_elements(cmn1)
      ; Get coordinates on the resampled/combined image grid
      out = trans_coo(cmn1.x,cmn1.y,filestr[i].resamptrans)
      newx = reform(out[0,*])
      newy = reform(out[1,*])
      cmn1.x = newx
      cmn1.y = newy
      if n_elements(allcmn) eq 0 then begin
        allcmn = cmn1
      endif else begin
        ; Remove any duplicates
        SRCMATCH,allcmn.x,allcmn.y,cmn1.x,cmn1.y,2.0,ind1,ind2,count=nmatch
        if nmatch gt 0 then begin
          if nmatch lt ncmn1 then remove,ind2,cmn1 else undefine,cmn1
        endif
        if n_elements(cmn1) gt 0 then push,allcmn,cmn1
      endelse
    endif
  endfor
  if n_elements(allcmn) gt 0 then begin
    WRITECOL,mchbase+'_comb.cmn.lst',allcmn.id,allcmn.x,allcmn.y,allcmn.mag,allcmn.err,allcmn.sky,allcmn.skysig,$
             allcmn.sharp,allcmn.round,allcmn.round2,fmt='(I7,2F9.2,3F9.3,F9.2,3F9.3)'
    WRITELINE,mchbase+'_comb.cmn.lst',[coohead1,''],/prepend  ; prepend the COO header
  endif else printlog,logf,'No common file to combine'
; Don't use common file 
endif else begin
  ; Make sure it doesn't exist otherwise it will be used
  FILE_DELETE,mchbase+'_comb.cmn.lst',/allow
endelse

BOMB:

if keyword_set(stp) then stop

end
