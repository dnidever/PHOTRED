;+
;
; PHOTRED_MAKETILES
;
; This program figures out the tiling scheme for a field.
;
; INPUTS:
;  fieldinput  The directory name and field name combined, e.g. "20100610/F1"
;  thisimager  The structure with information on the imager.
;  =pixscale   The pixel scale to use for the projection.  The default
;                is to use the mean pixel scale of all images rounded up to the
;                nearest 100th of a pixel.
;  =tilesize   The tile size (square) in pixels.  The default is 5,000 pixels.
;  =logfile    The name of the log file.
;  /stp        Stop at the end of the program.
; 
; OUTPUTS:
;  The tiling scheme will be written to a file called fieldinput+".tiling".
;  =tilehead   The header with the tiling scheme projection information.
;  =tilestr    The structure with information on each tile.
;  =error      The error message if one occurred.
;
; USAGE:
;  IDL>photred_maketiles,'F1',thisimager,logfile=logfile
;
; D.Nidever  Jan 2017
;-

pro photred_maketiles,fieldinput,thisimager,pixscale=inppixscale,tilesize=inptilesize,$
                      tilehead=tilehead,tilestr=tilestr,logfile=logfile,error=error,stp=stp

undefine,tilehead
undefine,tilestr
undefine,error

; Not enough inputs
if n_elements(fieldinput) eq 0 or n_elements(thisimager) eq 0 then begin
  print,'Syntax - photred_maketiles,fieldinput,thisimager,pixscale=pixscale,tilesize=tilesize,'
  print,'                           logfile=logfile,stp=stp,tilehead=tilehead,tilestr=tilestr,error=error'
  error = 'Not enough inputs'
  return
endif

; Logfile
if n_elements(logfile) eq 0 then logfile=-1

; Break "fieldinput" into directory and field name components
dir = FILE_DIRNAME(fieldinput)
cd,current=curdir
cd,dir
field = FILE_BASENAME(fieldinput)
  
printlog,logfile,'Generating TILING scheme for Field >>'+field+'<< in directory '+dir

;; Find all existing FITS files for this group
if thisimager.namps gt 1 then $
  fieldfiles = FILE_SEARCH(dir+'/'+field+'-*'+thisimager.separator+'*.fits',count=nfieldfiles) else $
  fieldfiles = FILE_SEARCH(dir+'/'+field+'-*.fits',count=nfieldfiles)

; Remove a.fits, s.fits, _comb.fits and other "temporary" files.
if nfieldfiles gt 0 then begin
  fbases = FILE_BASENAME(fieldfiles,'.fits')
  bad = where(stregex(fbases,'a$',/boolean) eq 1 or $ ; psf stars image
              stregex(fbases,'s$',/boolean) eq 1 or $ ; allstar subtracted file
              stregex(fbases,'_0$',/boolean) eq 1 or $ ; _0 head file from split
              stregex(fbases,'_comb$',/boolean) eq 1 or $ ; stacked field image
              stregex(fbases,'_comb.bpm$',/boolean) eq 1 or $ ; stacked field image
              stregex(fbases,'_comb_sub$',/boolean) eq 1 or $ ; allstar subtracted stacked image
              stregex(fbases,'j$',/boolean) eq 1 or $         ; allframe temp file
              stregex(fbases,'k$',/boolean) eq 1 or $         ; allframe temp file
              stregex(fbases,'jnk$',/boolean) eq 1,nbad)      ; daophot? temp file
  if nbad gt 0 then begin
    if nbad eq nfieldfiles then begin
      undefine,fieldfiles
      nfieldfiles = 0
    endif else begin
      REMOVE,bad,fieldfiles
      nfieldfiles = n_elements(fieldfiles)
    endelse
  endif                      ; some ones to remove
endif else begin              ; some fieldfiles
  error = 'No '+field+' files found in current directory'
  printlog,logfile,error
  return
endelse
printlog,logfile,strtrim(nfieldfiles,2)+' FITS files found for field '+field

;; Loop through the FITS files and get their information
printlog,logfile,'Gathering information on the files'
PHOTRED_GATHERFILEINFO,fieldfiles,filestr
; Only want good ones
gdfilestr = where(filestr.exists eq 1,ngdfilestr)
if ngdfilestr eq 0 then begin
  error = 'None of the files exist'
  printlog,logfile,error
  return
endif
filestr0 = filestr
filestr = filestr[gdfilestr]  ; only want the ones that exist

; Creating the overall tiling projection and scheme
;--------------------------------------------------
; this came from allframe_combine.pro
printlog,logfile,'Creating TILING scheme'
; The default projection is a tangent plane centered
; at halfway between the ra/dec min/max of all of
; the images.  The mean pixel scale is used.
;  near RA=0 line
if range(filestr.vertices_ra) gt 180 then begin
   vertices_ra = filestr.vertices_ra
   over = where(vertices_ra gt 180,nover,comp=under,ncomp=nunder)
   if nover gt 0 then vertices_ra[over]-=360
   rar = minmax(vertices_ra)
   cenra = mean(rar)
endif else begin
   rar = minmax(filestr.vertices_ra)
   cenra = mean(rar)
endelse
; RA/DEC ranges
decr = minmax(filestr.vertices_dec)
cendec = mean(decr)
; Pixel scale
if n_elements(inppixscale) gt 0 then begin
  pixscale = inppixscale
endif else begin
  pixscale = mean(filestr.pixscale)
  ; round up to nearest 100th
  pixscale = ceil(pixscale*100.)/100.
endelse
; Set up the tangent plane projection
step = pixscale/3600.0d0
delta_dec = range(decr)
delta_ra = range(rar)*cos(cendec/!radeg)
; Make the full tiling scheme extend slightly more than
;  full extent of the images
nx = ceil(delta_ra*1.01/step)
ny = ceil(delta_dec*1.01/step)
xref = nx/2
yref = ny/2
; Tile size in pixels
;  the default is 10,000 pixels
if n_elements(inptilesize) gt 0 then tilesize=inptilesize[0] else tilesize=5000L

; Make header with the WCS
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

; Number of tiling columns and rows
nxtile = ceil(nx/float(tilesize))
nytile = ceil(ny/float(tilesize))

; Print out the tiling information
printlog,logfile,'RA range = ['+strtrim(rar[0],2)+','+strtrim(rar[1],2)+'] deg'
printlog,logfile,'DEC range = ['+strtrim(decr[0],2)+','+strtrim(decr[1],2)+'] deg'
printlog,logfile,'Central RA = '+strtrim(cenra,2)
printlog,logfile,'Central DEC = '+strtrim(cendec,2)
printlog,logfile,'NX = '+strtrim(nx,2)
printlog,logfile,'NY = '+strtrim(ny,2)
printlog,logfile,'Tilesize = '+strtrim(tilesize,2)
printlog,logfile,'Tiles = '+strtrim(nxtile,2)+'x'+strtrim(nytile,2)


; Get information for each tile
tilestr = replicate({num:0L,name:'',x0:0L,x1:0L,nx:0L,y0:0L,y1:0L,ny:0L,nimages:0L},nxtile*nytile)
tilecnt = 0LL
for i=0,nxtile-1 do begin
  x0 = i*tilesize
  x1 = (x0+tilesize-1) < (nx-1)
  for j=0,nytile-1 do begin
    y0 = j*tilesize
    y1 = (y0+tilesize-1) < (ny-1)
    tilestr[tilecnt].num = tilecnt+1
    tilestr[tilecnt].x0 = x0
    tilestr[tilecnt].x1 = x1
    tilestr[tilecnt].nx = x1-x0+1
    tilestr[tilecnt].y0 = y0
    tilestr[tilecnt].y1 = y1
    tilestr[tilecnt].ny = y1-y0+1
    tilecnt++
  endfor
endfor
ntiles = n_elements(tilestr)

;; Figure out which images overlap the tiles
for i=0,nfieldfiles-1 do begin
  HEAD_ADXY,tilehead,filestr[i].vertices_ra,filestr[i].vertices_dec,fx,fy,/deg
  ;; Loop over the tiles
  for j=0,ntiles-1 do begin
    tx = [tilestr[j].x0,tilestr[j].x1,tilestr[j].x1,tilestr[j].x0]
    ty = [tilestr[j].y0,tilestr[j].y0,tilestr[j].y1,tilestr[j].y1] 
    overlap = DOPOLYGONSOVERLAP(fx,fy,tx,ty)
    if overlap eq 1 then tilestr[j].nimages++
 endfor  
endfor

; Only keep tiles that have images that overlap
gdtiles = where(tilestr.nimages gt 0,ngdtiles)
printlog,logfile,strtrim(ngdtiles,2)+' tiles have overlapping images'
tilestr0 = tilestr
tilestr = tilestr[gdtiles]
ntiles = ngdtiles
; Make final IDs and names
tilestr.num = lindgen(ntiles)+1
tilestr.name = field+'-T'+strtrim(tilestr.num,2)

; Make the tiling file
;---------------------
undefine,lines
; Lines with the tiling scheme first
push,lines,'CENRA  = '+strtrim(cenra,2)
push,lines,'NX     = '+strtrim(nx,2)
push,lines,'XSTEP  = '+strtrim(step,2)
push,lines,'XREF   = '+strtrim(xref+1,2)
push,lines,'CENDEC = '+strtrim(cendec,2)
push,lines,'NY     = '+strtrim(ny,2)
push,lines,'YSTEP  = '+strtrim(step,2)
push,lines,'YREF   = '+strtrim(yref+1,2)
push,lines,'NTILES = '+strtrim(ntiles,2)
; Then one file with info for each tile
for i=0,ntiles-1 do begin
  tilestr1 = tilestr[i]
  fmt = '(I-6,A8,6I8,I6)'
  ; Using IRAF indexing
  line = string(format=fmt,tilestr1.num,tilestr1.name,tilestr1.x0+1,tilestr1.x1+1,tilestr1.nx,$
                tilestr1.y0+1,tilestr1.y1+1,tilestr1.ny,tilestr1.nimages)
  push,lines,line
endfor
tilefile = dir+'/'+field+'.tiling'
printlog,logfile,'Writing tiling information to >>'+tilefile+'<<'
WRITELINE,tilefile,lines

if keyword_set(stp) then stop

end
