;+
;
; PHOTRED_MAKETILES
;
; This program figures out the tiling scheme for a field.
;
; INPUTS:
;  fieldinput  The directory name and field name combined, e.g. "20100610/F1"
;  thisimager  The structure with information on the imager.
;  =logfile    The name of the log file.
; 
; OUTPUTS:
;  The tiling scheme will be written to a file called fieldinput+".tiling".
;
; USAGE:
;  IDL>photred_maketiles,'F1',thisimager,logfile=logfile
;
; D.Nidever  Jan 2017
;-

pro photred_maketiles,fieldinput,thisimager,logfile=logfile

; Not enough inputs
if n_elements(fieldinput) eq 0 or n_elements(thisimager) eq 0 then begin
  print,'Syntax - photred_maketiles,fieldinput,thisimager,logfile=logfile'
  return
endif

; Break "fieldinput" into directory and field name components
dir = FILE_DIRNAME(fieldinput)
cd,current=curdir
cd,dir
field = FILE_BASENAME(fieldinput)
  
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
  printlog,logfile,'No ',field,' files found in current directory'
  return
endelse

;; Loop through the FITS files and get their information
filestr = replicate({field:'',file:'',nx:0L,ny:0L,filter:'',pixscale:0.0,cenra:0.0d0,cendec:0.0d0,$
                  vertices_ra:dblarr(4),vertices_dec:dblarr(4)},nfieldfiles)
For i=0,nfieldfiles-1 do begin
  head = HEADFITS(fieldfiles[i])
  filestr[i].field = field
  filestr[i].file = fieldfiles[i]
  filestr[i].filter = photred_getfilter(fieldfiles[i],/noupdate)
  filestr[i].nx = sxpar(head,'NAXIS1')
  filestr[i].ny = sxpar(head,'NAXIS2')
  GETPIXSCALE,'',pixscale,head=head
  filestr[i].pixscale = pixscale
  HEAD_XYAD,head,nx/2,ny/2,cenra,cendec
  filestr[i].cenra = cenra
  filestr[i].cendec = cendec
  HEAD_XYAD,head,[0,filestr[i].nx-1,filestr[i].nx-1,0],[0,0,filestr[i].ny-1,filestr[i].ny-1],vra,vdec,/degree
  filestr[i].vertices_ra = vra
  filestr[i].vertices_dec = vdec
  ; What do we do if there's no WCS???
Endfor

; Set central RA/DEC, based on 
; Deal with wraparound of RA=0

; this came from allframe_combine.pro
printlog,logf,'Creating TILE projection'
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
decr = minmax(filestr.vertices_dec)
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



stop

end
