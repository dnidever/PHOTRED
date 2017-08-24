;+
;
; PHOTRED_GATHERFILEINFO
;
; Gather information about files like their corner coordinates, etc.
; This is used by the various tiling-related programs.
;
; INPUTS:
;  file      The FITS file names.
;
; OUTPUTS:
;  filestr   The structure with information for each file.
;  =error    The error message if one occurred.
;
; USAGE:
;  IDL>photred_gatherfileinfo,file,filestr
;
; By D.Nidever  Jan. 2017
;-

pro photred_gatherfileinfo,files,filestr,error=error

undefine,filestr,error

; Not enough inputs
if n_elements(files) eq 0 then begin
  print,'Syntax - photred_gatherfileinfo,files,filestr,error=error'
  error = 'Not enough inputs'
  return
endif

nfiles = n_elements(files)
; Create structure
filestr = replicate({file:'',exists:0,nx:0L,ny:0L,filter:'',exptime:0.0,dateobs:'',pixscale:0.0,cenra:0.0d0,cendec:0.0d0,$
                  vertices_ra:dblarr(4),vertices_dec:dblarr(4)},nfiles)
filestr.file = files
; File loop
For i=0,nfiles-1 do begin
  info = FILE_INFO(files[i])
  filestr[i].exists = info.exists
  if filestr[i].exists eq 0 then goto,BOMB
  if strmid(files[i],6,7,/reverse_offset) eq 'fits.fz' then head=HEADFITS(files[i],exten=1) else $
    head = HEADFITS(files[i])
  filestr[i].file = files[i]
  filestr[i].filter = photred_getfilter(files[i],/noupdate,/silent,error=filterr)
    if n_elements(filterr) gt 0 then filestr[i].filter=sxpar(head,'filter')
  filestr[i].exptime = sxpar(head,'exptime')
  filestr[i].dateobs = sxpar(head,'date-obs')
  filestr[i].nx = sxpar(head,'NAXIS1')
  filestr[i].ny = sxpar(head,'NAXIS2')
  ;GETPIXSCALE,'',pixscale,head=head
  HEAD_XYAD,head,[0,1]+filestr[i].nx/2,[0,0]+filestr[i].ny/2,ra1,dec1,/deg  ; a little bit faster  
  pixscale = sphdist(ra1[0],dec1[0],ra1[1],dec1[1],/deg)*3600d0
  filestr[i].pixscale = pixscale
  HEAD_XYAD,head,filestr[i].nx/2,filestr[i].ny/2,cenra1,cendec1,/degree
  filestr[i].cenra = cenra1
  filestr[i].cendec = cendec1
  HEAD_XYAD,head,[0,filestr[i].nx-1,filestr[i].nx-1,0],[0,0,filestr[i].ny-1,filestr[i].ny-1],vra,vdec,/degree
  filestr[i].vertices_ra = vra
  filestr[i].vertices_dec = vdec
  ; What do we do if there's no WCS???
  BOMB:
Endfor

;stop

end
