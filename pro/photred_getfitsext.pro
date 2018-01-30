;+
;
; PHOTRED_GETFITSEXT
;
; Gets the fits extension (.fits or .fits.fz) and optionall
; the base name.
;
; INPUTS:
;  files      Scalar or array of file names.
;  /isfpack   Returns a boolean 1 if the file is a fpack compressed
;               FITS file and 0 if it.
;  /basename  Returns the basename only.
;  /full      Returns an array of basename and extension.  If a
;               single file is input then the output will a 2-element
;               array, otherwise it is [Nfiles,2].
;
; OUTPUTS:
;  output     The output of the program which depends on the inputs.
;               By default this is the FITS extension (.fits or .fits.fz)
;               and '' if none of those.
;               /isfpack   Returns a boolean 1 or 0 if the file(s) are
;                            fits.fz extensions.
;               /basename  Returns the file basename (without the
;                            extension or directory).
;               /full      Returns the file basename and extension in
;                            2-element or [Nfiles,2] array depending
;                            on how many files were input.
; 
;
; USAGE:
;  IDL>ext=photred_getfitsext('F2-04958679_01.fits.fz')
;
; By D. Nidever  August 2017
;-

function photred_getfitsext,files,isfpack=isfpack,basename=basename,full=full

; Not enough inputs
nfiles = n_elements(files)
if nfiles eq 0 then begin
  print,'Syntax - out = photred_getfitsext(files,/isfpack,/basename,/full)'
  return,-1
endif

; Cases
; 1 - extension only
; 2 - isfpack, boolean 1/0
; 3 - basename only
; 4 - full, basename and extension
usecase = 1                              ; extension by default
if keyword_set(isfpack) then usecase=2
if keyword_set(basename) then usecase=3
if keyword_set(full) then usecase=4

; Initializing the output array
undefine,out
Case usecase of
  ; Extension only
  1: out=strarr(nfiles)
  ; isfpack, boolean 1/0
  2: out=intarr(nfiles)
  ; basename only
  3: out=strarr(nfiles)
  ; full, basename and extension
  4: if nfiles eq 1 then out=strarr(2) else out=strarr(nfiles,2)
  else: ; Not supported
Case

; Loop through the files
For i=0,nfiles-1 do begin
  file = file_basename(strtrim(files[i],2))

  ; Get the extension
  ext = ''
  if strmid(file,6,7,/reverse_offset) eq 'fits.fz' then ext='.fits.fz'
  if strmid(file,3,4,/reverse_offset) eq 'fits' then ext='.fits'

  Case usecase of
    ; Extension only
    1: out[i] = ext
    ; isfpack, boolean 1/0
    2: if ext eq '.fits.fz' then out[i]=1
    ; basename only
    3: out[i] = file_basename(file,ext)
    ; full, basename and extension
    4: begin
         if nfiles eq 1 then begin
           out[0] = file_basename(file,ext)
           out[1] = ext
         endif else begin
           out[i,0] = file_basename(file,ext)
           out[i,1] = ext
         endelse
       end
    else: ; Not supported
  Endcase
Endfor

stop

return,out

end
