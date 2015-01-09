pro fiximages,input,satlevel=satlevel,nofix=nofix,stp=stp

;+
;
; FIXIMAGES
;
; This programs fixes bad pixels in images.
;
; INPUTS:
;  input      The filenames
;  =satlevel  The saturation level to use if NOT
;               in the FITS header.
;  /nofix     Don't fix the images.
;  /stp       Stop at the end of the program
;
; OUTPUTS:
;  The files are overwritten with the fixed images.
;  
; USAGE:
;  IDL>fiximages,'*.fits',satlevel=5e4
;
; By D.Nidever   August 2008
;-

; Not enough inputs
ninput = n_elements(input)
if ninput eq 0 then begin
  print,'Syntax - fiximages,input,satlevel=satlevel,stp=stp'
  return
endif

; Load the input
LOADINPUT,input,files,count=nfiles

; No files to process
if nfiles eq 0 then begin
  print,'No files to process'
  return
endif

; Loop through the files
For i=0,nfiles-1 do begin

  ifile = files[i]
  ibase = FILE_BASENAME(ifile,'.fits')
  idir = FILE_DIRNAME(ifile)

  ; Does the file exist
  test = FILE_TEST(ifile)
  if (test eq 0) then begin
    print,ifile,' NOT FOUND'
    goto,BOMB
  endif

  ; Read the FITS file
  undefine,im,head,error
  FITS_READ,ifile,im,head,message=error,/no_abort

  ; There a problem opening the file
  if error ne '' then begin
    print,'ERROR opening ',ifile
    print,error
  endif

  ; Add mask name to header
  maskname = idir+'/'+ibase+'.mask.fits'
  sxaddpar,head,'BPM',maskname,' BPM mask file'

  ; Do we have bad pixels?
  saturate = SXPAR(head,'SATURATE',count=nsaturate)
  if nsaturate eq 0 then begin
    if keyword_set(satlevel) then if satlevel gt 0 then saturate=satlevel
    if n_elements(saturate) eq 0 then saturate=50000.
  endif
  bd = where(im gt saturate,nbd)

  ; Fix bad pixels
  if (nbd gt 0) then begin

    ; Background value
    backg = median(im)

    ; Set all bad pixels to the background value
    first = im
    first[bd] = backg

    ; Do one smoothing
    smlen = 3 ;5
    sm = smooth(first,smlen,/nan,/edge_truncate)

    ; Use the convolved image for the bad pixels
    newim = im
    newim[bd] = sm[bd]


    ; Print
    print,ifile,'  ',strtrim(nbd,2),' bad pixels fixed'

    ; Write the output
    outname = ifile
    if not keyword_set(nofix) then $
      ;FITS_WRITE,outname,newim,head
      MWRFITS,newim,outname,head,/create,/silent

    ; Make the bad pixel mask image, 0-bad, 1-good
    mask = im*0.+1.0
    mask[bd] = 0.0
    ;maskname = idir+'/'+ibase+'.mask.fits'
    ;FITS_WRITE,maskname,mask,head
    MWRFITS,mask,maskname,head,/create,/silent

  ; NO bad pixels
  Endif else begin
    print,ifile,' has no bad pixels'

    ; Make the bad pixel mask image, 0-bad, 1-good
    mask = im*0.+1.0
    ;maskname = idir+'/'+ibase+'.mask.fits'
    ;FITS_WRITE,maskname,mask,head
    MWRFITS,mask,maskname,head,/create,/silent

  Endelse

  BOMB:

end


if keyword_set(stp) then stop

end
