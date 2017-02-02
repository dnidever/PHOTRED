pro addfakes

; Add artificial stars to a set of images for a given field

dir = '/datalab/users/dnidever/smash/cp/red/photred/addfakes/Field100/'

; Get list of artifical stars with the right bands
synth = IMPORTASCII(dir+'CMD_field100',/header)

; For Thomas' stars use the stellar locus to get the other
; three bands.

; Load the photometric transformation information
transfile = smashred_rootdir()+'cp/red/photred/stdred/smashred_transphot_eqns.fits'
trans_fitstr = MRDFITS(transfile,1,/silent)
trans_chipstr = MRDFITS(transfile,2,/silent)
trans_ntstr = MRDFITS(transfile,3,/silent)

; Get the images
files = file_search(dir+'F2/F2-*_??.fits',count=nfiles)
; Get the info for the files
PHOTRED_GATHERFILEINFO,files,filestr
filestr.filter = strmid(filestr.filter,0,1)
; Get spatial coverage of the data


; Create RA/DEC positions for the stars
stop

; Image loop
For i=0,nfiles-1 do begin

  ; Find artificial stars within the boundary of this image

  ; Use WCS to convert RA/DEC to X/Y and do final spatial cut

  ; Use photometric transformation equations to get instrumental
  ; magnitudes in this image (filter, exptime).

  ; Use daophot_addstar. pro to add artificial stars

Endfor 


stop

end
