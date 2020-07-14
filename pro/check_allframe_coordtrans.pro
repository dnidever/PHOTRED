;+
;
; CHECK_ALLFRAME_COORDTRANS
;
; Check the coordinates transformations used by ALLFRAME.
; This uses the combined .mch file, the .nmg file, the FITS headers
; and the .alf files.
;
; INPUTS:
;  mchfile  The combined mch filename (absolute if not in that directory).
;  /silent  No output to the screen.
;
; OUTPUTS:
;  out   Structure giving information for each input image file.
;
; USAGE:
;  IDL>out = check_allframe_coordtrans(mchfile)
;
; By D. Nidever  June 2020
;;-

function check_allframe_coordtrans,mchfile,silent=silent

;; Check the ALLFRAME coordinate transformation

;; Not enough inputs
if n_elements(mchfile) eq 0 then begin
  print,'Syntax - out = check_allframe_coordtrans(mchfile)'
  return,-1
endif

if file_test(mchfile) eq 0 then begin
  print,mchfile,' NOT FOUND'
  return,-1
endif

dir = file_dirname(mchfile)
base = file_basename(mchfile,'.mch')

LOADMCH,mchfile,files,trans
LOADALS,dir+'/'+base+'.nmg',nmg
n = n_elements(files)
head = photred_readfile(dir+'/'+base+'.fits.fz',exten=1,/header)
HEAD_XYAD,head,nmg.x-1,nmg.y-1,ra,dec,/deg
out = replicate({file:'',nalf:0L,decstd:-1.0},n)
for i=1,n-1 do begin
  alffile = dir+'/'+repstr(files[i],'.als','.alf')
  out[i].file = alffile
  LOADALS,alffile,alf,count=nalf
  out[i].nalf = nalf
  if nalf gt 0 then begin
    fitsfile = dir+'/'+repstr(files[i],'.als','.fits')
    head1 = photred_readfile(fitsfile,/header)
    HEAD_XYAD,head1,alf.x-1,alf.y-1,ra1,dec1,/deg
    MATCH,nmg.id,alf.id,ind1,ind2,/sort,count=nmatch
    if nmatch gt 0 then out[i].decstd = stddev(dec[ind1]-dec1[ind2])*3600
    if not keyword_set(silent) then print,alffile,out[i].nalf,out[i].decstd
  endif
endfor

return,out

end
