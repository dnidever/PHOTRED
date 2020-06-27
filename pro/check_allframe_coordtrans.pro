function check_allframe_coordtrans,mchfile

;; Check the ALLFRAME coordinate transformation

dir = file_dirname(mchfile)
base = file_basename(mchfile,'.mch')

loadmch,mchfile,files,trans
loadals,dir+'/'+base+'.nmg',nmg
n = n_elements(files)
head = photred_readfile(dir+'/'+base+'.fits.fz',exten=1,/header)
head_xyad,head,nmg.x-1,nmg.y-1,ra,dec,/deg
out = replicate({file:'',nalf:0L,decstd:-1.0},n)
for i=1,n-1 do begin
  alffile = dir+'/'+repstr(files[i],'.als','.alf')
  out[i].file = alffile
  loadals,alffile,alf,count=nalf
  out[i].nalf = nalf
  if nalf gt 0 then begin
    fitsfile = dir+'/'+repstr(files[i],'.als','.fits')
    head1 = photred_readfile(fitsfile,/header)
    head_xyad,head1,alf.x-1,alf.y-1,ra1,dec1,/deg
    match,nmg.id,alf.id,ind1,ind2,/sort,count=nmatch
    if nmatch gt 0 then begin
      out[i].decstd = stddev(dec[ind1]-dec1[ind2])*3600
    endif
    print,alffile,out[i].nalf,out[i].decstd
  endif
endfor

return,out

end
