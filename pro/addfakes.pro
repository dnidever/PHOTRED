pro addfakes,field,dir=dir,fakefield=fakefield

; Add artificial stars to a set of images for a given field

; Not enough inputs
if n_elements(field) eq 0 then begin
  print,'Syntax - addfakes,field,dir=dir,fakefield=fakefield'
  return
endif

;fakedir = '/datalab/users/dnidever/smash/cp/red/photred/addfakes/'
rootdir = smashred_rootdir()+'cp/red/photred/'
fakedir = rootdir+'addfakes/'

; Defaults
if n_elements(dir) eq 0 then begin
  cd,current=dir
  dir = dir+'/'
endif
;field = 'F2'
if n_elements(fakefield) eq 0 then fakefield = 'T1'

; Get the "global" field name from the "fields" file
readcol,dir+'fields',shname,lname,format='A,A',/silent
gshname = where(shname eq field,ngshname)
globalfield = lname[gshname[0]]

print,'Adding artificial stars to ',globalfield,' files in ',dir
print,'----------------------------------------------------------------'

; Get list of artifical stars with the right bands
;synth0 = IMPORTASCII(dir+'CMD_field100',/header)
;nsynth = n_elements(synth0)
;synth = replicate({id:0L,photid:0L,ra:0.0d0,dec:0.0d0,lon:0.0d0,lat:0.0d0,$
;                   u:0.0,uerr:0.0,g:0.0,gerr:0.0,r:0.0,rerr:0.0,$
;                   i:0.0,ierr:0.0,z:0.0,zerr:0.0},nsynth)
;synth.photid = lindgen(nsynth)+1
;synth.g = synth0.magnitude
;synth.i = synth0.magnitude - synth0.colour

; Make a "sheet" of artificial stars in CMD
gi_range = [-1.0,3.5]
g_range = [17.0,27.0]
nsynth = 1e5
synth = replicate({id:0L,photid:0L,ra:0.0d0,dec:0.0d0,lon:0.0d0,lat:0.0d0,$
                   u:0.0,uerr:0.0,g:0.0,gerr:0.0,r:0.0,rerr:0.0,$
                   i:0.0,ierr:0.0,z:0.0,zerr:0.0},nsynth)
synth.photid = lindgen(nsynth)+1
synth_gi = randomu(seed,nsynth)*range(gi_range)+gi_range[0]
synth.g = randomu(seed,nsynth)*range(g_range)+g_range[0]
synth.i = synth.g - synth_gi

; Use the stellar locus to get the other three bands.
restore,fakedir+'Field100_stellarlocus.dat'
; interpolate the stellar locus to get the other bands
INTERP,tstr.gibin,tstr.ui,synth.g-synth.i,ui
synth.u = ui + synth.i
INTERP,tstr.gibin,tstr.ri,synth.g-synth.i,ri
synth.r = ri + synth.i
INTERP,tstr.gibin,tstr.zi,synth.g-synth.i,zi
synth.z = zi + synth.i

; Get the images
files = file_search(dir+field+'/'+field+'-*_??.fits',count=nfiles)
if nfiles eq 0 then begin
  fzfiles = file_search(dir+field+'/'+field+'-*_??.fits.fz',count=nfzfiles)
  ; Need to uncompress the files
  if nfzfiles gt 0 then begin
    print,'Found ',strtrim(nfzfiles,2),' fits.fz files.  Uncompressing.'
    for i=0,nfzfiles-1 do spawn,['funpack',fzfiles[i]],/noshell
    files = file_dirname(fzfiles)+'/'+file_basename(fzfiles,'.fz')
    nfiles = n_elements(files)
  endif
endif

print,'Found ',strtrim(nfiles,2),' files to add artificial stars to'
if nfiles eq 0 then begin
  print,'No files found'
  return
endif
; Get the info for the files
PHOTRED_GATHERFILEINFO,files,filestr
add_tag,filestr,'dir','',filestr
filestr.dir = file_dirname(filestr.file)
add_tag,filestr,'base','',filestr
filestr.base = file_basename(filestr.file,'.fits')
add_tag,filestr,'field','',filestr
add_tag,filestr,'expnum','',filestr
add_tag,filestr,'chip',0L,filestr
for i=0,nfiles-1 do begin
  dum1 = strsplit(filestr[i].base,'-',/extract)
  filestr[i].field = dum1[0]
  dum2 = strsplit(dum1[1],'_',/extract)
  filestr[i].expnum = first_el(dum2)
  filestr[i].chip = long(first_el(dum2,/last))
endfor
filestr.filter = strmid(filestr.filter,0,1)
add_tag,filestr,'mjd',0L,filestr
for i=0,nfiles-1 do filestr[i].mjd=PHOTRED_GETMJD('','ctio',dateobs=filestr[i].dateobs)
add_tag,filestr,'airmass',0.0,filestr
for i=0,nfiles-1 do begin
  head = headfits(filestr[i].file)
  filestr[i].airmass = sxpar(head,'airmass')
endfor
; CAN'T I GET ALL OF THIS FROM THE FINAL PHOTRED SUMMARY FILE???

; Load the photometric transformation information
transfile = smashred_rootdir()+'cp/red/photred/stdred/smashred_transphot_eqns.fits'
trans_fitstr = MRDFITS(transfile,1,/silent)
trans_chipstr = MRDFITS(transfile,2,/silent)
trans_ntstr = MRDFITS(transfile,3,/silent)
trans_mjdchipfilt = strtrim(trans_fitstr.mjd,2)+'-'+strtrim(trans_fitstr.chip,2)+'-'+strtrim(trans_fitstr.filter,2)

; Load the CHIPS structure with the transformation equations
chstr = mrdfits(rootdir+'catalogs/final/v6/'+globalfield+'_combined_chips.fits.gz',1,/silent)
chstr.base = strtrim(chstr.base,2)

; Transfer the transformation info
add_tag,filestr,'colband','',filestr
add_tag,filestr,'colsign',0,filestr
add_tag,filestr,'zpterm',0.0,filestr
add_tag,filestr,'amterm',0.0,filestr
add_tag,filestr,'colterm',0.0,filestr
MATCH,chstr.expnum+'_'+string(chstr.chip,format='(i02)'),filestr.expnum+'_'+string(filestr.chip,format='(i02)'),ind1,ind2,/sort,count=nmatch
filestr[ind2].colband = chstr[ind1].colband
filestr[ind2].colsign = chstr[ind1].colsign
filestr[ind2].zpterm = chstr[ind1].zpterm
filestr[ind2].amterm = chstr[ind1].amterm
filestr[ind2].colterm = chstr[ind1].colterm
; it's mainly ZPTERM that will be different between the
; TRANS_FITSTR and CHSTR

;for i=0,nfiles-1 do begin
;  mjdchipfilt = strtrim(filestr[i].mjd,2)+'-'+strtrim(filestr[i].chip,2)+'-'+strtrim(filestr[i].filter,2)  
;  MATCH,trans_mjdchipfilt,mjdchipfilt,tind1,tind2,/sort,count=ntmatch
;  filestr[i].colband = trans_fitstr[tind1[0]].colband
;  filestr[i].colsign = trans_fitstr[tind1[0]].colsign
;  filestr[i].zpterm = trans_fitstr[tind1[0]].zpterm
;  filestr[i].amterm = trans_fitstr[tind1[0]].amterm
;  filestr[i].colterm = trans_fitstr[tind1[0]].colterm
;endfor

; For BAD exposures, just use the standard transformation equations
bd = where(filestr.colband eq '',nbd)
for i=0,nbd-1 do begin
  mjdchipfilt = strtrim(filestr[bd[i]].mjd,2)+'-'+strtrim(filestr[bd[i]].chip,2)+'-'+strtrim(filestr[bd[i]].filter,2)  
  MATCH,trans_mjdchipfilt,mjdchipfilt,tind1,tind2,/sort,count=ntmatch
  filestr[bd[i]].colband = trans_fitstr[tind1[0]].colband
  filestr[bd[i]].colsign = trans_fitstr[tind1[0]].colsign
  filestr[bd[i]].zpterm = trans_fitstr[tind1[0]].zpterm
  filestr[bd[i]].amterm = trans_fitstr[tind1[0]].amterm
  filestr[bd[i]].colterm = trans_fitstr[tind1[0]].colterm
endfor


; Load the aperture correction file
apcorfile = dir+'apcor.lst.orig'
apcor = IMPORTASCII(apcorfile,fieldnames=['name','value'],/noprint)
; Remove the 'a.del' endings for the names
apcor_orig = apcor
apcor.name = repstr(apcor.name,'a.del','')  ; base names 
add_tag,filestr,'apcor',0.0,filestr
MATCH,apcor.name,file_basename(filestr.file,'.fits'),aind1,aind2,/sort,count=namatch
filestr[aind2].apcor = apcor[aind1].value

; Get spatial coverage of the data
rar = minmax(filestr.vertices_ra)
decr = minmax(filestr.vertices_dec)
cenra = mean(rar)
cendec = mean(decr)
print,'CENRA = ',stringize(cenra,ndec=4)
print,'CENDEC = ',stringize(cendec,ndec=4)

; Use tangent pane projection
vra = (filestr.vertices_ra)(*)
vdec = (filestr.vertices_dec)(*)
ROTSPHCEN,vra,vdec,cenra,cendec,vlon,vlat,/gnomic
lonr = minmax(vlon)
latr = minmax(vlat)

; Create RA/DEC positions for the stars
sep = 7.0 ; 10.0  ; AST separation in arcsec
nlon = floor(range(lonr)/sep*3600.0)-1
nlat = floor(range(latr)/sep*3600.0)-1
ntot = nlon*nlat
print,'NLON = ',strtrim(nlon,2)
print,'NLAT = ',strtrim(nlat,2)
print,'NTOT = ',strtrim(ntot,2)

; Max/min g-band magnitude
gmin = 17.0
gmax = 27.0 ;26.2 ;25.0
gdsynth = where(synth.g ge gmin and synth.g le gmax,ngdsynth)
synth1 = synth[gdsynth]

; Randomly select the number of artificial stars that we need
;  with replacement
ind = floor(randomu(seed,ntot)*(n_elements(synth1)-1))
fsynth = synth1[ind]
lon1 = findgen(nlon)*(sep/3600.0d0) + 0.5*sep/3600.0 + lonr[0]
lon2 = lon1#replicate(1,nlat)
flon = (lon2)(*)
lat1 = findgen(nlat)*(sep/3600.0d0) + 0.5*sep/3600.0 + latr[0]
lat2 = replicate(1,nlon)#lat1
flat = (lat2)(*)
fsynth.lon = flon
fsynth.lat = flat
ROTSPHCEN,flon,flat,cenra,cendec,rr,dd,/reverse,/gnomic
fsynth.ra = rr
fsynth.dec = dd
; Make new ID
fsynth.id = lindgen(n_elements(fsynth))+1

; Image loop
stags = tag_names(fsynth)
undefine,fakestr
undefine,allsynth
For i=0,nfiles-1 do begin

  head = headfits(filestr[i].file)

  ; Find artificial stars within the boundary of this image
  rar1 = minmax(filestr[i].vertices_ra)
  decr1 = minmax(filestr[i].vertices_dec)
  off = 0.02
  ind1 = where(fsynth.ra ge rar1[0]-off and fsynth.ra le rar1[1]+off and $
               fsynth.dec ge decr1[0]-off and fsynth.dec le decr1[1]+off,nind1)
  fsynth1 = fsynth[ind1]

  ; Use WCS to convert RA/DEC to X/Y and do final spatial cut
  head_adxy,head,fsynth1.ra,fsynth1.dec,x,y,/deg
  ind2 = where(x ge 0 and x le filestr[i].nx-1 and y ge 0 and y le filestr[i].ny-1,nind2)
  fsynth2 = fsynth1[ind2]
  x2 = x[ind2]
  y2 = y[ind2]

  ; Use photometric transformation equations to get instrumental
  ; magnitudes in this image (filter, exptime).
  
  ;; calmag = instmag - t.zpterm - t.amterm*t.airmass - t.colterm*clr - t.amcolterm*t.airmass*clr - 
  ;;           t.colsqterm*clr*clr - t.apcor + 2.5*alog10(t.exptime)
  ;; instmag = calmag + t.zpterm + t.amterm*t.airmass + t.colterm*clr
  ;;             + t.apcor - 2.5*alog10(t.exptime)
  magind = where(stags eq strupcase(filestr[i].filter),nmagind)
  colind = where(stags eq strupcase(filestr[i].colband),ncolind)
  mag = fsynth2.(magind[0])
  ; Make the color
  if filestr[i].colsign eq 1 then begin
    clr = mag - fsynth2.(colind[0])
    ;clr = tempmag - colmag
  endif else begin
    clr = fsynth2.(colind[0]) - mag
    ;clr = colmag - tempmag
  endelse
  t = filestr[i]
  instmag = mag + t.zpterm + t.amterm*t.airmass + t.colterm*clr + $
              t.apcor - 2.5*alog10(t.exptime)

  ; Construct fake ALLSTAR catalog for addstar
  addcat = replicate({id:0L,x:0.0,y:0.0,mag:0.0},nind2)
  addcat.id = fsynth2.id
  addcat.x = x2+1  ; IRAF-like coordinates
  addcat.y = y2+1  ; IRAF-like coordinates
  addcat.mag = instmag

  ; Use daophot_addstar. pro to add artificial stars
  infile = filestr[i].dir+'/'+filestr[i].base
  fname = filestr[i].field+fakefield
  ;fname = filestr[i].field+'T1'
  outfile = filestr[i].dir+'/'+fname+'-'+filestr[i].expnum+'_'+string(filestr[i].chip,format='(i02)')+'.fits'
  print,'Adding ',strtrim(nind2,2),' artificial stars to ',file_basename(outfile)
  DAOPHOT_ADDSTAR,infile,addcat,outfile,/clobber

  ; Save the input catalog
  add_tag,fsynth2,'x',0.0,fsynth2  
  add_tag,fsynth2,'y',0.0,fsynth2  
  add_tag,fsynth2,'mag',0.0,fsynth2  
  fsynth2.x = addcat.x
  fsynth2.y = addcat.y
  fsynth2.mag = addcat.mag
  fakecat = filestr[i].dir+'/'+fname+'-'+filestr[i].expnum+'_'+string(filestr[i].chip,format='(i02)')+'_input.fits'
  MWRFITS,fsynth2,fakecat,/create

  ; Make structure of all artificial stars
  ;  deal with duplicates below
  PUSH,allsynth,fsynth2

  ; Create symoblic links to the PSF and opt files
  origbase = filestr[i].dir+'/'+filestr[i].base
  newbase = filestr[i].dir+'/'+fname+'-'+filestr[i].expnum+'_'+string(filestr[i].chip,format='(i02)')
  FILE_DELETE,newbase+['.psf','.opt','.als.opt'],/allow
  FILE_LINK,origbase+['.psf','.opt','.als.opt'],newbase+['.psf','.opt','.als.opt']

  ; Make FAKESTR structure
  fakestr1 = filestr[i]
  add_tag,fakestr1,'fake_ninput',0L,fakestr1
  fakestr1.fake_ninput = nind2
  add_tag,fakestr1,'fake_fits','',fakestr1
  fakestr1.fake_fits = filestr[i].dir+'/'+fname+'-'+filestr[i].expnum+'_'+string(filestr[i].chip,format='(i02)')+'.fits'
  add_tag,fakestr1,'fake_cat','',fakestr1
  fakestr1.fake_cat = filestr[i].dir+'/'+fname+'-'+filestr[i].expnum+'_'+string(filestr[i].chip,format='(i02)')+'_input.fits'
  PUSH,fakestr,fakestr1

  ;stop

Endfor 

; Get unique list of artificial stars input
ui = uniq(allsynth.id,sort(allsynth.id))
allsynth = allsynth[ui]
outsynthfile = dir+field+'/'+field+fakefield+'-input.fits'
;outsynthfile = dir+'F2/F2T1-input.fits'
print,'Saving artificial star list to ',outsynthfile
MWRFITS,allsynth,outsynthfile,/create

; Create a new "apcor.lst" list
names = file_basename(fakestr.fake_fits,'.fits')+'a.del'
maxlen = max(strlen(names))
fmt = '(A-'+strtrim(maxlen+3,2)+',F12.8)'
writecol,dir+'apcor.lst',file_basename(fakestr.fake_fits,'.fits')+'a.del',fakestr.apcor,fmt=fmt

; Save the FAKESTR structure
outfakestr = dir+field+'/'+field+'-fakestr.fits'
;outfakestr = dir+'F2/F2-fakestr.fits'
print,'Saving FAKESTR to ',outfakestr
MWRFITS,fakestr,outfakestr,/create

; Create new MCH files and symlinks for the combined files
mchfiles = file_search(dir+field+'/'+field+'-*_??.mch',count=nmchfiles)
;mchfiles = file_search(dir+'F2/F2-*_??.mch',count=nmchfiles)
for i=0,nmchfiles-1 do begin
  mchfile = mchfiles[i]
  mchdir = file_dirname(mchfile)
  mchbase = file_basename(mchfile,'.mch')
  newmchfile = mchdir+'/'+repstr(mchbase,field+'-',field+fakefield+'-')+'.mch'
  ;newmchfile = mchdir+'/'+repstr(mchbase,'F2-','F2T1-')+'.mch'
  READLINE,mchfile,mchlines
  mchlines2 = repstr(mchlines,"'"+field+"-","'"+field+fakefield+"-")
  ;mchlines2 = repstr(mchlines,"'F2-","'F2T1-")
  WRITELINE,newmchfile,mchlines2

  ; Need weights, .scale, .zero, _comb.psf, _shift.mch, .mag
  origbase = mchdir+'/'+mchbase
  newbase = mchdir+'/'+repstr(mchbase,field+'-',field+fakefield+'-')
  ;newbase = mchdir+'/'+repstr(mchbase,'F2-','F2T1-')
  FILE_DELETE,newbase+['.weights','.scale','.zero','_comb.psf','_comb.opt','_comb.als.opt','_shift.mch','.mag'],/allow
  FILE_LINK,origbase+['.weights','.scale','.zero','_comb.psf','_comb.opt','_comb.als.opt','_shift.mch','.mag'],$
            newbase+['.weights','.scale','.zero','_comb.psf','_comb.opt','_comb.als.opt','_shift.mch','.mag'],/allow

  ; Copy _comb.mch file (if it exists) and rename the files inside of it
  if file_test(origbase+'_comb.mch') eq 1 then begin
    READLINE,origbase+'_comb.mch',combmchlines
    combmchlines2 = repstr(combmchlines,field+'-',field+fakefield+'-')
    WRITELINE,newbase+'_comb.mch',combmchlines2
  endif

endfor

;stop

end
