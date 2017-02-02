pro getcomplete

; Figure out completeness from artificial stars

; Put together list of original detected source + input artificial
; stars. Then match those up with the final list of real sources
; and artificial stars.

resolve_routine,'photcalib',/compile_full_file

dir = '/datalab/users/dnidever/smash/cp/red/photred/addfakes/Field100/'
field = 'Field100'

origfile = 'F2-00423049_01.phot'
orig = IMPORTASCII(dir+'F2/'+origfile,/header)
synthfile = 'F2T1-input.fits'
synth = MRDFITS(dir+'F2/'+synthfile,1,/silent)
finalfile = 'F2T1-00423049_01.ast'
final = IMPORTASCII(dir+'F2/'+finalfile,/header)


fakestr = MRDFITS(dir+'F2/F2-fakestr.fits',1,/silent)
mchfile = dir+'F2/F2T1-00423049_01.mch'
loadmch,mchfile,mchlines
; put the FAKESTR elements in the right order
MATCH,file_basename(mchlines,'.als'),file_basename(fakestr.fake_fits,'.fits'),ind1,ind2,/sort
fakestr = fakestr[ind2]

; Load the final transformation equations for this field
reduxdir = smashred_rootdir()+'cp/red/photred/'
chstr = mrdfits(reduxdir+'catalogs/final/v4/'+field+'_combined_chips.fits.gz',1,/silent)

; Get the information we need for our 
add_tag,fakestr,'photometric',0B,fakestr
add_tag,fakestr,'band','',fakestr
fakestr.band = fakestr.filter
;add_tag,fakestr,'colband','',fakestr
;add_tag,fakestr,'colsign',0,fakestr
;add_tag,fakestr,'zpterm',0.0d0,fakestr
add_tag,fakestr,'zptermsig',0.0d0,fakestr
;add_tag,fakestr,'amterm',0.0d0,fakestr
add_tag,fakestr,'amtermsig',0.0d0,fakestr
;add_tag,fakestr,'colterm',0.0d0,fakestr
add_tag,fakestr,'coltermsig',0.0d0,fakestr
add_tag,fakestr,'amcolterm',0.0d0,fakestr
add_tag,fakestr,'amcoltermsig',0.0d0,fakestr
add_tag,fakestr,'colsqterm',0.0d0,fakestr
add_tag,fakestr,'colsqtermsig',0.0d0,fakestr
MATCH,fakestr.expnum+'_'+string(fakestr.chip,format='(i02)'),chstr.expnum+'_'+string(chstr.chip,format='(i02)'),ind1,ind2,/sort
fakestr[ind1].photometric = chstr[ind2].photometric
;fakestr[ind1].band = chstr[ind2].band
;fakestr[ind1].colband = chstr[ind2].colband
;fakestr[ind1].colsign = chstr[ind2].colsign
;fakestr[ind1].zpterm = chstr[ind2].zpterm
;fakestr[ind1].zptermsig = chstr[ind2].zptermsig
;fakestr[ind1].amterm = chstr[ind2].amterm
;fakestr[ind1].amtermsig = chstr[ind2].amtermsig
;fakestr[ind1].colterm = chstr[ind2].colterm
;fakestr[ind1].coltermsig = chstr[ind2].coltermsig
;fakestr[ind1].amcolterm = chstr[ind2].amcolterm
;fakestr[ind1].amcoltermsig = chstr[ind2].amcoltermsig
;fakestr[ind1].colsqterm = chstr[ind2].colsqterm
;fakestr[ind1].colsqtermsig = chstr[ind2].colsqtermsig
; MOST OF THESE ALREADY IN FAKESTR FROM WHEN THE CALIBRATED
; SYNTH MAGNITUDES WERE CONVERTED TO INSTRUMENTAL ONES!!!

; Need to properly calibrate to get REAL magnitudes.
; use photcalib.pro or smashred_apply_phottranseqn.pro
; can we use the phot trans file and make one "allnights.trans" to use
; with CALIB?  but how about non-photometric data??
; use the final CHSTR with it's phot trans info to apply my own calibration


;trans1 = {night:-1L,chip:-1,band:'',color:'',colband:'',colsign:0,zpterm:0.0d,amterm:0.0d,colterm:0.0d,$
;          amcolterm:0.0d,colsqterm:0.0d,zptermsig:0.0d,amtermsig:0.0d,coltermsig:0.0d,$
;          amcoltermsig:0.0d,colsqtermsig:0.0d}
;; one element per observation
;trans = replicate(trans1,numobs)

; Make INSTAR
nfinal = n_elements(final)
ntags = n_tags(final)
instar = dblarr(nfinal,ntags)
for i=0,ntags-1 do instar[*,i]=final.(i)

; Make TRANS
trans = fakestr

; Make INP
numobs = n_elements(fakestr)
inp = {magfile:'',outfile:'',night:lonarr(numobs),chip:lonarr(numobs),band:strarr(numobs),$
       airmass:dblarr(numobs),exptime:dblarr(numobs),apcorr:dblarr(numobs)}
inp.magfile = finalfile
inp.outfile = repstr(finalfile,'.ast','.phot')
inp.night = fakestr.mjd
inp.chip = fakestr.chip
inp.band = fakestr.band
inp.airmass = fakestr.airmass
inp.exptime = fakestr.exptime
inp.apcorr = fakestr.apcor

; Use SOLVESTAR from photcalib.pro
print,'Running SOLVESTAR'
SOLVESTAR,instar,trans,inp,outstar
; instar  id, x, y, unsolved magnitudes, chi, sharp  
;  inmag = instar[*,2*lindgen(numobs)+3]
;  inerr = instar[*,2*lindgen(numobs)+4]
; trans   transformation structure, one for each observation
;           night, chip, band, colband, colsign, zpterm, zptermsig, amterm,
;           amtermsig, colterm, coltermsig, amcolterm, amcoltermsig, colsqterm, colsqtermsig
; inp     the "input" structure with airmass, exptime, band, apcorr,
;           night, chip. one per magfile, one element per observation.
;         one element per photometry file, but night/chip/band/airmass/exptime/apcorr
;           are arrays with one element per observation.
; outstar   final output array
;  outstar[*,2*lindgen(numobs)+3] = tempmag
;  outstar[*,2*lindgen(numobs)+4] = temperr
outmag = outstar[*,2*lindgen(numobs)+3]
outerr = outstar[*,2*lindgen(numobs)+4]

obj = replicate({id:'',photid:'',ra:0.0d0,dec:0.0d0,u:0.0,uerr:0.0,ndetu:0L,$
                 g:0.0,gerr:0.0,ndetg:0L,r:0.0,rerr:0.0,ndetr:0L,$
                 i:0.0,ierr:0.0,ndeti:0L,z:0.0,zerr:0.0,ndetz:0L,$
                 chi:0.0,sharp:0.0,flag:0L,prob:0.0},nfinal)
obj.id = strtrim(final.id,2)
obj.ra = final.ra
obj.dec = final.dec
obj.chi = final.chi
obj.sharp = final.sharp
obj.flag = final.flag
obj.prob = final.prob
otags = tag_names(obj)

; Now we need to get average magnitudes for each band
filters = ['u','g','r','i','z']
nfilters = n_elements(filters)
for i=0,nfilters-1 do begin
  ind = where(trans.band eq filters[i],nind)
  AVERAGEMAG,outmag[*,ind],outerr[*,ind],avgmag,avgerr,/robust
  magind = where(otags eq strupcase(filters[i]),nmagind)
  errind = where(otags eq strupcase(filters[i])+'ERR',nerrind)
  detind = where(otags eq 'NDET'+strupcase(filters[i]),ndetind)
  obj.(magind) = avgmag
  obj.(errind) = avgerr
  obj.(detind) = total(outmag[*,ind] lt 50,2)
endfor

; if two stars land right on top of each other and only one is recovered then use the magnitude
; to match it up to the correct star (real or artificial)

; match obj to orig and obj to synth
SRCMATCH,obj.ra,obj.dec,orig.ra,orig.dec,0.5,ind1a,ind2a,/sph,count=nmatch1
SRCMATCH,obj.ra,obj.dec,synth.ra,synth.dec,0.5,ind1b,ind2b,/sph,count=nmatch2
; then deal with "duplicate" matches
; put everything in one structure with a FAKE column and columns for
; the recovered values

obj2 = obj[ind1b]
synth2 = synth[ind2b]

plot,synth2.g,obj2.g-synth2.g,ps=3,xr=[14,27],yr=[-3,3]
; Looks good

stop


end
