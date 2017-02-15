pro getcomplete,field,dir=dir

; Figure out completeness from artificial stars

; Put together list of original detected source + input artificial
; stars. Then match those up with the final list of real sources
; and artificial stars.

resolve_routine,'photcalib',/compile_full_file

fakedir = '/datalab/users/dnidever/smash/cp/red/photred/addfakes/'

; Not enough inputs
if n_elements(field) eq 0 then begin
  print,'Syntax - getcomplete,field,dir=dir'
  return
endif

if n_elements(dir) eq 0 then begin
  cd,current=dir
  dir = dir+'/'
endif

;globalfield = 'Field100'
;field = 'F2'
fakefield = 'T1'

; Load "fields" file to get "global" field name
READCOL,dir+'fields',shnames,lnames,format='A,A',/silent
globalfield = lnames[0]

print,'Getting completeness function for ',globalfield,' in ',dir+field
print,'------------------------------------------------------------'


;origfile = 'F2-00423049_01.phot'
;orig = IMPORTASCII(dir+'F2/'+origfile,/header)
origfile = file_search(dir+field+'/'+field+'-*_01.phot',count=norigfile)
orig = IMPORTASCII(origfile,/header,/silent)
norig = n_elements(orig)
print,'NORIG = ',strtrim(norig,2)
synthfile = field+fakefield+'-input.fits'
;synthfile = 'F2T1-input.fits'
synth = MRDFITS(dir+field+'/'+synthfile,1,/silent)
;synth = MRDFITS(dir+'F2/'+synthfile,1,/silent)
nsynth = n_elements(synth)
print,'NSYNTH = ',strtrim(nsynth,2)
;finalfile = 'F2T1-00423049_01.ast'
;final = IMPORTASCII(dir+'F2/'+finalfile,/header)
finalfile = file_search(dir+field+'/'+field+fakefield+'-*_01.ast',count=norigfile)
final = IMPORTASCII(finalfile,/header,/silent)
nfinal = n_elements(final)
print,'NFINAL = ',strtrim(nfinal,2)

fakestr = MRDFITS(dir+field+'/'+field+'-fakestr.fits',1,/silent)
;fakestr = MRDFITS(dir+'F2/F2-fakestr.fits',1,/silent)
mchfile = file_search(dir+field+'/'+field+fakefield+'-*_01.mch',count=nmchfile)
;mchfile = dir+'F2/F2T1-00423049_01.mch'
loadmch,mchfile,mchlines
; put the FAKESTR elements in the right order
MATCH,file_basename(mchlines,'.als'),file_basename(fakestr.fake_fits,'.fits'),ind1,ind2,/sort
fakestr = fakestr[ind2]

; Load the final transformation equations for this field
reduxdir = smashred_rootdir()+'cp/red/photred/'
chstr = mrdfits(reduxdir+'catalogs/final/v5/'+globalfield+'_combined_chips.fits.gz',1,/silent)

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
; MOST OF THESE ARE ALREADY IN FAKESTR FROM WHEN THE CALIBRATED
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

obj = replicate({id:'',photid:'',ra:0.0d0,dec:0.0d0,ndet:0L,u:0.0,uerr:0.0,ndetu:0L,$
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
  if nind gt 1 then begin
    AVERAGEMAG,outmag[*,ind],outerr[*,ind],avgmag,avgerr,/robust
    ndet = total(outmag[*,ind] lt 50,2)
  endif else begin
    avgmag = reform(outmag[*,ind])
    avgerr = reform(outerr[*,ind])
    ndet = long(avgmag lt 50)
  endelse
  magind = where(otags eq strupcase(filters[i]),nmagind)
  errind = where(otags eq strupcase(filters[i])+'ERR',nerrind)
  detind = where(otags eq 'NDET'+strupcase(filters[i]),ndetind)
  obj.(magind) = avgmag
  obj.(errind) = avgerr
  obj.(detind) = ndet
  obj.ndet += obj.(detind)
endfor

; if two stars land right on top of each other and only one is recovered then use the magnitude
; to match it up to the correct star (real or artificial)

; match obj to orig and obj to synth
SRCMATCH,obj.ra,obj.dec,orig.ra,orig.dec,0.5,oind1,oind2,/sph,count=nomatch
SRCMATCH,obj.ra,obj.dec,synth.ra,synth.dec,0.5,aind1,aind2,/sph,count=nastmatch
; then deal with "duplicate" matches
; put everything in one structure with a FAKE column and columns for
; the recovered values

; Deal with duplicates
;  all we care about is real objects falsely
;  identified as ASTs, so bad matches in AIND
objdblind = doubles([oind1,aind1],count=nobjdbl)
if nobjdbl gt 0 then begin
  print,'Resolving ',strtrim(nobjdbl,2),' duplicates'
  objdbl = ([oind1,aind1])(objdblind)
  flag = lonarr(nobjdbl)  ; 0-real, 1-ast
  origdbl = lonarr(nobjdbl)
  synthdbl = lonarr(nobjdbl)
  for i=0,nobjdbl-1 do begin
    objdbl1 = objdbl[i]  ; the OBJ index
    obj1 = obj[objdbl1]
    objmag = [obj1.u, obj1.g, obj1.r, obj1.i, obj1.z]

    ; ORIG
    MATCH,oind1,objdbl1,ind1,/sort
    origdbl[i] = ind1  ; index into OIND1/2
    orig1 = orig[oind2[ind1]]
    odist = sphdist(obj1.ra,obj1.dec,orig1.ra,orig1.dec,/deg)*3600.
    origmag = [orig1.umag, orig1.gmag, orig1.rmag, orig1.imag, orig1.zmag]
    gdorig = where(objmag lt 50 and origmag lt 50,ngdorig)
    if ngdorig gt 0 then omagdiff = mean(abs(origmag[gdorig]-objmag[gdorig])) else omagdiff=99.99
    ofinaldiff = sqrt(odist^2 + omagdiff^2)

    ; SYNTH
    MATCH,aind1,objdbl1,ind2,/sort
    synthdbl[i] = ind2  ; index into AIND1/2
    synth1 = synth[aind2[ind2]]
    adist = sphdist(obj1.ra,obj1.dec,synth1.ra,synth1.dec,/deg)*3600.
    synthmag = [synth1.u, synth1.g, synth1.r, synth1.i, synth1.z]
    gdsynth = where(objmag lt 50 and synthmag lt 50,ngdsynth)
    if ngdsynth gt 0 then amagdiff = mean(abs(synthmag[gdsynth]-objmag[gdsynth])) else amagdiff=99.99
    afinaldiff = sqrt(adist^2 + amagdiff^2)

    ; Which one is the match
    if ofinaldiff lt afinaldiff then begin
      com='  REAL'
      flag[i] = 0
    endif else begin
      com='  AST'
      flag[i] = 1
    endelse
    if keyword_set(verbose) then $
      print,strtrim(i+1,2),odist,omagdiff,ofinaldiff,' ',adist,amagdiff,afinaldiff,com
    ;if amagdiff lt omagdiff and adist gt odist then stop
  endfor

  ; Remove ASTs from the "ORIG" list
  bdomatch = where(flag eq 1,nbdomatch,ncomp=ngdomatch)  
  print,strtrim(nbdomatch,2),' are ASTs and ',strtrim(ngdomatch,2),' are REAL sources'
  if nbdomatch gt 0 then begin
    if nbdomatch lt nomatch then begin
      bdorigdbl = origdbl[bdomatch]
      REMOVE,bdorigdbl,oind1,oind2
      nomatch = n_elements(oind1)
    endif else begin
      undefine,oind1,oind2
      nomatch = 0
    endelse
  endif
endif ; duplicates

; "Prune" the real sources from the list of recovered sources
left = obj
if nomatch gt 0 then remove,oind1,left
; Now rematch SYNTH to the leftover sources
SRCMATCH,left.ra,left.dec,synth.ra,synth.dec,0.5,aind1,aind2,/sph,count=nastmatch
print,strtrim(nastmatch,2),' ASTs recovered'

; Put all AST information into one structure
schema = {astid:'',photid:'',recovered:-1,inp_ra:0.0d0,inp_dec:0.0d0,$
          inp_u:0.0,inp_g:0.0,inp_r:0.0,inp_i:0.0,inp_z:0.0,$
          id:'',ra:0.0d0,dec:0.0d0,ndet:0L,u:0.0,uerr:0.0,ndetu:0L,$
          g:0.0,gerr:0.0,ndetg:0L,r:0.0,rerr:0.0,ndetr:0L,$
          i:0.0,ierr:0.0,ndeti:0L,z:0.0,zerr:0.0,ndetz:0L,$
          chi:0.0,sharp:0.0,flag:0L,prob:0.0}
ast = replicate(schema,nsynth)
; Copy over the input synth data
ast.astid = synth.id
ast.photid = synth.photid
ast.inp_ra = synth.ra
ast.inp_dec = synth.dec
ast.inp_u = synth.u
ast.inp_g = synth.g
ast.inp_r = synth.r
ast.inp_i = synth.i
ast.inp_z = synth.z
; Copy over the recovered values
temp = ast[aind2]
STRUCT_ASSIGN,left[aind1],temp,/nozero
ast[aind2] = temp
;ast[aind1].ast = 1
ast[aind2].recovered = 1

; Write out the final file
outastfile = dir+globalfield+'_complete.fits'
print,'Writing final catalog to ',outastfile
MWRFITS,ast,outastfile,/create
if file_test(outastfile+'.gz') then file_delete,outastfile+'.gz'
spawn,['gzip',outastfile],/noshell

; Figure out the completeness
;gdall = where(comb.ast eq 1,ngdall)
;gdrecover = where(comb.ast eq 1 and comb.recovered eq 1,ngdrecover)
gdrecover = where(ast.recovered eq 1,ngdrecover)
dx = 0.2
dy = 0.4
xr = [-1,3.5]
yr = [17.0,27.0]
;hess,ast[gdall].inp_g-ast[gdall].inp_i,ast[gdall].inp_g,dum,imall,dx=dx,dy=dy,xr=xr,yr=yr,xarr=xarr,yarr=yarr,/noplot
hess,ast.inp_g-ast.inp_i,ast.inp_g,dum,imall,dx=dx,dy=dy,xr=xr,yr=yr,xarr=xarr,yarr=yarr,/noplot
hess,ast[gdrecover].inp_g-ast[gdrecover].inp_i,ast[gdrecover].inp_g,dum,imrec,dx=dx,dy=dy,xr=xr,yr=yr,xarr=xarr,yarr=yarr,/noplot
;displayc,float(imrec)/(imall>1),xarr,yarr,/yflip,xtit='g-i',ytit='g',tit='Completeness'

; Make some figures
if file_test(dir+'plots/',/directory) eq 0 then file_mkdir,dir+'plots/'
setdisp
;loadcol,3
!p.font = 0
; Input ASTS
file = dir+'plots/'+globalfield+'_input'
ps_open,file,/color,thick=4,/encap
device,/inches,xsize=8.5,ysize=10.5
displayc,imall,xarr,yarr,/yflip,xtit='g-i',ytit='g',tit='Input ASTs for '+globalfield,charsize=1.3
ps_close
ps2png,file+'.eps',/eps
spawn,['epstopdf',file+'.eps'],/noshell
; Recovered
file = dir+'plots/'+globalfield+'_recovered'
ps_open,file,/color,thick=4,/encap
device,/inches,xsize=8.5,ysize=10.5
displayc,imrec,xarr,yarr,/yflip,xtit='g-i',ytit='g',tit='Recovered ASTs for '+globalfield,charsize=1.3
ps_close
ps2png,file+'.eps',/eps
spawn,['epstopdf',file+'.eps'],/noshell
; Completeness
file = dir+'plots/'+globalfield+'_completeness'
ps_open,file,/color,thick=4,/encap
device,/inches,xsize=8.5,ysize=10.5
displayc,float(imrec)/(imall>1),xarr,yarr,/yflip,xtit='g-i',ytit='g',tit='Completeness for '+globalfield,charsize=1.3
ps_close
ps2png,file+'.eps',/eps
spawn,['epstopdf',file+'.eps'],/noshell
; Combine the figures
pdffiles = dir+'/plots/'+globalfield+'_'+['input','recovered','completeness']+'.pdf'
spawn,'gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile='+dir+'/plots/'+globalfield+'_complete.pdf '+strjoin(pdffiles,' ')

;stop

end
