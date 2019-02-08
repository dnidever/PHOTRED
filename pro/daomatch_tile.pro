;+
;
; DAOMATCH_TILE
;
; This is very similar to the DAOMATCH.PRO program that
; matches stars and finds the transformation but it does
; it using the "tiling" coordinate system.
;
; INPUTS:
;  files     Array of ALS files,  The first file will be used
;            as the reference.  If /fake set then the first ALS file
;            in "files" should already have an associated MCH file.
;  tilestr   The tiling scheme structure.
;  =maxshift Constraints on the initial X/Y shifts.
;  /usewcs   Use the WCS for initial matching.  This is the default.
;  /fake     Run for artificial stars.  The MCH file should be input
;              and daomaster run to create raw/tfr files.
;  /verbose  Verbose output
;  /stp      Stop at the end of the program
;  =hi       Not used anymore.
;
; OUTPUTS:
;  An MCH, TFR and RAW file with basename of the first file.
;  =error    The error message if one occurred.
;
; USAGE:
;  IDL>daomatch_tile,['obj1034_1.als','obj1035_1.als','obj1036_1.als']
;
; Add options to freeze the scale and/or rotation.
; 
; By D. Nidever   December 2006
;-

pro daomatch_dummy
FORWARD_FUNCTION test_trans, trans_coo, trans_coo_dev
end

;---------------------------------------------------------------

function test_trans,trans

; This function tests if a transformation equation from
; daomatch is good or not.  The scale should be nearly 1.0
; and the rotation should be near 0.
;
; Return value:
;  1  Good
;  0  Bad
; -1  Fail
;

test = -1

if n_elements(trans) eq 0 then return,-1

; The test mainly looks at the rotation/scale values
; and not the xoff/yoff values.

sz = size(trans)
if sz[0] ne 2 or sz[1] ne 2 or sz[2] ne 6 then return,-1


xoff = (reform(trans[1,0]))(0)
yoff = (reform(trans[1,1]))(0)
c = (reform(trans[1,2]))(0)
e = (reform(trans[1,3]))(0)
d = (reform(trans[1,4]))(0) 
f = (reform(trans[1,5]))(0)
; from ccdpck.txt
;              x(1) = A + C*x(n) + E*y(n)
;              y(1) = B + D*x(n) + F*y(n)

; C=F~1 and D=E~0
test = 1
if abs(c-f) gt 0.1 then test=0
if abs(d-e) gt 0.1 then test=0
if abs(c-1.0) gt 0.1 then test=0
if abs(e) gt 0.1 then test=0

return,test

end

;---------------------------------------------------------------

pro daomatch_tile,files,tilestr,groupstr,mchbase,stp=stp,verbose=verbose,logfile=logfile,$
             error=error,maxshift=maxshift,fake=fake

t0 = systime(1)

undefine,error

nfiles = n_elements(files)
if nfiles eq 0 then begin
  print,'Syntax - daomatch_tile,files,tilestr,groupstr,mchbase,fake=fake,stp=stp,verbose=verbose'
  error = 'Not enough inputs'
  return
end

; Logfile
if keyword_set(logfile) then logf=logfile else logf=-1

; Compile MATCHSTARS.PRO
RESOLVE_ROUTINE,'matchstars',/compile_full_file

; Current directory
CD,current=curdir

dir = FILE_DIRNAME(files[0])
CD,dir

bases = FILE_BASENAME(files,'.als')

; FAKE, running for artificial star tests
if keyword_set(fake) then begin
  ; Check that MCH file exists
  if file_test(mchbase+'.mch') eq 0 then begin
    error = '/fake set but '+mchbase+'.mch NOT FOUND'
    printlog,logf,error
    return
  endif
  ;; Skip the MCH creation process
  goto,rundaomaster
endif


; Gather information on all of the files
printlog,logf,'Gathering file information'
fitsfiles = bases+'.fits'
bdfits = where(file_test(fitsfiles) eq 0,nbdfits)
if nbdfits gt 0 then fitsfiles[bdfits]+='.fz'
PHOTRED_GATHERFILEINFO,fitsfiles,filestr
ntrans = 6
add_tag,filestr,'catfile','',filestr
add_tag,filestr,'resamptrans',dblarr(ntrans),filestr
add_tag,filestr,'resamptransrms',0.0,filestr
filestr.catfile = bases+'.als'


; Creating MCH file in the tile coordinate system
for i=0,nfiles-1 do begin
  ; Get the header
  if strmid(filestr[i].file,6,7,/reverse_offset) eq 'fits.fz' then begin
    fhead = PHOTRED_READFILE(filestr[i].file,exten=1,/header)
    ; Fix the NAXIS1/2 in the header
    sxaddpar,fhead,'NAXIS1',sxpar(fhead,'ZNAXIS1')
    sxaddpar,fhead,'NAXIS2',sxpar(fhead,'ZNAXIS2')
  endif else begin
    fhead = PHOTRED_READFILE((filestr[i].file,/header)
  endelse

  ; Convert X/Y of this system into the combined reference frame
  ;  The pixel values are 1-indexed like DAOPHOT uses.
  ;  Use a 2D grid of points in the image and the WCS to get
  ;  the transformation.
  ngridbin = 50
  nxgrid = filestr[i].nx / ngridbin
  nygrid = filestr[i].ny / ngridbin
  xgrid = (lindgen(nxgrid)*ngridbin+1)#replicate(1,nygrid)
  ygrid = replicate(1,nxgrid)#(lindgen(nygrid)*ngridbin+1)
  HEAD_XYAD,fhead,xgrid-1,ygrid-1,ragrid,decgrid,/deg
  HEAD_ADXY,tilestr.head,ragrid,decgrid,refxgrid,refygrid,/deg
  refxgrid += 1  ; convert 0-indexed to 1-indexed
  refygrid += 1
  ; Convert to tile X/Y values
  refxgrid -= groupstr.x0
  refygrid -= groupstr.y0

  ; Now fit the transformation
  xdiff = refxgrid-xgrid
  ydiff = refygrid-ygrid
  xmed = median([xdiff],/even)
  ymed = median([ydiff],/even)
  ; Fit rotation with linear fits if enough points
  coef1 = robust_poly_fitq(ygrid,xdiff,1)  ; fit rotation term
  coef1b = dln_poly_fit(ygrid,xdiff,1,measure_errors=xdiff*0+0.1,sigma=coef1err,/bootstrap)
  coef2 = robust_poly_fitq(xgrid,ydiff,1)  ; fit rotation term
  coef2b = dln_poly_fit(xgrid,ydiff,1,measure_errors=ydiff*0+0.1,sigma=coef2err,/bootstrap)
  ;theta = mean([-coef1[1],coef2[1]])
  WMEANERR,[-coef1[1],coef2[1]],[coef1err[1],coef2err[1]],theta,thetaerr

  ; [xoff, yoff, cos(th), sin(th), -sin(th), cos(th)]
  ;trans = [xmed, ymed, 1.0, 0.0, 0.0, 1.0]
  trans = [xmed, ymed, 1.0-theta^2, theta, -theta, 1.0-theta^2]
  ; Adjust Xoff, Yoff with this transformation
  xyout = trans_coo(xgrid,ygrid,trans)
  trans[0] += median([refxgrid - xyout[0,*]],/even)
  trans[1] += median([refygrid - xyout[1,*]],/even)

  ; Fit full six parameters if there are enough stars
  fa = {x1:(refxgrid)(*),y1:(refygrid)(*),x2:(xgrid)(*),y2:(ygrid)(*)}
  initpar = trans
  fpar = MPFIT('trans_coo_dev',initpar,functargs=fa, perror=perror, niter=iter, status=status,$
               bestnorm=chisq, dof=dof, autoderivative=1, /quiet)
  trans = fpar

  diff = trans_coo_dev(fpar,x1=refxgrid,y1=refygrid,x2=xgrid,y2=ygrid)
  rms = sqrt(mean(diff^2.))
  filestr[i].resamptrans = trans
  filestr[i].resamptransrms = rms

  ; The output is:
  ; filename, xshift, yshift, 4 trans, mag offset, magoff sigma
  format = '(A2,A-30,A1,2A10,4A12,F9.3,F8.4)'
  ; In daomaster.f the translations are 10 digits with at most 4
  ; decimal places (with a leading space), the transformation
  ; coefficients are 12 digits with at most 9 decimal places.
  ; Need a leading space to separate the numbers.
  strans = ' '+[strtrim(string(trans[0:1],format='(F30.4)'),2),$
               strtrim(string(trans[2:5],format='(F30.9)'),2)]
  newline = STRING("'",filestr[i].catfile,"'", strans, 0.0, rms, format=format)
  PUSH,mchfinal,newline

  ; Printing the transformation
  printlog,logf,format='(A-22,2A10,4A12,F9.3,F8.4)',filestr[i].catfile,strans,0.0,rms
endfor

; Write to the new MCH file
mchfile = mchbase+'.mch'
WRITELINE,mchfile,mchfinal


;#####################
; Create the RAW file
;#####################
; I can't use daomaster because it won't work for non-overlapping
; images, and it always makes the first image the reference.
rundaomaster:

; Create the Schema
schema = {id:0L,x:0.0,y:0.0}
for i=0,nfiles-1 do $
  schema=create_struct(schema,'MAG'+strtrim(i+1,2),99.99,'ERR'+strtrim(i+1,2),9.99)
schema = create_struct(schema,'chi',99.99,'sharp',99.99,'nobs',0)
rawtags = tag_names(schema)

; Number of sources in each als
nalsarr = lonarr(nfiles)
for i=0,nfiles-1 do nalsarr[i]=file_lines(bases[i]+'.als')-3

; Initialize the RAW structure
raw = REPLICATE(schema,total(nalsarr))

; Loop over the images
;---------------------
tfr = lonarr(n_elements(raw),nfiles)   ; tfr array
count = 0LL
For i=0,nfiles-1 do begin
  if keyword_set(verbose) then printlog,logf,strtrim(i+1,2),' ',bases[i]

  ; Get the header
  if strmid(filestr[i].file,6,7,/reverse_offset) eq 'fits.fz' then begin
    fhead = PHOTRED_READFILE((filestr[i].file,exten=1,/header)
    ; Fix the NAXIS1/2 in the header
    sxaddpar,fhead,'NAXIS1',sxpar(fhead,'ZNAXIS1')
    sxaddpar,fhead,'NAXIS2',sxpar(fhead,'ZNAXIS2')
  endif else begin
    fhead = PHOTRED_READFILE((filestr[i].file,/header)
  endelse

  ; Load the ALS file
  LOADALS,bases[i]+'.als',als,alshead
  nals = n_elements(als)
  alsind = lindgen(nals)+1
  if i eq 0 then rawhead=alshead

  ; Convert to tile coordinates
  HEAD_XYAD,fhead,als.x-1,als.y-1,ra,dec,/deg
  HEAD_ADXY,tilestr.head,ra,dec,xref,yref,/deg
  xref += 1  ; convert 0-indexed to 1-indexed
  yref += 1
  ; Convert to tile group X/Y values
  xref -= groupstr.x0
  yref -= groupstr.y0

  ; Get mag/err columns
  magind = where(rawtags eq 'MAG'+strtrim(i+1,2),nmagind)
  errind = where(rawtags eq 'ERR'+strtrim(i+1,2),nerrind)

  ; First image
  if i eq 0 then begin
    raw[0:nals-1].id = lindgen(nals)+1
    raw[0:nals-1].x = xref
    raw[0:nals-1].y = yref
    raw[0:nals-1].(magind) = als.mag
    raw[0:nals-1].(errind) = als.err
    raw[0:nals-1].chi = als.chi
    raw[0:nals-1].sharp = als.sharp 
    raw[0:nals-1].nobs++
    ; TFR
    tfr[0:nals-1,i] = alsind
    count += nals

  ; Second and later images, crossmatch
  endif else begin

    ; Cross-match
    SRCMATCH,raw[0:count-1].x,raw[0:count-1].y,xref,yref,2,ind1,ind2,count=nmatch
    if keyword_set(verbose) then printlog,logf,strtrim(nmatch,2),' matches'

    ; Some matches, add data to existing records for these sources
    if nmatch gt 0 then begin
      raw[ind1].(magind) = als[ind2].mag
      raw[ind1].(errind) = als[ind2].err
      raw[ind1].chi += als[ind2].chi      ; cumulative sum
      raw[ind1].sharp += als[ind2].sharp  ; cumulative sum
      raw[ind1].nobs++
      ; TFR
      tfr[ind1,i] = ind2
      ; Remove stars
      if nmatch lt nals then REMOVE,ind2,als,xref,yref,alsind else undefine,als
      nals = n_elements(als)
    endif

    ; Add leftover sources
    if nals gt 0 then begin
      if keyword_set(verbose) then printlog,logf,'Adding ',strtrim(nals,2),' leftover sources'
      raw[count:count+nals-1].id = lindgen(nals)+1+count
      raw[count:count+nals-1].x = xref
      raw[count:count+nals-1].y = yref
      raw[count:count+nals-1].(magind) = als.mag
      raw[count:count+nals-1].(errind) = als.err
      raw[count:count+nals-1].chi += als.chi
      raw[count:count+nals-1].sharp += als.sharp
      raw[count:count+nals-1].nobs++
      ; TFR
      tfr[count:count+nals-1,i] = alsind
      count += nals
    endif
  endelse
Endfor
; Trim extra elements
raw = raw[0:count-1]
tfr = tfr[0:count-1,*]
; Calculate average chi/sharp
raw.chi /= (raw.nobs>1)
raw.sharp /= (raw.nobs>1)
nraw = n_elements(raw)


; Write out the RAW file
;-----------------------
openw,unit,/get_lun,mchbase+'.raw'
; Header
printf,unit,rawhead[0]
printf,unit,rawhead[1]
printf,unit,''
; Create MAG/ERR array
magarr = fltarr(nraw,nfiles*2)
for i=0,nfiles-1 do begin
  magind = where(rawtags eq 'MAG'+strtrim(i+1,2),nmagind)
  errind = where(rawtags eq 'ERR'+strtrim(i+1,2),nerrind)
  magarr[*,i*2] = raw.(magind)
  magarr[*,i*2+1] = raw.(errind)
endfor
; Only 12 mag/err/chi/sharp columns per row
nrows = ceil((nfiles*2+2)/12.)
; Loop over raw elements
for i=0,nraw-1 do begin
  ; The floating point numbers, MAG, ERR, CHI, SHARP
  arr = [reform(magarr[i,*]),raw[i].chi,raw[i].sharp]
  narr = n_elements(arr)
  ; Loop over rows for this object
  for j=0,nrows-1 do begin
    if narr gt 12 then begin
      thisarr = arr[0:11]
      arr = arr[12:*]
      narr = n_elements(arr)
    endif else begin
      thisarr = arr
      undefine,arr
      narr = 0
    endelse     
    if j eq 0 then begin
      format = '(I7,2F9.3,'+strtrim(n_elements(thisarr),2)+'F9.4)'
      printf,unit,raw[i].id,raw[i].x,raw[i].y,thisarr,format=format
    endif else begin
      ; 27 leading spaces
      format = '(A25,'+strtrim(n_elements(thisarr),2)+'F9.4)'
      printf,unit,'',thisarr,format=format
    endelse
  endfor
endfor
close,unit
free_lun,unit


; Write out TFR file
;-------------------
openw,unit,/get_lun,mchbase+'.tfr'
for i=0,nfiles-1 do printf,unit,'',bases[i]+'.als',99.9999,9.9999,format='(A1,A-30,F9.4,F9.4)'
printf,unit,' =============================='
format = '(I7,F9.2,F9.2,'+strtrim(nfiles,2)+'I7)'
for i=0,nraw-1 do printf,unit,raw[i].id,raw[i].x,raw[i].y,reform(tfr[i,*]),format=format
close,unit
free_lun,unit

; Back to the original directory
CD,curdir

if keyword_set(stp) then stop

end
