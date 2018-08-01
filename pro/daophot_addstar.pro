;+
;
; DAOPHOT_ADDSTAR
;
; Add a fake star to an image
;
; INPUTS:
;  filebase  The base of the filename for the ".psf", ".fits" and ".opt" files,
;             INCLUDING the directory path (otherwise it is assumed that they
;             are in the current directory).
;  cat       The input catalog structure of sources to add.  It must have
;              AT LEAST these columns:
;              X, Y     x/y-coordinates of sources in IRAF/DAOPHOT format (1-indexed)
;              MAG      magnitude of star
;  outfile   The name of the output FITS file.  If it already exists
;              then the star is added to the existing image (unless /clobber set).
;  /blank    Add the artificial stars to a blank image.  The default
;              is to add them to the original image.
;  /nonoise  Don't any Poisson noise to the image, by default ADDSTAR adds Poisson noise.
;  /clobber  Delete the output file if it already exists.
;  /silent   Don't print anything to the screen.
;  /stp      Stop at the end of the program.
;
; OUTPUTS:
;  The artificial stars are added to an image.  By default, they are
;  added to the original image unless /blank is set and then they are
;  added to a blank image.
;  =error     The error message if one occurred, otherwise undefined.
;
; USAGE:
;  IDL>daophot_addstar,'F1-0087651_01',cat,'fake.fits',error=error
;
; By D.Nidever  March 2016
;-

pro daophot_addstar,filebase,cat,outfile,blank=blank,nonoise=nonoise,clobber=clobber,error=error,silent=silent,stp=stp

undefine,error

; Do we have enough inputs
if n_elements(filebase) eq 0 or n_elements(cat) eq 0 or n_elements(outfile) eq 0 then begin
  error = 'Not enough inputs'
  print,'Syntax - daophot_addstar,filebase,cat,outfile,blank=blank,nonoise=nonoise,clobber=clobber,error=error,silent=silent,stp=stp'
  return
endif

; Make the psf, opt and fits filenames, with absolute paths
dir = ( file_search(file_dirname(filebase),/fully_qualify)+'/' )[0]
psffile = dir+file_basename(filebase)+'.psf'
optfile = dir+file_basename(filebase)+'.opt'
fitsfile = dir+file_basename(filebase)+'.fits'


;------------------------
; Perform error checking
;========================

; Does the PSF file exist
if file_test(psffile) eq 0 then begin
  error = psffile+' NOT FOUND'
  if not keyword_set(silent) then print,error
  return
endif
; Does the OPT file exist
if file_test(optfile) eq 0 then begin
  error = optfile+' NOT FOUND'
  if not keyword_set(silent) then print,error
  return
endif
; Does the FITS file exist
if file_test(fitsfile) eq 0 then begin
  error = fitsfile+' NOT FOUND'
  if not keyword_set(silent) then print,error
  return
endif

; Check that the PSF file is okay
READLINE,psffile,psflines,count=npsflines
if npsflines eq 0 then begin
  error = psffile+' IS EMPTY'
  if not keyword_set(silent) then print,error
  return
endif
; Check that there are no NaNs in the PSF file
bd = where(stregex(psflines[0],' NaN ',/boolean,/fold_case) eq 1,nbd)
if nbd gt 0 then begin
  error = psffile+' has NaN in first line'
  if not keyword_set(silent) then print,error
  return
endif

; Check the FITS file and header
head = headfits(fitsfile,errmsg=errmsg)
if errmsg ne '' then begin
  error = 'Error reading '+fitsfile+' header. '+errmsg
  if not keyword_set(silent) then print,error
  return
endif
; check that there's an image
nx = sxpar(head,'NAXIS1')
ny = sxpar(head,'NAXIS2')
if nx le 0 or ny le 0 then begin
  error = fitsfile+' has no image.'
  if not keyword_set(silent) then print,error
  return
endif

; Load the opt file and get GAIN
READLINE,optfile,optlines,count=noptlines
if noptlines le 0 then begin
  error = optfile+' is EMPTY.'
  if not keyword_set(silent) then print,error
  return
endif
gdgain = where(strmid(optlines,0,4) eq 'GA =',ngdgain)
if ngdgain eq 0 then begin
  error = 'No GAIN in '+optfile
  if not keyword_set(silent) then print,error
  return
endif
gain = float(first_el(strsplit(optlines[gdgain[0]],'=',/extract),/last))  ; gain (e/ADU)

; Do we have DAOPHOT?
SPAWN,['which','daophot'],out,errout,/noshell
daophotfile = FILE_SEARCH(out,count=ndaophotfile)
if (ndaophotfile eq 0) then begin
  error = 'DAOPHOT PROGRAM NOT AVAILABLE'
  if not keyword_set(silent) then print,error
  return
endif

; Is CAT a structure?
type = size(cat,/type)
if type ne 8 then begin
  error = 'CAT must be a STRUCTURE'
  if not keyword_set(silent) then print,error
  return
endif
; Does the catalog file have the correct tags?
tags = tag_names(cat)
reqtags = ['X','Y','MAG']
for i=0,n_elements(reqtags)-1 do begin
  dum = where(tags eq reqtags[i],nreqtags)
  if nreqtags eq 0 then begin
    error = 'CAT must have '+reqtags[i]+' column'
    if not keyword_set(silent) then print,error
    return
  endif
endfor

; Does the output file already exist
if file_test(outfile) eq 1 then begin
  if not keyword_set(clobber) then begin
    error = outfile+' ALREADY EXISTS and /clobber not set.'
    if not keyword_set(silent) then print,error
    return
  endif else begin
    FILE_DELETE,outfile,/allow   ; deleting existing file
  endelse
endif


;------------------
; Now run ADDSTAR
;==================

; Go to the directory with the original files
CD,current=curdir
CD,dir

; Make sure daophot.opt file exists
if FILE_TEST('daophot.opt') eq 0 then FILE_COPY,optfile,'daophot.opt',/over

; Base for temporary files
tmpbase = (MKTEMP('psf'))[0]

; Use blank sky
if keyword_set(blank) then begin
  FITS_READ,fitsfile,im,head
  inptmpfile = tmpbase+'.blank.fits'
  FILE_DELETE,inptmpfile,/allow
  MWRFITS,im*0,inptmpfile,head,/create
  inputimage = inptmpfile

; Add to original file
endif else begin
  inputimage = fitsfile
endelse


; Create input data file for DAOPHOT
;  DAOPHOT "header", I don't think the values here actually matter for ADDSTAR
undefine,hdr
push,hdr,' NL    NX    NY  LOWBAD HIGHBAD  THRESH     AP1  PH/ADU  RNOISE    FRAD'
push,hdr,'  1  2046  4094   117.7 38652.0   13.12    0.00    3.91    1.55    6.53'
; make the "fake" ALS catalog file, copied from wcsfit_daomatch.pro
ncat = n_elements(cat)
dum = {id:0L,x:0.0,y:0.0,mag:0.0,err:0.0,sky:10.0,niter:1,chi:1.0,sharp:0.0}
catals = replicate(dum,ncat)
catals.id = lindgen(ncat)+1
catals.x = cat.x
catals.y = cat.y
catals.mag = cat.mag
catals.err = 0.0
tcatfile = tmpbase+'.add'
FILE_DELETE,tcatfile,/allow
WRITEALS,tcatfile,catals,hdr

; Make output DAOPHOT script
push,lines,'#!/bin/sh'
push,lines,'export image=${1}'
push,lines,'daophot << END_DAOPHOT >> '+file_basename(tmpbase)+'.log'
push,lines,'OPTIONS'
push,lines,'${image}.opt'
push,lines,''
push,lines,'ATTACH '+file_basename(inputimage)
push,lines,'ADDSTAR'
push,lines,'${image}.psf'
push,lines,'1'
if keyword_set(nonoise) then push,lines,'99999.' else $   ; no Poisson noise added
  push,lines,stringize(gain,ndec=5)           ; photons per ADU
push,lines,file_basename(tcatfile)            ; input data file
push,lines,file_basename(outfile)             ; output picture name, keep in same directory for now
push,lines,''
push,lines,'EXIT'
push,lines,'END_DAOPHOT'
scriptfile = tmpbase+'.sh'
WRITELINE,scriptfile,lines
FILE_CHMOD,scriptfile,'755'o
file_delete,filebase+'.add.log',/allow

; Run the program
SPAWN,[scriptfile,file_basename(filebase)],out,errout,/noshell

; Delete temporary files, absolute filenames
FILE_DELETE,[scriptfile,tcatfile,tmpbase,tmpbase+'.log'],/allow
if keyword_set(blank) then FILE_DELETE,inptmpfile,/allow

; Back to original directory
CD,curdir

; If outfile exists, move to it's originally intended location
info = file_info(dir+file_basename(outfile))
if info.exists eq 1 and info.size gt 0 then begin
  FILE_MOVE,dir+file_basename(outfile),outfile,/allow  ; nothing will happen if they are the same filename

; Output file not found
endif else begin
  error = outfile+' NOT FOUND'
  if not keyword_set(silent) then print,error
  return
endelse

; There was an error
if errout ne '' then begin
  error = 'There was an error in running DAOPHOT/ADDSTAR. '+errout
  if not keyword_set(silent) then print,error
  return
endif

if keyword_set(stp) then stop

end
