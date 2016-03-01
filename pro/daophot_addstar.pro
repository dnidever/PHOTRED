;+
;
; DAOPHOT_ADDSTAR
;
; Add a fake star to an image
;
; INPUTS:
;  psffile   The DAOPHOT PSF file.
;  cat       The input catalog structure of sources to add.  It must have
;              AT LEAST these columns:
;              X, Y     x/y-coordinates of sources in IRAF/DAOPHOT format (1-indexed)
;              MAG      magnitude of star
;  outfile   The name of the output FITS file.  If it already exists
;              then the star is added to the existing image (unless /clobber set).
;  /clobber  Delete the output file if it already exists.
;  /silent   Don't print anything to the screen.
;
; OUTPUTS:
;  The fake PSF is added to the output image.  If the output file
;  already exists then the sources are added to it (unless /clobber set).
;  =error     The error message if one occurred, otherwise undefined.
;
; USAGE:
;  IDL>daophot_addstar,'F1-0087651_01.psf',cat,'fake.fits',error=error
;
; By D.Nidever  March 2016
;-

pro daophot_addstar,filebase,cat,outfile,clobber=clobber,error=error,silent=silent

undefine,error

; Do we have enough inputs
if n_elements(filebase) eq 0 or n_elements(cat) eq 0 or n_elements(outfile) eq 0 then begin
  error = 'Not enough inputs'
  print,'Syntax -  daophot_addstar,filebase,cat,outfile,clobber=clobber,error=error,silent=silent'
  return
endif

; Make the psf, opt and fits filenames
psffile = filebase+'.psf'
optfile = filebase+'.opt'
fitsfile = filebase+'.fits'

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
ny = sxpar(head,'NAXIS1')
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
  if strpos(tags,reqtags[i]) eq -1 then
    error = 'CAT must have '+reqtags[i]+' column'
    if not keyword_set(silent) then print,error
    return
  endif
endfor

; Make sure daophot.opt file exists
if FILE_TEST('daophot.opt') eq 0 then FILE_COPY,optfile,'daophot.opt',/over

; Base for temporary files
tmpbase = MKTEMP('psf')

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
;push,lines,'rm ${image}.temp.log      >& /dev/null'
;push,lines,'rm ${image}.temp.coo      >& /dev/null'
;push,lines,'rm ${image}.temp.ap       >& /dev/null'
push,lines,'daophot << END_DAOPHOT >> ${image}.add.log'
push,lines,'OPTIONS'
push,lines,'${image}.opt'
push,lines,''
push,lines,'ATTACH ${image}.fits'
push,lines,'ADDSTAR'
push,lines,'${image}.psf'
push,lines,'1'
push,lines,stringize(gain,ndec=5) ; photons per ADU
push,lines,tmpfile                ; output picture name
push,lines,''
push,lines,'EXIT'
push,lines,'END_DAOPHOT'
scriptfile = tmpbase+'.sh'
WRITELINE,scriptfile,lines
FILE_CHMOD,scriptfile,'755'o

; Run the program
SPAWN,[scriptfile,filebase],out,errout,/noshell
;FILE_DELETE,scriptfile   ; delete the temporary script

;; Test the coo and ap file
;cootest = FILE_TEST(base+'.temp.coo')
;if cootest eq 1 then coolines=FILE_LINES(base+'.temp.coo') else coolines=0
;aptest = FILE_TEST(base+'.temp.ap')
;if aptest eq 1 then aplines=FILE_LINES(base+'.temp.ap') else aplines=0

; Output file, are we deleting?

  
stop

end
