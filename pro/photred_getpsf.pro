;+
;
; PHOTRED_GETPSF
;
; This runs the DAOPHOT routines to create the PSF for an image.
;
; INPUTS:
;  base      The base name of the FITS and MCH files.
;  /fake     Run for artificial star tests.
;  =logfile  A logfile to print the output to.
;
; OUTPUTS:
;  The DAOPHOT psf file with name BASE.psf.
;  =error    The error message if one occurred.
;
; USAGE:
;  IDL>photred_getpsf,'F1-23430911_10'
;
; By D.Nidever  Feb 2019
;-

pro photred_getpsf,base,fake=fake,error=error,logfile=logfile

;; Not enough inputs
if n_elements(base) eq 0 then begin
  error = 'Not enough inputs'
  print,'Syntax - photred_getpsf,base,fake=fake,error=error,logfile=logfile'
  return
endif

if n_elements(logfile) eq 0 then logfile=-1

;; Make SExtractor output parameter file
secols = ['NUMBER','X_IMAGE','Y_IMAGE','MAG_APER(1)','MAGERR_APER(1)','MAG_AUTO','MAGERR_AUTO','BACKGROUND',$
          'THRESHOLD','ISOAREA_IMAGE','A_IMAGE','B_IMAGE','THETA_IMAGE','ELLIPTICITY','FWHM_IMAGE','FLAGS',$
          'IMAFLAGS_ISO(1)','NIMAFLAGS_ISO(1)','CLASS_STAR']
WRITELINE,base+'.param',secols

;; Create the SExtractor config file
PHOTRED_MKSEXCONFIG,base+'.fits',base+'.sex',base+'.cat',flagfile=base+'.bpm.fits',$
                    paramfile=base+'.param',cattype='FITS_1.0'

;; Run SExtractor for detection
SPAWN,['sex',base+'.fits','-c',base+'.sex'],out,errout,/noshell
hd = PHOTRED_READFILE(base+'.cat',exten=1,/header)
nsources = sxpar(hd,'NAXIS',1)
if nsources lt 1 then begin
  error = 'Only '+strtrim(nsources,2)+' sources. Need at least 1 to create a PSF'
  printlog,logfile,error
  return
endif

;; Apply cuts to get good stars
sex = MRDFITS(base+'.cat',1,/silent)
if n_tags(sex) eq 1 then sex = MRDFITS(catfile,2,/silent)
nsex = n_elements(sex)
si = sort(sex.fwhm_image)
fwhm80 = sex[si[floor(nsex*0.80)]].fwhm_image
bd = where(sex.imaflags_iso gt 0 or sex.class_star lt 0.5 or $
           sex.ellipticity gt 0.8 or sex.fwhm_image gt fwhm80 or $
           ((sex.flags and 8) eq 8) or ((sex.flags and 16) eq 16) or $
           1.087/sex.magerr_aper le 10,nbd,comp=gd,ncomp=ngd)
;; Not enough stars, remove class_star cut and raise S/N cut
if ngd lt 50 then $
  bd = where(sex.imaflags_iso gt 0 or $
             sex.ellipticity gt 0.8 or sex.fwhm_image gt fwhm80 or $
             ((sex.flags and 8) eq 8) or ((sex.flags and 16) eq 16) or $
             1.087/sex.magerr_aper le 7,nbd,comp=gd,ncomp=ngd)
;; No sources left
if ngd eq 0 then begin
  error = 'No sources left after stellar cuts'
  print,error
  return
endif
sex = sex[gd]
nsources = ngd
MWRFITS,sex,base+'.cat',/create

;; Convert to DAOPHOT coo format
SEX2DAOPHOT,base+'.cat',base+'.fits',base+'.coo'

;; Check that there are enough stars for our VA setting
;;   VA=2 quadratic spatial variations, 6 stars minimum
;;   VA=1 linear spatial variations, 3 stars minimum
;;   VA=0 empirical corrections but no spatial variations, 1 star
;;   minimum
;;   VA=-1 analytical, 1 star minimum
READLINE,base+'.opt',optlines
optarr = strtrim(strsplitter(optlines,'=',/extract),2)
vaind = where(reform(optarr[0,*]) eq 'VA',nvaind)
vaval = float(reform(optarr[1,vaind]))
minsources = [1,1,3,6]
if nsources lt minsources[long(vaval)+1] then begin
  ;; Get new VA value, largest allowed for this number of sources
  if nsources ge 3 then newvaval=1.0 else newvaval=0.0
  printlog,logfile,'Only '+strtrim(nsources,2)+' sources and need '+strtrim(minsources[long(vaval)+1],2)+' for VA='+strtrim(vaval,2)+'. Lower to VA='+strtrim(newvaval,2)
  newoptlines = optlines
  newoptlines[vaind] = 'VA = '+string(newvaval,format='(F8.2)')
  WRITELINE,base+'.opt',newoptlines
endif

;; Sometimes the filenames get too long for DAOPHOT
;; use temporary files and symlinks
tbase = (file_basename(MKTEMP('cmb',/nodot)))[0]  ; create base, leave so other processes won't take it
tfits = tbase+'.fits'    &  file_delete,tfits,/allow  &  file_link,base+'.fits',tfits
topt = tbase+'.opt'      &  file_delete,topt,/allow   &  file_link,base+'.opt',topt
taopt = tbase+'.als.opt' &  file_delete,taopt,/allow  &  file_link,base+'.als.opt',taopt
tcoo = tbase+'.coo'      &  file_delete,tcoo,/allow   &  file_link,base+'.coo',tcoo


;; Get the PSF of the combined image
SPAWN,['./getpsfnofind.sh',tbase],/noshell

;; If getpsf failed, change to VA=0
info = file_info(base+'.psf')
if info.exists eq 0 or info.size eq 0 then begin
  READLINE,base+'.opt',optlines
  vaind = where(strmid(optlines,0,2) eq 'VA',nvaind)
  vaval = float(reform(optarr[1,vaind]))
  newoptlines = optlines
  newoptlines[vaind] = 'VA = '+string(0.0,format='(F8.2)')
  WRITELINE,base+'.opt',newoptlines
  printlog,logfile,'getpsfnofind.sh failed.  Changing to VA=0.  Trying again.'
  SPAWN,['./getpsfnofind.sh ',tbase],/noshell
endif

;; Delete the temporary symlinks
FILE_DELETE,[tbase,tfits,topt,taopt,tcoo],/allow

;; No PSF file found
info = file_info(base+'.psf')
if info.exists eq 0 or info.size eq 0 then begin
  error = 'Could not create PSF for '+base
  printlog,logfile,error
  return
endif

end
