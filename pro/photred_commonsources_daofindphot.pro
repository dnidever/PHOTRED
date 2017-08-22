;+
;
; PHOTRED_COMMONSOURCES_DAOFINDPHOT
;
; Find sources and get simple aperture photometry using DAOPHOT.
;
; INPUTS:
;  file      The FITS file name.
;  /clobber  Overwrite any existing output files.
;
; OUTPUTS:
;  The DAOPHOT "find" .cmn.coo and "photometry" .cmn.ap files
;   are created.
;  =error    The error message if one occurred.
;
; USAGE:
;  IDL>photred_commonsources_daofindphot,'F1-00126342_10.fits',/clobber
;
; By D.Nidever  May 2015
;-

pro photred_commonsources_daofindphot,file,clobber=clobber,error=error

; Not enough inputs
if n_elements(file) eq 0 then begin
  print,'syntax - photred_commonsources_daofindphot,file,clobber=clobber,error=error'
  error = 'Not enough inputs'
  return
endif

; Is this a FITS file
if (strmid(file,3,4,/reverse_offset) ne 'fits' and strmid(file,6,7,/reverse_offset) ne 'fits.fz') then begin
  print,file,' is NOT a FITS file'
  error = file+' is NOT a FITS file'
  return
endif

; Does the FITS file exist
if FILE_TEST(file) eq 0 then begin
  print,file,' NOT FOUND'
  error = file+' NOT FOUND'
  return
endif

if strmid(file[0],6,7,/reverse_offset) eq 'fits.fz' then begin
  fpack = 1
  base = FILE_BASENAME(file[0],'.fits.fz')
endif else begin
  fpack = 0
  base = FILE_BASENAME(file[0],'.fits')
endelse
dir = FILE_DIRNAME(file[0])

cd,current=curdir
cd,dir  ; go to the directory in which we need to be

; Do the CMN.COO and CMN.AP files already exist?
coofile = base+'.cmn.coo'
cootest = FILE_TEST(coofile)
if cootest eq 1 then coolines=FILE_LINES(coofile) else coolines=0
apfile = base+'.cmn.ap'
aptest = FILE_TEST(apfile)
if aptest eq 1 then aplines=FILE_LINES(apfile) else aplines=0

; Run DAOPHOT FIND and PHOTOMETRY
if (coolines lt 4 or aplines lt 4 or keyword_set(clobber)) then begin

  ; RUN DAOPHOT
  print,'Getting sources for ',base,' using DAOPHOT'

  ; Make a .opt file
  optfile = base+'.opt'
  opttest = FILE_TEST(optfile)
  if opttest eq 0 then begin
    print,'No opt files.  Making them.'
    PHOTRED_MKOPT,file,error=mkopterror
    ; redo the test
    optfile = base+'.opt'
    opttest = FILE_TEST(optfile)
  endif

  ; Problems with .opt file
  if FILE_TEST(optfile) eq 0 then begin
    print,'No opt file.  Cannot proceed.'
    goto,BOMB_DAOPHOT
  endif

  ; Load the opt file
  READLINE,optfile,optlines
  fwhmind = where(stregex(optlines,'FW =',/boolean) eq 1,nfwhmind)
  fwhmarr = strsplit(optlines[fwhmind[0]],'=',/extract)
  fwhm = float(fwhmarr[1])
 
  ; Copy the .opt file daophot.opt 
  if FILE_TEST('daophot.opt') eq 0 then $
    FILE_COPY,base+'.opt','daophot.opt',/over,/allow
  
  ; Make temporary photo.opt file
  tphotofile = MKTEMP('photo')
  undefine,photlines
  push,photlines,'A1 = '+STRING(fwhm*3.0,format='(F7.4)')
  push,photlines,'IS = 45.0000'
  push,photlines,'OS = 50.0000'
  WRITELINE,tphotofile,photlines
 
  ; Make a temporary script to run FIND
  undefine,lines
  push,lines,'#!/bin/sh'
  push,lines,'export image=${1}'
  push,lines,'rm ${image}.cmn.log      >& /dev/null'
  push,lines,'rm ${image}.cmn.coo      >& /dev/null'
  push,lines,'rm ${image}.cmn.ap       >& /dev/null'
  push,lines,'daophot << END_DAOPHOT >> ${image}.cmn.log'
  push,lines,'OPTIONS'
  push,lines,'${image}.opt'
  push,lines,''
  push,lines,'ATTACH ${image}.fits'
  push,lines,'FIND'
  push,lines,'1,1'
  push,lines,'${image}.cmn.coo'
  push,lines,'y'
  push,lines,'PHOTOMETRY'
  push,lines,FILE_BASENAME(tphotofile)
  push,lines,''
  ; Erase any PSF file for this frame
  psffile = base+'.psf'
  if FILE_TEST(psffile) eq 1 then begin
    FILE_DELETE,psffile,/allow,/quiet
    ;if FILE_LINES(psffile) gt 0 then $
    ;  push,lines,'END-OF-FILE'  ; don't use profile-fitting
  endif
  push,lines,'${image}.cmn.coo'
  push,lines,'${image}.cmn.ap'
  push,lines,'EXIT'
  push,lines,'END_DAOPHOT'
  ;tempfile = maketemp('dao','.sh')
  tempfile = MKTEMP('dao')    ; absolute path
  WRITELINE,tempfile,lines
  FILE_CHMOD,tempfile,'755'o

  ; If this is a fpack-compressed file then temporarily uncompress it
  if fpack eq 1 then spawn,['funpack',ifile],/no_shell

  ; Run the program
  SPAWN,[tempfile,base],out,errout,/noshell

  ; Delete funpacked file
  if fpack eq 1 then FILE_DELETE,ibase+'.fits',/allow

  ; Test the coo and ap files
  cootest = FILE_TEST(base+'.cmn.coo')
  if cootest eq 1 then coolines=FILE_LINES(base+'.cmn.coo') else coolines=0
  if coolines gt 0 then print,base+'.cmn.coo created with '+strtrim(coolines,2)+' lines' else $
    print,base+'.cmn.coo NOT FOUND'
  aptest = FILE_TEST(base+'.cmn.ap')
  if aptest eq 1 then aplines=FILE_LINES(base+'.cmn.ap') else aplines=0
  if aplines gt 0 then print,base+'.cmn.ap created with '+strtrim(aplines,2)+' lines' else $
    print,base+'.cmn.ap NOT FOUND'

  ; Remove the temporary files
  FILE_DELETE,tempfile,/allow    ; delete the temporary script
  ;FILE_DELETE,'daophot.opt',/allow
  ; LEAVE IT THERE OTHERWISE THERE COULD BE "RACE" ISSUES
  FILE_DELETE,tphotofile,/allow
  junk = FILE_TEST(base+'jnk.fits')
  if junk eq 1 then FILE_DELETE,base+'jnk.fits'

; Files already exist and /clobber not set
endif else begin
  print,'Photometry files already exist and /clobber NOT set'
endelse

BOMB_DAOPHOT:

cd,curdir  ; back to original directory

;stop

end
