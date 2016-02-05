pro photred_mkopt,input,hilimit=inp_hilimit,va=inp_va,fwhm=fwhm,fitradius_fwhm=inp_fitradius_fwhm,$
                  error=error,verbose=verbose,stp=stp

;+
;
; PHOTRED_MKOPT
;
; This makes opt files for FITS files to be used with
; DAOPHOT and ALLSTAR in the PHOTRED pipeline
;
; INPUTS:
;  input     Input files. Three formats can be used (1) Name of file
;              with a list of input filenames.  Must start with an '@';
;              (2) A name with wildcard characters, such as '*';
;              (3) An array of filenames.
;  =hilimit  The saturation upper limit, 64,000 by default.
;  =va       The spatial variable PSF setting to use
;  =fitradius_fwhm  The value to use for the fitting radius (FI), in
;                      units of the FWHM.
;  /verbose Output information about what is happening.
;  /stp     Stop at the end of the program.
;
; OUTPUTS:
;  Makes .opt and .als.opt files for each FITS file, in the same
;  directory that the FITS file is in.
; 
;  =fwhm    The image FWHM
;  =error   The error message if anything went wrong.
;
; EXAMPLE:
;  IDL>photred_mkopt,'mkopt.lst'
;
; Very similar to Tony Sohn's mkopt.f fortran program
; but estimate FWHM automatically with IMFWHM.PRO
; and gets RDNOISE and GAIN directly from the image
; headers.
;
; By D.Nidever  May 2008   basically a copy of MKOPT.PRO
;                           which was copied from Tony's mkopt
;                           fortran program
;-

undefine,error

; Not enough inputs
ninput = n_elements(input)
if ninput eq 0 then begin
  print,'Syntax -  photred_mkopt,input,hilimit=hilimit,stp=stp,fwhm=fwhm,'
  print,'                        verbose=verbose,error=error'
  return
endif


; Loading input
LOADINPUT,input,files,count=nfiles

; Not enough inputs
if nfiles eq 0 then begin
  print,'No files'
  return
endif

; More than one name input
if nfiles gt 1 then begin
  fwhm = fltarr(nfiles)
  for i=0,nfiles-1 do begin
    PHOTRED_MKOPT,files[i],hilimit=inp_hilimit,va=inp_va,fitradius_fwhm=inp_fitradius_fwhm,fwhm=fwhm1,verbose=verbose
    fwhm[i] = fwhm1
    if keyword_set(verbose) then print,''
  endfor
  return
endif

file = strtrim(files[0],2)

; Default settings
if n_elements(inp_hilimit) gt 0 then hilimit=inp_hilimit else hilimit=6.4e4
if n_elements(inp_va) gt 0 then VA=inp_va>0 else VA=2 ; PSF varies quadratically in the frame
if n_elements(inp_fitradius_fwhm) gt 0 then begin
  if inp_fitradius_fwhm gt 0.0 then fitradius_fwhm=inp_fitradius_fwhm else fitradius_fwhm=1.0  ; can't be >=0
endif else fitradius_fwhm=1.0


; Processing ONE file
;--------------------

test = file_test(file)
if test eq 0 then begin
  print,file,' NOT FOUND'
  fwhm = 99.99
  return
endif

if keyword_set(verbose) then print,'Running PHOTRED_MKOPT on ',file


base = strtrim(FILE_BASENAME(file,'.fits'),2)
dir = FILE_DIRNAME(file)



; Get the FITS header
head = HEADFITS(file,errmsg=errmsg)
if errmsg ne '' then begin
  print,'Error reading header for ',file
  print,errmsg
  error = errmsg
  return
endif


; We need GAIN, READNOISE, FWHM, and Hi-limit
;--------------------------------------------

; Getting GAIN
gain = PHOTRED_GETGAIN(file)
if gain lt 0 then begin
  print,'GAIN NOT found in header'
  error = 'GAIN NOT found in header'
  return
endif

; Getting READNOISE
rdnoise = PHOTRED_GETRDNOISE(file)
if rdnoise lt 0 then begin
  print,'READNOISE NOT found in header'
  error = 'READNOISE NOT found in header'
  return
endif

; Run IMFWHM to get the FWHM
undefine,im,fwhm
IMFWHM,file,fwhm,im=im,/silent

if fwhm gt 90.0 then begin
  print,'Error with FWHM'
  error = 'ERROR with FWHM'
  return
endif

; Load the image
message=''
FITS_READ,file,im,head,/no_abort,message=message
  
; Fits_read error   
if (message ne '') then begin
  error = 'ERROR reading '+file
  print,error
  return
endif


; Getting saturation limit from the header
lolimit = 10000.0                                   ; just in case
saturate = SXPAR(head,'SATURATE',count=nsaturate,/silent)
;if nsaturate eq 0 then saturate=(max(im) < hilimit)  ; if not found
if nsaturate eq 0 then saturate=lolimit > (max(im)-1000) < hilimit  ; if not found

; Minimum of all saturation levels
;;hi = lolimit > ( (saturate - 4000.0) < hilimit )
;hi = lolimit > ( (saturate - 1000.0) < hilimit )
; Don't constrain the saturation value that is input
hi = saturate

; Verbose output
if keyword_set(verbose) then begin
  print,'gain = ',strtrim(gain,2)
  print,'rdnoise = ',strtrim(rdnoise,2)
  print,'fwhm = ',strtrim(fwhm,2)
  print,'saturation = ',strtrim(hi,2)
endif




;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;% MAKING THE OPT FILES


; (1) DAOPHOT parameters
;
; LO    : Low good datum (7. works fine on most imags)
; TH    : Threshold (3.5 works fine)
; LS,HS : Low and high sharpness (default : 0.2 - 1.0)
; LR,HR : Low roundness and high roundness (default : -1.0 - 1.0)
; WA    : Watch progress
; VA    : Variable PSF
; AN    : Analytic model PSF
; EX    : Extra PSF cleaning passes
; PE    : Percent error
; PR    : Profile error

; (2) ALLSTAR parameters
;
; CR    : Clipping range (leave it)
; CE    : Clipping exponent (leave it)
; MA    : Maximum group size
; RED   : Redetermine centroid (0 = no, 1 = yes)

; Frame-specific parameters.
;
; GA    : gain (e/ADU)
; RD    : readout noise (e)
; RE    : readout noise (ADU)
; FW    : FWHM
; HI    : hi good datum in ADU - saturation level
; FI    : fitting radius
; PS    : PSF radius
; IS,OS : inner and outer sky annalus

LO =  7.0
TH =  3.5
LS =  0.2
HS =  1.0
LR = -1.0
HR =  1.0
WA = -2
; VA  defined above
AN = -6     ; It will try all PSF models (#1-6) and use the one with the lowest chi value
EX =  5     ; extra PSF passes
PE =  0.75
PR =  5.00
CR =  2.5
CE =  6.0
MA = 50.
RED = 1.0
WA2 = 0.0


; Frame specific parameters
GA = gain
RD = rdnoise
FW = fwhm
;HI = hi


; Calculating some things
FW = FW < 20             ; daophot won't accept anything higher than 20
RE = RD/GA
FI = fitradius_fwhm*FW < 51                    ; daophot won't accept anything higher than 51
PS = (4.0*FW) < 51       ; daophot won't accept anything higher than 51
IS = (FI - 1.0) < 35     ; daophot won't accept anything higher than 35
OS = (PS + 1.0) < 100    ; daophot won't accept anything higher than 100

; Writing the DAOPHOT parameter
;------------------------------
;
; RE    : readout noise (ADU)
; GA    : gain (e/ADU)
; LO    : Low good datum (7. works fine on most imags)
; HI    : hi good datum in ADU - saturation level
; FW    : FWHM
; TH    : Threshold (3.5 works fine)
; LS,HS : Low and high sharpness (default : 0.2 - 1.0)
; LR,HR : Low roundness and high roundness (default : -1.0 - 1.0)
; WA    : Watch progress
; FI    : fitting radius
; PS    : PSF radius
; VA    : Variable PSF
; AN    : Analytic model PSF
; EX    : Extra PSF cleaning passes
; PE    : Percent error
; PR    : Profile error

outarr = [RE,GA,LO,HI,FW,TH,LS,HS,LR,HR,WA,FI,PS,VA,AN,EX,PE,PR]
anotarr = ['RE','GA','LO','HI','FW','TH','LS','HS','LR','HR','WA','FI','PS','VA','AN','EX','PE','PR']
anotarr = anotarr+' = '
nanot = n_elements(anotarr)

openw,unit,/get_lun,dir+'/'+base+'.opt'
for j=0,nanot-1 do begin
  form = '(A5,F8.2)'
  if anotarr[j] eq 'HI = ' then form = '(A5,I8)'
  printf,unit,format=form,anotarr[j],outarr[j]
endfor
printf,unit,''
close,unit
free_lun,unit


; Writing the ALLSTAR parameter file
;-----------------------------------

; FI    : fitting radius
; IS    :  ??
; OS    :  ??
; RED   : Redetermine centroid (0 = no, 1 = yes)
; WA2   : Watch progress
; PE    : Percent error
; PR    : Profile error
; CR    : Clipping range (leave it)
; CE    : Clipping exponent (leave it)
; MA    : Maximum group size

outarr2 = [FI,IS,OS,RED,WA2,PE,PR,CR,CE,MA]
anotarr2 = ['FI','IS','OS','RE','WA','PE','PR','CR','CE','MA']
anotarr2 = anotarr2+' = '
nanot2 = n_elements(anotarr2)
form = '(A5,F8.2)'

openw,unit,/get_lun,dir+'/'+base+'.als.opt'
for j=0,nanot2-1 do printf,unit,format=form,anotarr2[j],outarr2[j]
printf,unit,''
close,unit
free_lun,unit

; Verbose output
if keyword_set(verbose) then begin
  print,'Created ',dir+'/'+base+'.opt'
  print,'Created ',dir+'/'+base+'.als.opt'
endif


if keyword_set(stp) then stop

end
