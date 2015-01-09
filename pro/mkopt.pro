pro mkopt,input,hilimit=hilimit,satlevel=satlevel,stp=stp,fwhm=fwhm

;+
;
; MKOPT
;
; This makes opt files for FITS files to be used with
; DAOPHOT and ALLSTAR
;
; INPUTS:
;  input     Input files. Three formats can be used (1) Name of file
;            with a list of input filenames.  Must start with an '@';
;            (2) A name with wildcard characters, such as '*';
;            (3) An array of filenames.
;  =satlevel The saturation level to use.
;  =hilimit  The saturation upper limit, 64,000 by default.
;  /stp      Stop at the end of the program.
;
; OUTPUTS:
;  Makes .opt and .als.opt files for each FITS file
;  =fwhm   The image FWHM
;
; EXAMPLE:
;  IDL>mkopt,'mkopt.lst'
;
; Very similar to Tony Sohn's mkopt.f fortran program
; but estimate FWHM automatically with IMFWHM.PRO
; and gets RDNOISE and GAIN directly from the image
; headers.
;
; By D.Nidever  Oct. 2006  (copied from Tony's mkopt fortran program)
;-

; (1) DAOPHOT parameters
;
; LO    : Low good datum (7. works fine on most imags)
;           this is number of standard deviations below the mean sky level
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
VA =  2     ; PSF varies quadratically in the frame
AN = -6     ; It will try all PSF models (#1-6) and use the one with the lowest chi value
EX =  5     ; extra PSF passes
PE =  0.75
PR =  5.00
CR =  2.5
CE =  6.0
MA = 50.
RED = 1.0
WA2 = 0.0

if n_elements(hilimit) eq 0 then hilimit = 6.4e4
lolimit = 10000.0       ; just in case

; Not enough inputs
ninput = n_elements(input)
if ninput eq 0 then begin
  print,'Syntax -  mkopt,input,hilimit=hilimit,stp=stp,fwhm=fwhm'
  return
endif

; Loading input
LOADINPUT,input,files
nfiles = n_elements(files)

; Not enough inputs
if nfiles eq 0 then begin
  print,'No files'
  return
endif

; Starting arrays
filarr = strarr(nfiles)
gaarr = fltarr(nfiles)
rdarr = fltarr(nfiles)
fwarr = fltarr(nfiles)
hiarr = fltarr(nfiles)


;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;% GETTING THE DATA

; Looping through the files
for i=0,nfiles-1 do begin

  dum = file_search(files[i])
  if dum eq '' then begin
    print,files(i),' DOES NOT EXIST'
    goto,skip
  end

  FITS_READ,files[i],im,head
  ;head = headfits(files(i))

  ; We need file name (without extension), exptime, gain, rdnoise
  ; FWHM, Hi-limit and filter

  file = strtrim(file_basename(files(i),'.fits'),2)

  ; Getting GAIN
  gain = photred_getgain(files[i])
  ;gain = sxpar(head,'GAIN')
  ;if gain eq 0 then gain=sxpar(head,'EGAIN')          ; imacs
  if gain le 0 then begin
    print,'gain NOT found in header. Using gain=1.0'
    gain = 1.0
  endif

  ; Getting READNOISE
  rdnoise = photred_getrdnoise(files[i])
  ;rdnoise = sxpar(head,'RDNOISE')
  ;if rdnoise eq 0 then rdnoise=sxpar(head,'READNOIS') ; Swope
  ;if rdnoise eq 0 then rdnoise=sxpar(head,'ENOISE')   ; imacs
  if rdnoise le 0 then begin
    print,'readnoise NOT found in header. Using rdnoise=0.01'
    rdnoise = 0.01
  endif

  ; Run IMFWHM to get the FWHM
  undefine,fwhm
  ;IMFWHM,files[i],fwhm,im=im,/silent
  IMFWHM,files[i],fwhm,/silent

  ; Get Hi-limit
  ;fits_read,files(i),im
  ;hilim = max(im) - 4000.0
  ;hilim = floor(hilim/100.)*100  ; round down to the nearest 100

  ; Getting saturation limit from the header
  saturate = sxpar(head,'SATURATE')
  if saturate eq 0 then saturate=(max(im) < hilimit)  ; if not found

  ; Minimum of all saturation levels
  hi = lolimit > ( (saturate - 4000.0) < hilimit )
  if n_elements(satlevel) gt 0 then hi=satlevel   ; saturation level fixed

  ; Putting data in the arrays
  filarr[i] = file
  gaarr[i] = gain
  rdarr[i] = rdnoise
  fwarr[i] = fwhm
  hiarr[i] = hi

  SKIP:

end


;stop


;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;% USE A MEDIAN VALUE FOR EACH MOSAIC FRAME

filebase = strarr(nfiles)
for i=0,nfiles-1 do $
  filebase[i] = first_el(strsplit(filarr[i],'_',/extract))
ui = uniq(filebase,sort(filebase))
nui = n_elements(ui)
fils = filebase(ui)

fwarr_orig = fwarr

; Looping through the MOSAIC files
;print,''
for i=0,nui-1 do begin
  g = where(filebase eq fils[i] and fwarr lt 40.0,ng)
  gall = where(filebase eq fils[i],ngall)

  if (ng gt 0) then begin
    medfw = median(fwarr[g])
    resistant_mean,fwarr[g],2.0,mnfw,sigma_mnfw
    fwarr[gall] = mnfw            ; every frame gets this value!!
  endif else begin

    fwarr[gall] = 99.99
    mnfw = 99.99
    print,'*** NO GOOD FWHM FOR ',fils[i],' ***'

    ;stop   ; to enter the FWHM manually

  endelse

  print,'FWHM for ',fils[i],' = ',string(mnfw,format='(F7.3)')

end


;stop


;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;% MAKING THE OPT FILES

; Looping through the files
for i=0,nfiles-1 do begin

  file = filarr[i]
  GA = gaarr[i]
  RD = rdarr[i]
  FW = fwarr[i]
  HI = hiarr[i]

  ; Calculating some things
  FW = FW < 20             ; daophot won't accept anything higher than 20
  RE = RD/GA
  FI = FW < 51             ; daophot won't accept anything higher than 51
  PS = (4.0*FW) < 51       ; daophot won't accept anything higher than 51
  IS = (FI - 1.0) < 35     ; daophot won't accept anything higher than 35
  OS = (PS + 1.0) < 100    ; daophot won't accept anything higher than 100

  ; Writing the DAOPHOT parameter
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
  form = '(A5,F8.2)'

  openw,unit,/get_lun,file+'.opt'
  for j=0,nanot-1 do printf,unit,format=form,anotarr[j],outarr[j]
  printf,unit,''
  close,unit
  free_lun,unit


  ; Writing the ALLSTAR parameter file

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

  openw,unit,/get_lun,file+'.als.opt'
  for j=0,nanot2-1 do printf,unit,format=form,anotarr2[j],outarr2[j]
  printf,unit,''
  close,unit
  free_lun,unit

  ;stop

end

if keyword_set(stp) then stop

end
