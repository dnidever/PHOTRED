;+
;
; PHOTRED_ASTROM_SINGLE
;
; This gets coordinates from the FITS WCS and puts them
; in the photometry file.
;
; INPUTS:
;  longfile    Full path to the file (mch or mag).
;  =catformat  Type of format to use for output catalog.
;  =ndetmin    Minimum number of detections to impose.
;  /redo       Redo files that were already done.
;  /stp        Stop at the end of the program.
;
; OUTPUTS:
;  The final calibrated photometry with accurate astrometry
;
; By D.Nidever  Mar 2008
;-

pro photred_astrom_single,longfile,catformat=catformat,ndetmin=ndetmin,$
                          redo=redo,stp=stp

  ;; Defaults
  if n_elements(catformat) eq 0 then catformat='FITS'

  file = FILE_BASENAME(longfile)
  filedir = FILE_DIRNAME(longfile)

  print,''
  print,'=========================================='
  print,'Getting coordinates for '+longfile
  print,'=========================================='


  CD,current=origdir
  CD,filedir

  ending = first_el(strsplit(longfile,'.',/extract),/last)

  ; Wrong input file
  if (ending ne 'mch' and ending ne 'mag') then begin
    print,file+' ERROR - Input files must end in ".mag" or ".mch"'
    CD,origdir
    return
  endif

  ; ALLFRAME output
  ;------------------
  If (ending eq 'mag') then begin
    base = FILE_BASENAME(file,'.mag')
    mchfile = base+'.mch'
    magfile = base+'.mag'
    ; Switched to the stacked/comb file for the WCS, 10/23/16
    ;  in the original combination procedure these two files ahd
    ;  the identical WCS, in the new combination procedure they
    ;  are very different but the "reference frame" is the 
    ;  combined frame.
    ;fitsfile = base+'.fits'
    fitsfile = base+'_comb.fits'
    photfile = magfile
  Endif ; 'mag' file

  ; DAOPHOT output
  ;------------------
  If (ending eq 'mch') then begin
    base = FILE_BASENAME(file,'.mch')
    mchfile = file
    rawfile = base+'.raw'
    fitsfile = base+'.fits'
    if file_test(fitsfile) eq 0 then fitsfile=base+'.fits.fz'
    photfile = rawfile
  Endif

  ;; Check if the output ast file already exists
  astfile = filedir+'/'+base+'.ast'
  if file_test(astfile) eq 1 and not keyword_set(redo) then begin
    print,astfile,' already exists and /redo not set'
    CD,origdir
    return
  endif

  ; Check that the MAG/RAW, MCH and FITS files exist
  ;-------------------------------------------------
  mchtest = FILE_TEST(mchfile)
  if mchtest eq 1 then nmchlines = FILE_LINES(mchfile) else nmchlines=0
  phottest = FILE_TEST(photfile)
  if phottest eq 1 then nphotlines = FILE_LINES(photfile) else nphotlines=0
  fitstest = FILE_TEST(fitsfile)
  if (nmchlines eq 0) or (nphotlines eq 0) or (fitstest eq 0) then begin
    PUSH,failurelist,longfile
    if mchtest eq 0 then print,mchfile+' NOT FOUND'
    if mchtest eq 1 and nmchlines eq 0 then print,mchfile+' HAS 0 LINES'
    if phottest eq 0 then print,photfile+' NOT FOUND'
    if phottest eq 1 and nphotlines eq 0 then print,photfile+' HAS 0 LINES'
    if fitstest eq 0 then print,fitsfile+' NOT FOUND'
    CD,origdir
    return
  endif


  ; Load the MCH file
  ;------------------
  LOADMCH,mchfile,alsfiles
  nalsfiles = n_elements(alsfiles)    
  numobs = nalsfiles

  ; Load the photometry file
  ;-------------------------  
  phot0 = PHOTRED_READFILE(photfile)
  nphot = n_elements(phot0)
  schema = phot0[0]
  struct_assign,{dum:''},schema
  schema = CREATE_STRUCT(schema,'RA',0.0d0,'DEC',0.0d0)
  phot = REPLICATE(schema,nphot)
  struct_assign,phot0,phot
  undefine,phot0
  print,'Nstars = '+strtrim(nphot,2)

  ; Load the FITS header
  if strmid(fitsfile,6,7,/reverse_offset) eq 'fits.fz' then begin
    head = PHOTRED_READFILE(fitsfile,exten=1,/header)
    ; Fix the NAXIS1/2 values in the header
    sxaddpar,head,'NAXIS1',sxpar(head,'ZNAXIS1')
    sxaddpar,head,'NAXIS2',sxpar(head,'ZNAXIS2')
  endif else begin
    head = PHOTRED_READFILE(fitsfile,/header)
  endelse

  ; Checking that the header has a WCS
  EXTAST,head,astr
  nastr = n_elements(astr)
  if (nastr eq 0) then begin
    print,fitsfile,' DOES NOT HAVE A WCS'
    CD,origdir
    return
  endif

  ; Check that the structure has X/Y
  tags = TAG_NAMES(phot)
  xgd = where(tags eq 'X',nxgd)
  ygd = where(tags eq 'Y',nygd)
  if (nxgd eq 0) or (nygd eq 0) then begin
    print,file,' DOES NOT HAVE X/Y COORDINATES'
    CD,origdir
    return
  endif


  ; Converting to IDL X/Y convention, starting at (0,0)
  ; DAOPHOT has X/Y start at (1,1)
  x = phot.x - 1.0
  y = phot.y - 1.0


  ; Get RA/DEC coordinates for X/Y
  HEAD_XYAD,head,x,y,ra,dec,/degree

  ;; Add the RA/DEC tags
  ;ragd = where(tags eq 'RA',nragd)
  ;if (nragd) eq 0 then ADD_TAG,phot,'RA',0.0d0,phot
  ;decgd = where(tags eq 'DEC',ndecgd)
  ;if (ndecgd) eq 0 then ADD_TAG,phot,'DEC',0.0d0,phot

  ; Put RA/DEC into the structure
  phot.ra = ra
  phot.dec = dec


  ; Apply minimum number of detections, NDETMIN
  if n_elements(ndetmin) gt 0 then begin
    print,'Applying NDETMIN = '+strtrim(ndetmin,2)
    ; ID X Y MAG1 MAG1ERR MAG2 MAG2ERR MAG3 MAG3ERR MAG4 MAG4ERR
    magind = where(stregex(tags,'^MAG',/boolean) eq 1 and stregex(tags,'ERR$',/boolean) eq 0,nmagind)
    print,strtrim(nmagind,2)+' magnitude columns'
    ndet = lonarr(n_elements(phot))
    for i=0,nmagind-1 do ndet += (phot.(magind[i]) lt 50)
    gddet = where(ndet ge ndetmin,ngddet)
    if ngddet gt 0 then begin
      print,'Keeping '+strtrim(ngddet,2)+' of '+strtrim(n_elements(phot),2)+' sources'
    endif else begin
      print,'NO SOURCES PASS.  Keeping the first source.'      
      gddet = 0
      ngddet = 1
    endelse
    phot = phot[gddet]
  endif
  
  
  ; Output the structure to the AST file
  astfile = base+'.ast'
  print,'File with RA/DEC coordinates is: ',filedir+'/'+astfile

  if catformat eq 'FITS' then begin
    MWRFITS,phot,astfile,/create,/silent
  endif else begin   ; ASCII
    PRINTSTR,phot,astfile,/silent
  endelse
    
  ; Check that the file AST file is there
  asttest = FILE_TEST(astfile)
  if (asttest eq 0) then begin
    print,astfile,' NOT FOUND'
  endif

  if keyword_set(stp) then stop

  CD,origdir

end
