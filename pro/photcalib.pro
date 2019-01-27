;+
;
; PHOTCALIB
;
; This is a modified version of Mike Siegel's MAGMA program that
; calibrates raw photometry.  MAGMA is very versatile and can
; average multiple exposures per band, offsets, variables stars
; and the more.  PHOTCALIB just calibrates raw photometry
; non-interactively.  If you need more than that then you should
; probably use MAGMA instead.
; This is an updated version of MSCMAGMA.PRO
;
; While the original MAGMA took minutes to run, PHOTCALIB has
; been "vectorized" and runs in just a few seconds.
;
; INPUTS:
;  inpfile    This is an input file that lists information about
;             the raw photometry files.  You can calibrate many
;             raw photometry files at once.  Each file needs its
;             own line in the input file that has the following
;             information:
;
;             raw photometry filename, Band1 name, Band1 airmass,
;             Band1 exptime, Band1 aperture correction, Band2 ...
;             The aperture corrections need to be POSITIVE!!
;
; This is an example inpfile:
; obj1110_1.raw  I  1.7270  90.0  0.0141  M  1.6960  60.0  0.0081  D  1.7040  420.0  0.0079
;
;             The input files can easily be generated with the
;             MSCMAGMA_PREP.PRO program.
;
;             The raw photometry files are assumed to be in
;             the DAOMASTER format:
;             -The first three lines are a header, and are skipped
;             -Each star has a separate line with the format:
;              ID, X, Y, Band1, Band1 error, Band2, Band2 error, ...
;               chi, sharp
;
;             The raw files can have as many bands as you like,
;             but otherwise they must conform to this format.
;
;
;  transfile  This gives the transformation information needed
;             to calibrate the raw photometry.  Normally there
;             is a different one of these for every night.
;
;             There need to be two lines per band.
;             FIRST line:  band name,  color name, transformation
;             equation coefficients (zero point, airmass, color
;             airmass*color, color^2)
;             SECOND line: errors in the transformation equation
;             coefficients
;
;     This is an example transfile:
;     M    M-T  -0.9990    0.1402     -0.1345    0.0000   0.0000
;               1.094E-02  5.037E-03  2.010E-03  0.0000   0.0000
;     T    M-T  -0.0061    0.0489     0.0266     0.0000   0.0000
;               6.782E-03  3.387E-03  1.374E-03  0.0000   0.0000
;     D    M-D  1.3251     0.1403     -0.0147    0.0000   0.0000
;               1.001E-02  5.472E-03  2.653E-02  0.0000   0.0000
;
;  =inptrans  The transformation structure to use.
;  /average  This averages multiple frames in the same band, but
;              also keeps the individual frames.  The tag names will
;              be:  "M1, M1err, M2, M2err, ..., M, Merr"
;  /onlyaverage  This averages multiples frames, but does not keep
;                the individual observations.  (setting /combine
;                  does the same thing).
;  /keepinstrumental  Keep the instrumental magnitudes and errors.
;                     They will have names of "I_MAG","I_MAGERR".
;  /header  The column names are in the first line.
;  =catformat  Catalog format, FITS or ASCII.  ASCII by default.
;  /silent  Don't print anything out
;  /stp     Stop at the end of the program
;
; OUTPUTS:
;  A calibrated photometry file will be output for each raw photometry
;  file in the input file.  The ".phot" extension will be used for
;  the calibrated photometry files.
;
;
; EXAMPLE:
;  IDL>photcalib,'field1.input','n1.trans'
; 
;
; History:  1989:       getmags.for written by Majewski
;           1997-2003:  Magma.pro written by Mike Siegel
;           2007:       Mscmagma.pro by David Nidever
;           Mar 2008:   Renamed photcalib.pro
;-

FUNCTION simplerr,inerr,am,cl,cler,t
;=====================================================================
;
;   This propogates photometric error through the transformation equations
;   This version takes into account error in transformation constants, although
;   it assumes perfect airmass
;
; sigma(V)=sqrt[ sigma(mV)^2 + + sigma(v1)^2 + (XV*sigma(v2))^2+((B-V)*sigma(v3))^2 
;		+ (XV*(B-V)*sigma(v4))^2 + (B-V*sigma(v5))^2]
;                +(v3+v4*XV*+v5*2*(B-V))^2*sigma(B-V)^2               
;
;=====================================================================

temper = inerr^2 + t.zptermsig^2 + (am*t.amtermsig)^2
temper = temper + (cl*t.coltermsig)^2 + (am*cl*t.amcoltermsig)^2
temper = temper + (cl*cl*t.colsqtermsig)^2
temper = temper + (t.colterm+t.amcolterm*am+2*t.colsqterm*cl)^2*cler^2
outerr = sqrt(temper)

; Set bad ones to 9.9999
bd = where(inerr gt 9.,nbd)
if nbd gt 0 then outerr[bd] = 9.9999

return,outerr

end

;------------------------------------------------------------


FUNCTION simplestar,inmag,am,colr,apcorr,exp,t
;=====================================================================
;
;   This uses the solved colors and other terms to find the magnitude in each band
;
;   V = mV - v1 - v2 * XV - v3 * (B-V) - v4 * XV*(B-V) - v5 * (B-V) * (B-V)
;         +(aperture correction) + (time correction)
; 
;   The aperture corrections need to be POSITIVE
;
;=====================================================================

outmag = inmag - t.zpterm - t.amterm*am - t.colterm*colr
outmag = outmag - t.amcolterm*am*colr-t.colsqterm*colr*colr
outmag = outmag - apcorr 
;outmag = outmag + apcorr

; Correct for exposure time
if (exp gt 0) then begin
   outmag = outmag+2.5*alog10(exp)
endif

; Set bad ones to 99.9999
bd = where(inmag gt 90.,nbd)
if nbd gt 0 then outmag[bd] = 99.9999

return,outmag

end

;------------------------------------------------------------

PRO solvestar,instar,trans,inp,outstar
;=====================================================================
;
;   The heart of the program.  This iteratively solves for each star.  It applies 
;   the transformation equation assuming a color of zero.  It then averages the 
;   passbands that are common, solves for the color and resolves for each
;   magnitude given the new color, gradually iterating until convergence.
;   
;   The solved magnitude is in the form: (V in the example)
;   V = mV - v1 - v2 * XV - v3 * (B-V) - v4 * XV*(B-V) - v5 * (B-V) * (B-V)
;         +(aperture correction)+(time correction)
;
;  It sends back to the code outstar, an array of average values in each
;  passband and tempstar, an array of the individual solved magnitudes
;
;=====================================================================

; Setting up some important arrays
numobs = n_elements(inp.band)
numstar = n_elements(instar[*,0])
passband = trans.band
colband = trans.colband
colsign = trans.colsign

; The input magnitudes and errors
inmag = instar[*,2*lindgen(numobs)+3]
inerr = instar[*,2*lindgen(numobs)+4]

; Initializing some arrays
clr = inmag*0.
clrerr = inmag*0.
tempmag = inmag*0.
temperr = inmag*0.
laststar = inmag*0.


;#############################
;# First we set the color terms to zero, then solve for the magnitudes using
;#   simplestar

; Loop through the exposures
for a=0,numobs-1 do begin

  ; Assume an initial color and color error of zero
  clr[*,a] = 0
  clrerr[*,a] = 0

  ; Run simplestar to calculate the tranformed magnitudes
  newmag = SIMPLESTAR(inmag[*,a],inp.airmass[a],clr[*,a],inp.apcorr[a],inp.exptime[a],trans[a])
  tempmag[*,a] = newmag

  ; Run simplerr to calculate the error in the transformed magnitudes
  newerr = SIMPLERR(inerr[*,a],inp.airmass[a],clr[*,a],clrerr[*,a],trans[a])
  temperr[*,a] = newerr

endfor


;##################################
;# Now begin the iteration loop
;##################################

niter = 0
converge = 0

WHILE (converge eq 0) do begin


  ; #############################
  ; First set the color term.
  ; Passbands with an indefinite color will have it set to zero.

  ; Loop through the bands/exposures
  for d=0,numobs-1 do begin

    ; Index of the passband to use for the color
    ind = where(passband eq colband[d],nind)
    ;ind = first_el(where(passband eq colband[d],nind))

    ; More than one color band, average them
    ;----------------------------------------
    if (nind gt 1) then begin

      ; Average the mags
      AVERAGEMAG,tempmag[*,ind],temperr[*,ind],colmag,colerr,/robust

      ; Making the color
      ; Color sign = 1
      if (colsign[d] eq 1) then begin
        clr[*,d] = tempmag[*,d]-colmag

      ; Color sign = 2
      endif else begin
        clr[*,d] = colmag-tempmag[*,d]
      endelse

      ; To avoid the color^2 recursion, color error is taken from instrumental errors
      gd = where(colerr lt 9.0,ngd)
      if ngd gt 0 then clrerr[gd,d]  = colerr[gd]
      ;gd = where(inerr[*,ind] lt 9.0,ngd)
      ;if ngd gt 0 then clrerr[gd,d]  = inerr[gd,ind]

      ; Bad photometry, use color=0
      bd = where( (colmag gt 90.) OR (tempmag[*,d] gt 90.) ,nbd)
      if nbd gt 0 then begin
        clr[bd,d] = 0.0  
        clrerr[bd,d] = 0.0
      endif

    ; One or No color band
    ;----------------------
    endif else begin

      ; We have a valid color index
      if (nind eq 1) then begin
        ind = ind[0]

        ; Color sign = 1
        if (colsign[d] eq 1) then begin
          clr[*,d] = tempmag[*,d]-tempmag[*,ind]
    
        ; Color sign = 2
        endif else begin
          clr[*,d] = tempmag[*,ind]-tempmag[*,d]
        endelse

        ; To avoid the color^2 recursion, color error is taken from instrumental errors
        gd = where(inerr[*,ind] lt 9.0,ngd)
        if ngd gt 0 then clrerr[gd,d]  = inerr[gd,ind]

        ; Bad photometry, use color=0
        bd = where( (tempmag[*,ind] gt 90.) OR (tempmag[*,d] gt 90.) ,nbd)
        if nbd gt 0 then begin
          clr[bd,d] = 0.0
          clrerr[bd,d] = 0.0
        endif

      ; No valid color index, set color=0
      endif else begin
        clr[*,d] = 0.0
        clrerr[*,d] = 0.0
      endelse

    endelse  ; one or no color band

  endfor ; looping through the bands/exposures


  ; ############################
  ; Now resolve the star

  ; Loop through the stars
  for f=0,numobs-1 do begin
    ; Run simplestar to calculate the tranformed magnitudes
    tempmag[*,f] = SIMPLESTAR(inmag[*,f],inp.airmass[f],clr[*,f],inp.apcorr[f],inp.exptime[f],trans[f])
  
    ; Run simplerr to calculate the error in the transformed magnitudes
    temperr[*,f] = SIMPLERR(inerr[*,f],inp.airmass[f],clr[*,f],clrerr[*,f],trans[f])
  endfor


  ; #######################
  ; Check for convergence
  ; the iteration loop recycles until the solution is good
  ; Every star+band must not change at the 0.002 level in order to stop
  ; OR niter>30

  converge = 1  ; assume good at first

  maxiter = 50      ; 30
  maxdiff = 0.0001  ;0.002
  absdiff = abs(tempmag-laststar)
  bd = where( absdiff gt maxdiff,nbd)
  if nbd gt 0 then converge=0

  ; Copying current solution to "last" solution
  if (converge eq 0) then laststar = tempmag

  ; Go up to 50 iterations, send out a warning message if it doesn't converge
  if (niter gt maxiter) then begin
    bd2 = array_indices(tempmag,bd)
    nbdstars = n_elements(bd2[0,*])
    ;print,strtrim(nbdstars,2),'/',strtrim(numstar,2),' failed to converge.'
    print,strtrim(nbdstars,2),'/',strtrim(numstar,2),' failed to converge after ',strtrim(maxiter,2),' iterations.  Max differences of ',strtrim(max(absdiff),2),' mag'
    converge = 1
  endif

  ; Increment
  niter = niter+1

ENDWHILE

; Putting together the output array
outstar = instar*0.
; Transfer over the id,position,chi and sharp
outstar[*,0] = instar[*,0]
outstar[*,1] = instar[*,1]
outstar[*,2] = instar[*,2]
;outstar[*,2*numobs+3] = instar[*,2*numobs+3]
;outstar[*,2*numobs+4] = instar[*,2*numobs+4]
; Transferring over all other columns
outstar[*,2*numobs+3:*] = instar[*,2*numobs+3:*]

; Transfer the final magnitudes and errors
outstar[*,2*lindgen(numobs)+3] = tempmag
outstar[*,2*lindgen(numobs)+4] = temperr

end



;----------------------------------------------------------



PRO  photcalib,inpfile,transfile,silent=silent,stp=stp,average=average,$
               onlyaverage=onlyaverage,keepinstrumental=keepinstrumental,$
               combine=combine,header=header,logfile=logfile,$
               catformat=catformat,inptrans=inptrans

;=====================================================================
;
;   This is the main program.  It reads in the data from magfile and trans
;   after asking for user input, initializes variables and then starts 
;   solving for each star.
;
;=====================================================================

; Not enough inputs
if n_params() lt 2 then begin
  print,'Syntax - photcalib,inpfile,transfile'
  return
endif

; Logfile
if keyword_set(logfile) then logf=logfile else logf=-1
; Catalog format
if n_elements(catformat) eq 0 then catformat='ASCII'

; Testing the files
test = file_test(inpfile)
if test eq 0 then begin
  print,'FILE ',inpfile,' DOES NOT EXIST'
  return
endif

; Using input transformation structure
ninptrans = n_elements(inptrans)
if ninptrans gt 0 then begin

  printlog,logf,'USING INPUT TRANSFORMATION EQUATIONS'
  trans = inptrans
  numbands = n_elements(trans)
  if tag_exist(trans,'night') eq 0 then add_tag,trans,'night',-1,trans
  if tag_exist(trans,'chip') eq 0 then add_tag,trans,'chip',-1,trans
  if tag_exist(trans,'file') eq 0 then add_tag,trans,'file','',trans
  
  printlog,logf,' TRANSFORMATION EQUATIONS'
  printlog,logf,'--------------------------------------------------------------------------------'
  printlog,logf,'  NIGHT/CHIP/FILE  BAND COLOR ZERO-POINT  AIRMASS   COLOR     AIR*COL   COLOR^2 '
  printlog,logf,'--------------------------------------------------------------------------------'
  for i=0,ninptrans-1 do begin
    form1 = '(I10,I6,A6,A7,F10.4,F10.4,F10.4,F10.4,F10.4)'
    form1f = '(A-16,A6,A7,F10.4,F10.4,F10.4,F10.4,F10.4)'
    ; FILENAME
    if trans[i].file ne '' then $
      printlog,logf,format=form1f,trans[i].file,'  '+trans[i].band,trans[i].color,trans[i].zpterm,$
                        trans[i].amterm,trans[i].colterm,trans[i].amcolterm,trans[i].colsqterm
    ; NO Filename
    if trans[i].file eq '' then $
      printlog,logf,format=form1,trans[i].night,trans[i].chip,'  '+trans[i].band,trans[i].color,trans[i].zpterm,$
                        trans[i].amterm,trans[i].colterm,trans[i].amcolterm,trans[i].colsqterm
    form2 = '(A29,F10.4,F10.4,F10.4,F10.4,F10.4)'
    printlog,logf,format=form2,'',trans[i].zptermsig,trans[i].amtermsig,trans[i].coltermsig,$
                      trans[i].amcoltermsig,trans[i].colsqtermsig
  endfor
  printlog,logf,'--------------------------------------------------------------------------------'
  printlog,logf,''
  
; Loading transformation equations from file
endif else begin

  test2 = file_test(transfile)
  if test2 eq 0 and ninptrans eq 0 then begin
    print,'FILE ',transfile,' DOES NOT EXIST'
    return
  endif


  ;# #####################################################
  ;# READ THE TRANSFORMATION FILE
  ;READ_TRANS,transfile,trans,logfile=logf,silent=silent
  READ_TRANS,transfile,trans,logfile=logf,/silent
  numbands = n_elements(trans)

endelse

; Do we have NIGHT information
transnightinfo = 0
if tag_exist(trans,'NIGHT') then begin
  gdtransnight = where(trans.night ge 0,ngdtransnight)
  if ngdtransnight gt 0 then transnightinfo=1
endif
; Do we have CHIP information
transchipinfo = 0
if tag_exist(trans,'CHIP') then begin
  gdtranschip = where(trans.chip ge 0,ngdtranschip)
  if ngdtranschip gt 0 then transchipinfo=1
endif
; Do we have FILE information
transfileinfo = 0
if tag_exist(trans,'FILE') then begin
  gdtransfile = where(trans.file ne '',ngdtransfile)
  if ngdtransfile gt 0 then transfileinfo=1
endif


;###############################
;# READ THE INPUT FILE
;# Read in the input file with metadata information for each
;#   band/exposure in the RAW input photometry files
inparr = importascii(inpfile,/noprint)
ninp = n_elements(inparr)

tags = tag_names(inparr)
ntags = n_elements(tags)

; ---- Transferring to a more user-friendly structure ----
;
;   Old or New format, check if the second value is
;     an integer (new format, NIGHT) or character (old format, band/filter name)
;    All lines need to use the same format (old or new)
newformat = valid_num(inparr[0].(1),/integer)

; --- NEW Format ----
; added NIGHT and CHIP for each band/exposure
; photometry filename, NIGHT, CHIP, filter, airmass, exptime, aperture correction, Band2 ...
if newformat then begin

  numobs = (ntags-1)/6
  dum = {magfile:'',outfile:'',night:lonarr(numobs),chip:lonarr(numobs),band:strarr(numobs),$
         airmass:dblarr(numobs),exptime:dblarr(numobs),apcorr:dblarr(numobs)}
  input = replicate(dum,ninp)
  input.magfile = strtrim(inparr.(0),2)
  for i=0,numobs-1 do begin
    input.night[i] = strtrim(inparr.(1+i*6),2)
    input.chip[i] = strtrim(inparr.(2+i*6),2)
    input.band[i] = strtrim(inparr.(3+i*6),2)
    input.airmass[i] = double(inparr.(4+i*6))
    input.exptime[i] = double(inparr.(5+i*6))
    input.apcorr[i] = double(inparr.(6+i*6))
  endfor

; --- OLD Format ----
; photometry filename, Band1 name, Band1 airmass, Band1 exptime, Band1 aperture correction, Band2 ...
endif else begin

  numobs = (ntags-1)/4
  dum = {magfile:'',outfile:'',night:lonarr(numobs),chip:lonarr(numobs),band:strarr(numobs),$
         airmass:dblarr(numobs),exptime:dblarr(numobs),apcorr:dblarr(numobs)}
  input = replicate(dum,ninp)
  input.magfile = strtrim(inparr.(0),2)
  for i=0,numobs-1 do begin
    input.night[i] = 1   ; dummy value
    input.chip[i] = 1  ; dummy value
    input.band[i] = strtrim(inparr.(1+i*4),2)
    input.airmass[i] = double(inparr.(2+i*4))
    input.exptime[i] = double(inparr.(3+i*4))
    input.apcorr[i] = double(inparr.(4+i*4))
  endfor
endelse


; Making the output filename
ext = 'phot'
for i=0,ninp-1 do begin
  magfile = input[i].magfile
  arr = strsplit(magfile,'.',/extract)
  narr = n_elements(arr)
  outfile = strjoin(arr[0:narr-2],'.')+'.'+ext
  input[i].outfile = outfile
endfor

if not keyword_set(silent) then begin
  printlog,logf,'Running PHOTCALIB on ',strtrim(ninp,2),' input files'
  printlog,logf,''
endif


; Looping through the magnitude files
FOR i=0L,ninp-1 do begin

  inp = input[i]
  magfile = inp.magfile
  magbase = strsplit(magfile,'.',/extract)    ; base name
  if n_elements(magbase) gt 1 then $
    magbase = strjoin(magbase[0:n_elements(magbase)-2],'.')

  ; Testing the file
  test = file_test(magfile)
  if test eq 0 then begin
    printlog,logf,'FILE ',magfile,' DOES NOT EXIST'
    goto,BOMB
  endif

  ; Load the MCH file as well, need this to get the individual
  ; exposure/observation file names
  mchfile = magbase+'.mch'
  if file_test(mchfile) eq 0 then begin
    printlog,logf,'MCH FILE ',mchfile,' DOES NOT EXIST'
    goto,BOMB
  endif
  LOADMCH,mchfile,mfiles
  ; they should be in the same order as the magnitude columns
  ; in the mag file
  
  ; Stars in this file
  if file_isfits(magfile) eq 0 then begin
    numstar = file_lines(magfile)-3L
    ; For 12+ files DAOMASTER starts writing on a second line
    if (numobs ge 12) then numstar = numstar / 2L
    ; Not DAOPHOT file, Header line
    if keyword_set(header) then numstar = file_lines(magfile)-1L
  endif else begin
    hd1 = headfits(magfile,exten=1,/silent)
    numstar = sxpar(hd1,'NAXIS1') 
  endelse
    
  ; Print file info
  if not keyword_set(silent) then begin
    printlog,logf,format='(A-9,A-20)','FILE ',input[i].magfile
    if newformat then begin
      printlog,logf,format='(A-9,'+strtrim(numobs,2)+'I-7)','NIGHT',input[i].night
      printlog,logf,format='(A-9,'+strtrim(numobs,2)+'I-7)','CHIP',input[i].chip
    endif
    printlog,logf,format='(A-9,'+strtrim(numobs,2)+'A-7)','BAND',input[i].band
    printlog,logf,format='(A-9,'+strtrim(numobs,2)+'F-7.4)','AIRMASS',input[i].airmass
    printlog,logf,format='(A-9,'+strtrim(numobs,2)+'F-7.1)','EXPTIME',input[i].exptime
    printlog,logf,format='(A-9,'+strtrim(numobs,2)+'F-7.4)','APCORR',input[i].apcorr
    printlog,logf,format='(A-9,I-8)','NSTARS',numstar
    printlog,logf,''
  endif


  ;#############################
  ;# READING IN THE PHOTOMETRY
  ;#############################
  phot = PHOTRED_READFILE(magfile)
  numstar = n_elements(phot)
  tags = tag_names(phot)
  ncol = n_tags(phot)
  nextra = ncol - 2*numobs
  mastable = dblarr(numstar,2*numobs+nextra)
  for j=0,ncol-1 do mastable[*,j] = phot.(j)

  raarray = reform(mastable[*,1])
  decarray = reform(mastable[*,2])

  ; Initializing the arrays
  goodstar = dblarr(2*numbands+nextra)
  manystars = dblarr(2*numobs)
  startable = dblarr(numstar,2*numbands+nextra)
  indystar = dblarr(numstar,2*numobs)


  ;########################
  ;# Making the transformation structure for the bands of this input file
  mastrans = replicate(trans[0],numobs)

  ; Associate transformation equations with each observation
  for j=0,numobs-1 do begin

    ; Get the filename of each observation, get from mch file
    obsfile = file_basename(mfiles[j],'.als')
     
    ; There are four options for matching the file:
    ; 1) No chip/night/filename information
    ; 2) Only chip given
    ; 3) Night+chip given
    ; 4) Filename given
    ; This information can be different for each file, i.e. multiple
    ; formats can be used in a given trans file.  Try MOST specific
    ; (#4) to LEAST specific (#1).

    nmatch = 0
    ; Try filename + band
    if transfileinfo eq 1 then $
       MATCH,trans.file+':'+trans.band,obsfile+':'+inp.band[j],ind1,ind2,/sort,count=nmatch
    ; Try night+chip + band
    ;  only want to match lines with file=''
    if nmatch eq 0 and transchipinfo eq 1 and transnightinfo eq 1 then $
       MATCH,trans.file+':'+strtrim(trans.night,2)+':'+strtrim(trans.chip,2)+':'+trans.band,$
             ':'+strtrim(inp.night[j],2)+':'+strtrim(inp.chip[j],2)+':'+inp.band[j],ind1,ind2,/sort,count=nmatch
    ; Try chip + band
    ;  only want to match lines with file='' and night=-1
    if nmatch eq 0 and transchipinfo eq 1 then $
       MATCH,trans.file+':'+strtrim(trans.night,2)+':'+strtrim(trans.chip,2)+':'+trans.band,$
             ':-1:'+strtrim(inp.chip[j],2)+':'+inp.band[j],ind1,ind2,/sort,count=nmatch
    ; Try just the band
    ;  only want to match lines with file='', night=-1 and chip=-1
    if nmatch eq 0 then $
       MATCH,trans.file+':'+strtrim(trans.night,2)+':'+strtrim(trans.chip,2)+':'+trans.band,$
             ':-1:-1:'+inp.band[j],ind1,ind2,/sort,count=nmatch

    ; Found the transformation for this file+band
    if (nmatch gt 0) then begin
      mastrans[j] = trans[ind1[0]]
    endif else begin
      printlog,logf,'NO TRANSFORMATION INPUT FOR  REFILE=',magbase,' OBSFILE=',obsfile,' NIGHT=',strtrim(inp.night[j],2),' CHIP=',$
                strtrim(inp.chip[j],2),' FILTER=',inp.band[j]
      return
    endelse

    ; Check that the color exists
    gdcol = where(inp.band eq mastrans[j].colband,ngdcol)
    if (ngdcol eq 0) then begin
      printlog,logf,mastrans[j].colband,' BAND NOT FOUND. CANNOT FORM ',mastrans[j].color,' COLOR FOR BAND ',inp.band[j]
      return
    endif
  endfor


  ;##################
  ;# CALIBRATING
  ;##################
  ; Ah, the meat of the program.  One by one, each star is popped off of
  ; the observation file, thrown into the solution engine (SOLVESTAR), solved
  ; and then brought out as an average solution (goodstar) and individual measures
  ; (indystar).  After that, frame to frame residuals are calculated
  ; indystar is where the individual solved magnitudes will be stored
  SOLVESTAR,mastable,mastrans,inp,goodstar


  ;########################
  ;# PREPARING THE OUTPUT
  ;########################

  ; Getting the unique passbands
  ui = uniq(mastrans.band,sort(mastrans.band))
  ui = ui[sort(ui)]
  ubands = mastrans[ui].band
  nubands = n_elements(ubands)

  ;--------------
  ; Head columns 
  ;--------------
  finalstar = goodstar[*,0:2]
  headline = '     ID       X         Y     '
  ;;format = '(2X,I5,2F9.3'
  ;format = '(2X,I5,2F10.3'
  format = '(2X,I6,2F10.3'


  ;--------------------------------------------
  ; Instrumental Individual Magnitudes columns
  ;--------------------------------------------
  if keyword_set(keepinstrumental) then begin

    ; mastable is: id, x, y, unsolved magnitudes, chi, sharp
    instrstar = mastable[*,3:numobs*2+2]
    
    ; Add to the final array
    finalstar = [ [finalstar], [instrstar] ]
    
    ; Looping through the unique bands
    instroutband = strarr(numobs)
    instrouterr = strarr(numobs)
    for j=0,nubands-1 do begin
      gdbands = where(mastrans.band eq ubands[j],ngdbands)      
      ; More than one observation in this band
      if (ngdbands gt 1) then begin
        instroutband[gdbands] = 'I_'+ubands[j]+strtrim(indgen(ngdbands)+1,2)
        instrouterr[gdbands] = 'I_'+ubands[j]+strtrim(indgen(ngdbands)+1,2)+'ERR'

      ; Only ONE obs in this band
      endif else begin
        instroutband[gdbands[0]] = 'I_'+ubands[j]
        instrouterr[gdbands[0]] = 'I_'+ubands[j]+'ERR'
      endelse
      
    end  ; looping through unique bands

    ; header
    instrbandspace = strarr(numobs)
    instrerrspace = strarr(numobs)
    for j=0,numobs-1 do instrbandspace[j] = string(bytarr( (10-strlen(instroutband[j])) > 1)+32B)
    for j=0,numobs-1 do instrerrspace[j] = string(bytarr( (8-strlen(instrouterr[j])) > 1)+32B)
    for j=0,numobs-1 do headline=headline+instroutband[j]+instrbandspace[j]+instrouterr[j]+instrerrspace[j]

    ; format
    format = format+','+strtrim(2*numobs,2)+'F9.4'
  endif


  ;------------------------------------------
  ; Calibrated Individual Magnitudes columns
  ;------------------------------------------
  if not keyword_set(onlyaverage) then begin
    indivstar = goodstar[*,3:(numobs*2)+2]

    ; Add to the final array
    finalstar = [ [finalstar], [indivstar] ]

    ; Looping through the unique bands
    outband = strarr(numobs)
    outerr = strarr(numobs)
    for j=0,nubands-1 do begin
      gdbands = where(mastrans.band eq ubands[j],ngdbands)
      
      ; More than one observation in this band
      if (ngdbands gt 1) then begin
        outband[gdbands] = ubands[j]+'MAG'+strtrim(indgen(ngdbands)+1,2)
        outerr[gdbands] = ubands[j]+strtrim(indgen(ngdbands)+1,2)+'ERR'

      ; Only ONE obs in this band
      endif else begin
        outband[gdbands[0]] = ubands[j]+'MAG'
        outerr[gdbands[0]] = ubands[j]+'ERR'
      endelse
    endfor  ; looping through unique bands

    ; header
    bandspace = strarr(numobs)
    errspace = strarr(numobs)
    for j=0,numobs-1 do bandspace[j] = string(bytarr( (10-strlen(outband[j])) > 1)+32B)
    for j=0,numobs-1 do errspace[j] = string(bytarr( (8-strlen(outerr[j])) > 1)+32B)
    for j=0,numobs-1 do headline=headline+outband[j]+bandspace[j]+outerr[j]+errspace[j]

    ; format
    format = format+','+strtrim(2*numobs,2)+'F9.4'
  endif


  ;---------------------------------------
  ; Calibrated Average Magnitudes columns
  ;---------------------------------------
  if keyword_set(average) or keyword_set(combine) or keyword_set(onlyaverage) then begin

    ; Looping through the unique bands
    undefine,multibands,multierr
    for j=0,nubands-1 do begin
      gdbands = where(mastrans.band eq ubands[j],ngdbands)
      ; More than one observation in this band
      if (ngdbands gt 1) then begin
        PUSH,multibands,ubands[j]+'MAG'
        ;PUSH,multibands,ubands[j]
        PUSH,multierr,ubands[j]+'ERR'

      ; Only ONE obs in this band
      endif else begin
        if keyword_set(onlyaverage) then begin
          PUSH,multibands,ubands[j]+'MAG'
          PUSH,multierr,ubands[j]+'ERR'
        endif
      endelse
    endfor  ; looping through unique bands
    nmultibands = n_elements(multibands)

    ; We have some multi observations
    if (nmultibands gt 0) then begin

      ; Combining the data
      combstar = fltarr(numstar,nmultibands*2)

      ; Looping through the unique bands
      for j=0,nmultibands-1 do begin
        gdband = where(mastrans.band+'MAG' eq multibands[j],ngdband)
        ; Multiple exposures in this band
        if (ngdband gt 1) then begin
          ; Average the mags
          AVERAGEMAG,goodstar[*,3+2*gdband],goodstar[*,4+2*gdband],newmag,newerr,/robust
          ; Now put in the final output array
          combstar[*,j*2] = newmag
          combstar[*,j*2+1] = newerr

        ; One exposure in this band
        endif else begin
          combstar[*,j*2] = goodstar[*,3+2*gdband[0]]    ; transfer mag
          combstar[*,j*2+1] = goodstar[*,4+2*gdband[0]]    ; transfer error
        endelse
      endfor ; multiband loop

      ; Add to the final array
      finalstar = [ [finalstar], [combstar] ]

      ; header
      multibandspace = strarr(nmultibands)
      multierrspace = strarr(nmultibands)
      for j=0,nmultibands-1 do multibandspace[j] = string(bytarr( (10-strlen(multibands[j])) > 1)+32B)
      for j=0,nmultibands-1 do multierrspace[j] = string(bytarr( (8-strlen(multierr[j])) > 1)+32B)
      for j=0,nmultibands-1 do headline=headline+multibands[j]+multibandspace[j]+multierr[j]+multierrspace[j]

      ; format
      format = format+','+strtrim(2*nmultibands,2)+'F9.4'

    endif   ; nmultibands gt 0
  endif  ; combine, average or onlyaverage
  

  ;---------------
  ; Extra columns
  ;---------------
  extrastar = goodstar[*,numobs*2+3:numobs*2+nextra-1]

  ; Add to the final array
  finalstar = [ [finalstar], [extrastar] ]

  ; We have header from input
  if keyword_set(header) then begin

    extratags = tags[2*numobs+3:*]
    headline = headline + strjoin(extratags,'   ')

    ; format
    for j=2*numobs+3,ncol-1 do begin
      form = 'F9.4'
      type = SIZE(phot[0].(j),/type)
      if type eq 1 then form='I5'  ; byte
      if type eq 2 then form='I9'  ; int
      if type eq 3 then form='I12'  ; long
      if type eq 4 then form='F11.4'  ; float
      if type eq 5 then form='F13.6'  ; double
      if type eq 7 then form='A20'  ; string

      if tags[j] eq 'CHI' and type eq 4 then form='F9.4'
      if tags[j] eq 'SHARP' and type eq 4 then form='F9.4'
      if tags[j] eq 'FLAG' and type eq 3 then form='I5'
      if tags[j] eq 'PROB' and type eq 4 then form='F7.2'

      format = format+','+form
    endfor

    format = format+')'

  ; Figure out the names for the Extra columns
  endif else begin

    ; header
    headline = headline+' CHI      SHARP'
    ; We have FLAG/PROB colums
    if (nextra eq 7) then begin
      headline = headline+'  FLAG  PROB'
    endif
    ; We have other columns
    if nextra gt 5 and nextra ne 7 then begin
      for j=0,nextra-6 do headline = headline+'    EXTRA'+strtrim(j+1,2)
    endif

    ; format
    format = format+',2F9.4'     ; chi and sharp
    ; We have FLAG/PROB colums   
    if nextra eq 7 then begin
      format = format+',I5,F7.2'
    endif
    ; We have other columns
    if nextra gt 5 and nextra ne 7 then begin
      for j=0,nextra-6 do format = format+','+strtrim(nextra-5,2)+'F9.4'
    endif
    format = format+')'
  endelse



  ;----------------
  ; WRITE the file
  ;----------------

  ; Printing results to output file
  ; Header information is printed (an index of columns) and then
  ; the stars, one by one
  outfile = input[i].outfile

  ;; ASCII file
  if catformat eq 'ASCII' then begin
stop
    OPENW,unit,/get_lun,outfile  
    ; Print the header
    printf,unit,headline
    ; Loop through the stars
    for d=0l,numstar-1 do printf,unit,format=format,finalstar[d,*]
    CLOSE,unit
    FREE_LUN,unit

  ;; FITS file
  endif else begin
    ; Get column types
    ncolumns = n_elements(finalstar[0,*])
    formatarr1 = strsplit(strmid(format,1,strlen(format)-2),',',/extract)
    bdf = where(stregex(formatarr1,'X',/boolean,/fold_case) eq 1,nbdf)
    if nbdf gt 0 then REMOVE,bdf,formatarr1
    ;; Deal with repeats, eg.g 5F9.4
    undefine,formatarr
    for j=0,n_elements(formatarr1)-1 do begin
      if valid_num(strmid(formatarr1[j],0,1)) eq 1 then begin
        vnum = intarr(strlen(formatarr1[j]))
        fbytes = byte(formatarr1[j])
        for k=0,strlen(formatarr1[j])-1 do vnum[k]=valid_num(fbytes[k])
        hi = first_el(where(vnum eq 0))
        fnum = strmid(formatarr1[j],0,hi)
        fmt1 = strmid(formatarr1[j],hi)
        push,formatarr,replicate(fmt1,fnum)
      endif else push,formatarr,formatarr1[j]
    endfor
    ; Convert format codes to column types
    types = lonarr(ncolumns)
    for j=0,ncolumns-1 do begin
       let = strmid(formatarr[j],0,1)
       dum = strsplit(strmid(formatarr[j],1),'.')
       ndig = dum[0]
       if n_elements(dum) gt 1 then ndec=dum[1] else ndec=0
       case let of
        'A': types[j]=7                                    ; string
        'I': begin
              if ndig lt 9 then types[j]=1                 ; byte
              if ndig ge 9 and ndig lt 12 then types[j]=2  ; int
              if ndig ge 12 then types[j]=3                ; long
           end
        'F': begin
              types[j]=4                                   ; float
              if ndig ge 13 or ndec ge 6 then types[j]=5   ; double
           end
        else:
      endcase
    endfor
    ;; Get the column names
    colnames = strsplit(headline,' ',/extract)
    if colnames[0] eq '#' then colnames=colnames[1:*]
    ;; Create the schema
    schema = create_struct(colnames[0],fix('',types[0]))
    for j=1,ncolumns-1 do schema=create_struct(schema,colnames[j],fix('',types[j]))
    ;; Copy in the data
    final = create_struct(schema,numstar)
    for j=0,ncolumns-1 do final.(j)=finalstar[*,j]
stop
    ;; Write the file
    MWRFITS,final,outfile,/create
  endelse ; FITS

  printlog,logf,'Final Photometry File is = ',outfile
  
  BOMB:

ENDFOR ; loop through the magnitude files

if not keyword_set(silent) then printlog,logf,'PHOTCALIB FINISHED'

if keyword_set(stp) then stop

end
