;+
;
; APCORRECT
;
; This applies the DAOGROW aperture correction to DAOPHOT
; aperture photometry.  The output is almost identical to
; the TOT files that DAOGROW outputs.  The final error is
; slightly different by ~0.0002 mags, not sure why, maybe
; from round-off errors in the DAOGROW growth-curves.
;
; The DAOGROW growth-curve, cumulative aperture corrections,
; and the cumulative aperture correction errors need to be
; input (or the "grofile" and "gronum" input).  The "adopt",
; "cum", and "cumsig" arrays can be gotten from the .GRO file
; that DAOGROW outputs.  DAOGROW should be run on ALL of the
; frames for a given night.  For each frame the .GRO output
; looks like this:
;
; --> 0001 Dn1-ccd544a.ap       ...           <sky> =      5.2     Ro =     1.4580381
;
; --> 0001             Model  -0.1222  -0.0613  -0.0312  -0.0190  -0.0126  -0.0093 ...
; --> 0001
; --> 0001          Observed  -0.1188  -0.0683  -0.0373  -0.0215  -0.0131  -0.0092 ...
; --> 0001          Sigma(O)   0.0012   0.0012   0.0012   0.0012   0.0013   0.0015 ...
; --> 0001
; --> 0001           Adopted  -0.1188  -0.0680  -0.0365  -0.0208  -0.0128  -0.0092 ...
; --> 0001          Sigma(A)   0.0119   0.0069   0.0038   0.0023   0.0016   0.0012 ...
; --> 0001
; --> 0001    Cumulative  -0.3142  -0.1954  -0.1274  -0.0909  -0.0701  -0.0573 ...
; --> 0001      Sigma(C)   0.0148   0.0087   0.0053   0.0037   0.0029   0.0024 ...
;
;
;  2398  791.80 1919.51:  13.067   13.067   13.068   13.066   13.065   13.064 ...
;                          0.0152   0.0094   0.0064   0.0050   0.0046   0.0044 ...
;  1818 1010.07 1523.88:  13.150   13.145   13.143   13.142   13.141   13.141 ...
;                          0.0152   0.0098   0.0069   0.0057   0.0052   0.0051 ...
;  2183 1960.59 1774.89:  13.242   13.244   13.244   13.243   13.243   13.242 ...
;                          0.0153   0.0095   0.0065   0.0053   0.0048   0.0046 ...
;
; It goes on to output the "good" magnitudes and errors for all of the
; stars for that frame.  If there are Nap number of apertures then there
; will be Nap-1 elements in the growth curves.  There will be Nap
; cumulative aperture corrections (one for each aperture) but only
; sigma errors for the first Nap-1 corrections.
;
; INPUTS:
;  aper      The aperture photometry structure loaded by
;              LOADAPER.PRO.
;  adopt     The adopted curve-of-growth array that DAOGROW found
;              for *THIS* frame.  This changes from frame to
;              frame, so be careful.
;  cum       The adopted cumulative aperture corrections array for
;              this frame.
;  cumsig    The errors for "cum".  Normally this is one element
;              less than "cum".
;  =grofile  The filename of the .GRO to use for the aperture
;              correction.  "apername" must also be input.
;  =apername The name of the aperture file, e.g. "gn10-00278875_50a.ap".
;  /silent   Don't print anything to the screen.
;  /stp      Stop at the end of the program.
;
; OUTPUTS:
;  final     The structure with the final cumulative photometry.
;              Same output as the .TOT files that DAOGROW
;              creates.
;  =error    Error message if any errors occured.
;
; USAGE:
;  IDL>apcorrect,aper,final,adopt,cum,cumsig,grofile=grofile,apername=apername,
;                silent=silent,stp=stp
;
; By D.Nidever  May 2008
;-

pro apcorrect,aper0,final,adopt,cum,cumsig,grofile=grofile,$
              apername=apername,silent=silent,stp=stp,error=error

undefine,error,final

naper = n_elements(aper0)
nadopt = n_elements(adopt)
ncum = n_elements(cum)
ncumsig = n_elements(cumsig)
ngrofile = n_elements(grofile)
napername = n_elements(apername)

; Not enough inputs
if NOT (naper gt 0 and nadopt gt 0 and ncum gt 0 and ncumsig gt 0) and $
   NOT (naper gt 0 and ngrofile gt 0 and napername gt 0) then begin
  print,'Syntax - apcorrect,aper,final,adopt,cum,cumsig,grofile=grofile,apername=apername,'
  print,'                   silent=silent,stp=stp'
  error = 'Not enough inputs'
  return
endif


; Grab the cumulative aperture correction and the errors
;-------------------------------------------------------
if (ngrofile gt 0 and napername gt 0) then begin

  ; Load the .GRO file
  grotest = FILE_TEST(grofile[0])
  if (grotest eq 0) then begin
    print,grofile[0],' NOT FOUND'
    error = grofile[0]+' NOT FOUND'
    return
  endif
  READLINE,grofile,grolines,count=ngrolines
  if not keyword_set(silent) then print,'DAOGROW .GRO FILE = ',grofile[0]

  if (ngrolines eq 0) then begin
    print,grofile[0],' IS EMPTY'
    error = grofile[0]+' IS EMPTY'
    return
  endif

  ; Frame number in .GRO file
  if VALID_NUM(apername) eq 1 then begin
    print,'APERNAME='+apername+' MUST BE A STRING'
    error = 'APERNAME='+apername+' MUST BE A STRING'
    return
  endif

  ; The frame number is not reliable, it can change with respect to
  ; the INF list if there are "bad" frames, i.e. no good stars.
  ; Use the frame name to identify the right frame number
  ;   --> 4254 gn10-00278875_50a.ap 
  frameind = where(stregex(grolines,'--> ',/boolean) eq 1 and $
                   stregex(grolines,apername,/boolean) eq 1,nframeind)
  if nframeind eq 0 then begin
    error = apername+' NOT FOUND in '+grofile
    print,error
    return
  endif
  gronum = long((strsplit(grolines[frameind[0]],' ',/extract))[1])  ; pull out the GRO frame number
  framegrolines = grolines[frameind[0]:(frameind[0]+12)<(ngrolines-1)]

  if not keyword_set(silent) then print,'Frame number = ',strtrim(gronum,2),'   ',apername
  num = string(gronum,format='(I04)')

  ; --> 0001           Adopted  -0.1188  -0.0680  -0.0365  -0.0208  -0.0128  -0.0092  -0.0077 ...
  ; --> 0001          Sigma(A)   0.0119   0.0069   0.0038   0.0023   0.0016   0.0012   0.0011 ...
  ; --> 0001
  ; --> 0001    Cumulative  -0.3142  -0.1954  -0.1274  -0.0909  -0.0701  -0.0573  -0.0480 ...
  ; --> 0001      Sigma(C)   0.0148   0.0087   0.0053   0.0037   0.0029   0.0024   0.0021 ...
  ; "Sigma(C)" has one less element than "Cumulative"

  ; Getting adopted growth curve
  adoptind = where(stregex(framegrolines,'--> '+num,/boolean) eq 1 and $
                   stregex(framegrolines,'Adopted',/boolean) eq 1,nadoptind)
  if (nadoptind eq 0) then begin
    print,'Adopted growth-curve for frame '+gronum+' NOT FOUND'
    error = 'Adopted growth-curve for frame '+gronum+' NOT FOUND'
    return
  endif
  stradopt = framegrolines[adoptind[0]]
  stradoptarr = strsplit(stradopt,' ',/extract)
  adopt = float(stradoptarr[3:*])
  nadopt = n_elements(adopt)
  if not keyword_set(silent) then $
    print,'Adopted ',adopt,format='(A-10,'+strtrim(nadopt,2)+'F9.4)'

  ; Getting cumulative aperture correction
  cumind = where(stregex(framegrolines,'--> '+num,/boolean) eq 1 and $
                 stregex(framegrolines,'Cumulative',/boolean) eq 1,ncumind)
  if (ncumind eq 0) then begin
    print,'Cumulative aperture corrections for frame '+gronum+' NOT FOUND'
    error = 'Cumulative aperture corrections for frame '+gronum+' NOT FOUND'
    return
  endif
  strcum = framegrolines[cumind[0]]
  strcumarr = strsplit(strcum,' ',/extract)
  cum = float(strcumarr[3:*])
  ncum = n_elements(cum)
  if not keyword_set(silent) then $
    print,'Cumulative ',cum,format='(A-10,'+strtrim(ncum,2)+'F9.4)'

  ; Getting cumulative aperture correction errors
  csigind = where(stregex(framegrolines,'--> '+num,/boolean) eq 1 and $
                  stregex(framegrolines,'Sigma\(C\)',/boolean) eq 1,ncsigind)
  if ncsigind eq 0 then begin  ;; newer daogrow use Sigma_C
    csigind = where(stregex(framegrolines,'--> '+num,/boolean) eq 1 and $
                    stregex(framegrolines,'Sigma_C',/boolean) eq 1,ncsigind)
  endif
  if (ncsigind eq 0) then begin
    print,'Cumulative aperture correction errors for frame '+gronum+' NOT FOUND'
    error = 'Cumulative aperture correction errors for frame '+gronum+' NOT FOUND'
    return
  endif
  strcsig = framegrolines[csigind[0]]
  strcsigarr = strsplit(strcsig,' ',/extract)
  cumsig = float(strcsigarr[3:*])
  ncumsig = n_elements(cumsig)
  if not keyword_set(silent) then $
    print,'Sigma(C) ',cumsig,format='(A-10,'+strtrim(ncumsig,2)+'F9.4)'

endif  ; getting aperture correction



; How many stars were input
ntotstars = n_elements(aper0)
if not keyword_set(silent) then print,'Nstars = ',strtrim(ntotstars,2)

; How many apertures
napertures = n_elements(aper0[0].mag)
if not keyword_set(silent) then print,strtrim(napertures,2),' apertures'

; Get only stars with "good" photometry
aper = aper0
;gdphot = where(aper.err[0] lt 9.99,ngdphot)
gdphot = where(aper.err[0] lt 9.0,ngdphot)
; No stars with good photometry
if (ngdphot eq 0) then begin
  if not keyword_set(silent) then print,'NO stars with good photometry'
  dum = {ID:0L,x:0.0,y:0.0,mag:99.999,err:9.9999,sky:0.0,magfap:99.999,apcorr:99.999,finalap:-1L}
  final = REPLICATE(dum,ntotstars)
  final.id = aper0.id
  final.x = aper0.x
  final.y = aper0.y
  final.sky = aper0.sky
  return
endif
aper = aper[gdphot]
nstars = n_elements(aper)
if not keyword_set(silent) then print,strtrim(nstars,2),' stars have good photometry'




; Compute "curve-of-growth" for ever star
; compare to "adopted" curve of growth
; and compute weight and SIGSQ
; also get best aperture and its error
;dmag = fltarr(nstars,napertures)+9.9999
;diff = fltarr(nstars,napertures)+9.9999
;weight = fltarr(nstars,napertures)
;sigsq = fltarr(nstars,napertures)+999999.

bestap = fltarr(nstars)
finalmag = fltarr(nstars)   
sfinal = fltarr(nstars)   
apcorr = fltarr(nstars)   
totmag = fltarr(nstars)   
     
; Looping through all stars
for k=0,nstars-1 do begin

  gd = where(aper[k].err lt 9.0,ngd)
  ;gd = gd[0:(ngd-1)<(napertures-1)]   ; max napertures 
  nap = ngd < (napertures-1)

  sigmin = 1e38

  mag = aper[k].mag[0:nap-1]
  sig = aper[k].err[0:nap-1]

  ; ADOPT, WCUM and CUM are shifted by one compared to the DAOGROW code
  ; WCUM[J] instead of WCUM(J+1).
  ; ADOPT[J-1] instead of ADOPT(J)

  wcum = cumsig^2.
  w = fltarr(nap)
  obs = fltarr(nap)
  wobs = fltarr(nap)
  sigsq = fltarr(nap)
  For j=0,nap-1 do begin

    if (j ne 0) then begin
      diff = mag[j] - mag[j-1] - adopt[j-1]
      W[J] = 1.0/( 1.0 + ( DIFF/(2.*SIG[J]) )^2.0 )
      W[J] = MIN( [ W[J], W[J-1] ] )
    endif else begin
      W[J] = 1.0
    endelse
    ;SIGSQ = WCUM[J+1] + SIG[J]^2.0/W[J]
    ;OBS[J] = MAG[J]+CUM[J+1]
    SIGSQ[J] = WCUM[J] + (SIG[J]^2.0)/W[J]
    OBS[J] = MAG[J]+CUM[J]
    WOBS[J] = MIN( [9.999, SQRT(SIGSQ[J])] )
    if (sigsq[j] lt sigmin) then begin
      FINAL = MAG[J]
      SFINAL1 = WOBS[J]
      JFINAL = J
      SIGMIN = SIGSQ[J]
    endif

  Endfor  ; aperture loop

  bestap[k] = JFINAL+1
  finalmag[k] = FINAL
  sfinal[k] = SFINAL1
  apcorr[k] = cum[JFINAL]
  totmag[k] = finalmag[k] + apcorr[k]

endfor

;  ; More than 1 good ap
;  if (nap gt 1) then begin
;
;    mag = aper[k].mag[0:nap-1]
;    err = aper[k].err[0:nap-1]
;    err2 = err[1:*]             ; want the error for the outer ap
;
;    ; We need to calculate how different the star's
;    ; growth-curve is from the adopted frame growth-curve
;    ; This is NOT done for the 1st aperture, only the 2nd and up.
;
;    ; Calculate the star's empirical growth-curve
;    magout = mag[1:*]
;    magin = mag[0:nap-2]
;    dmag = magout-magin
;
;    ; Calculate the difference between the star's
;    ; empirical growth-curve and the adopted growth-curve
;    ; for this frame
;    diff = dmag - adopt[0:nap-2]
;
;    ; Calculate the weight
;    wt1 = 1.0/( 1.0 + ( diff/(2.0*err2) )^2.0 )   ; error for outer ap
;
;    ; Get minimum of current and previous aperture's weight
;    if (nap gt 2) then begin
;      wt2 = [1.0, wt1[0:nap-3]]           ; the weight of the previous ap
;      wt = MIN( [[wt1],[wt2]], dim=2)
;    endif else wt=MIN([1.0,wt1])
;
;    ; Add weight for the first aperture
;    wt = [1.0,wt]
;
;    ; Calculate the total error for each aperture
;    sigsq = cumsig[0:nap-1]^2.0 + (err^2.0)/wt
;
;
;  ; Only one good aperture
;  endif else begin
;
;    mag = aper[0].mag[0]
;    err = aper[k].err[0]
;    wt = 1.0
;
;    ; SIGSQ = WCUM(J+1) + SIG(J,ISTAR)**2/W(J)
;    sigsq = cumsig[0]^2.0 + (err^2.0)/wt
;
;  endelse
;
;  ; Find aperture with lowest sigsq
;  bestind = first_el(minloc(sigsq))
;
;  bestap[k] = bestind+1
;  finalmag[k] = mag[bestind]
;  sfinal[k] = MIN([9.999,SQRT(sigsq[bestind])])
;  apcorr[k] = cum[bestind]
;  totmag[k] = finalmag[k] + apcorr[k]
;
; end


; Make the final structure
dum = {ID:0L,x:0.0,y:0.0,mag:99.999,err:9.9999,sky:0.0,magfap:99.999,apcorr:99.999,finalap:-1L}
final = REPLICATE(dum,ntotstars)
final.id = aper0.id
final.x = aper0.x
final.y = aper0.y
final[gdphot].mag = totmag
final[gdphot].err = sfinal
final.sky = aper0.sky
final[gdphot].magfap = finalmag
final[gdphot].apcorr = apcorr
final[gdphot].finalap = bestap



if keyword_set(stp) then stop

end



;C At this point, the arrays MAG and SIG contain the raw magnitudes and
;C standard errors of the aperture photometry for all the stars; the 
;C array IFILE tells which of the input files each star belongs to.
;C
;C Now we will compute the magnitude differences:  MAG(1) will still be
;C the magnitude in the first aperture, but now for each star MAG(2)
;C will contain aperture 2 - aperture 1, MAG(3) will contain aperture 3 -
;C aperture 2, and so on.  For simplicity, the raw standard error of the
;C larger aperture in each case will be taken to be the standard error
;C of the difference.
;C
;      DO I=1,NSTAR
;         IF (NAP(I) .GT. 1) THEN
;            DO J=NAP(I),2,-1
;               MAG(J,I) = MAG(J,I)-MAG(J-1,I)
;            END DO
;         END IF
;      END DO
;.....

;DLN: ADOPT is the adopted curve of growth for this frame
;DLN: DIFF must be the difference between this star's curve of growth
;DLN:   and the adopted curve of growth for this frame
;DLN:   It probably accounts for any deviations from the growth curve
;DLN:    due to other stars in the aperture.
;
;            DO J=1,NAP(ISTAR)
;               IF (J .NE. 1) THEN
;                  DIFF = MAG(J,ISTAR) - ADOPT(J)
;                  MAG(J,ISTAR) = MAG(J-1,ISTAR) + MAG(J,ISTAR)
;C
;C MAG contains the magnitudes again.
;C
;                  W(J) = 1./( 1. + ( DIFF/(2.*SIG(J,ISTAR)) )**2 )
;                  W(J) = MIN( W(J), W(J-1) )
;               ELSE
;                  W(J) = 1.
;               END IF
;               SIGSQ = WCUM(J+1) + SIG(J,ISTAR)**2/W(J)
;               OBS(J) = MAG(J,ISTAR)+CUM(J+1)
;               WOBS(J) = MIN(9.999, SQRT(SIGSQ))
;               IF (SIGSQ .LT. SIGMIN) THEN
;                  FINAL = MAG(J,ISTAR)
;                  SFINAL = WOBS(J)
;                  JFINAL = J
;                  SIGMIN = SIGSQ
;               END IF
;            END DO
;            S1 = RNDOFF (X(ISTAR), 9, 3)
;            S2 = RNDOFF (Y(ISTAR), 9, 3)
;            S3 = RNDOFF (SKY(ISTAR), 9, 3)
;            WRITE (2,222) ID(ISTAR), S1, S2, FINAL+CUM(JFINAL+1), 
;     .           SFINAL, S3, FINAL, CUM(JFINAL+1), JFINAL
