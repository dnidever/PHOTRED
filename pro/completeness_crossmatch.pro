;+
;
; COMPLETENESS_CROSSMATCH
;
; This performs the cross-match between the recovered and original
; and artificial stars for completeness.pro.
;
; INPUTS:
;  final     The final recovered photometry structure.
;  orig      The original photometry structure.
;  synth     The injected artificial stars information structure.
;  =logfile  The logfile name.
;
; OUTPUTS:
;  find      The index for recovered artificial stars in "final".
;  sind      The index for recovered artificial stars in "synth".
;  =nmatch   The number of artificial stars recovered.
;  =error    The error message if one occurred.
;
; USAGE:
;  IDL>completeness_crossmatch,final,orig,synth,find,sind,logfile=logfile
;
; By D. Nidever  August 2017
;-

pro completeness_crossmatch,final,orig,synth,find,sind,nmatch=nmatch,logfile=logfile,error=error

undefine,find,sind
undefine,error
nmatch = 0

; Not enough inputs
if n_elements(final) eq 0 or n_elements(orig) eq 0 or n_elements(synth) eq 0 then begin
  print,'Syntax - completeness_crossmatch,final,orig,synth,find,sind,nmatch=nmatch,logfile=logfile,error=error'
  error = 'Not enough inputs'
  return 
endif
; Logfile
if keyword_set(logfile) then logf=logfile else logf=-1

; Do the CROSS-MATCHING between the three lists (final/orig/synth)
;-----------------------------------------------------------------
; If two stars land right on top of each other and only one is recovered then use the magnitude
; to match it up to the correct star (real or artificial)
; We don't want to think we recovered an artificial star when it's
; actually a real star (likely bright one).
; We are using the instrumental photomery (.mag/.raw) to do this matching.

; Crossmatch final to orig and final to synth
SRCMATCH,final.x,final.y,orig.x,orig.y,2,oind1,oind2,count=nomatch
SRCMATCH,final.x,final.y,synth.x,synth.y,2,aind1,aind2,count=nastmatch

; Deal with duplicates
;  all we care about is real objects falsely
;  identified as ASTs, so bad matches in AIND
finaldblind = doubles([oind1,aind1],count=nfinaldbl)
if nfinaldbl gt 0 then begin
  printlog,logfile,'  Resolving '+strtrim(nfinaldbl,2)+' duplicates'
  finaldbl = ([oind1,aind1])(finaldblind)
  flag = lonarr(nfinaldbl)  ; 0-real, 1-ast
  origdbl = lonarr(nfinaldbl)
  synthdbl = lonarr(nfinaldbl)
  ; Get magnitude column indices for the three structures
  finaltags = tag_names(final)
  finalmagind = where(stregex(finaltags,'^MAG',/boolean) eq 1,nfinalmagind)
  origtags = tag_names(orig)
  origmagind = where(stregex(origtags,'^MAG',/boolean) eq 1,norigmagind)
  synthtags = tag_names(synth)
  synthmagind = where(stregex(synthtags,'^MAG',/boolean) eq 1,nsynthmagind)
  if nfinalmagind ne norigmagind or nfinalmagind ne nsynthmagind then begin
    printlog,logfile,'  FINAL/ORIG/SYNTH photometry files have different magnitude columns'
    error = 'FINAL/ORIG/SYNTH photometry files have different magnitude columns'
    nmatch = 0
    return
  endif
  ; Loop over duplicates
  for k=0,nfinaldbl-1 do begin
    finaldbl1 = finaldbl[k]  ; the FINAL index
    final1 = final[finaldbl1]
    finalmag = fltarr(nfinalmagind)
    for l=0,nfinalmagind-1 do finalmag[l]=final1.(finalmagind[l])

    ; ORIG
    MATCH,oind1,finaldbl1,ind1,/sort
    origdbl[k] = ind1  ; index into OIND1/2
    orig1 = orig[oind2[ind1]]
    odist = sqrt( (final1.x-orig1.x)^2 + (final1.y-orig1.y)^2 )
    origmag = fltarr(norigmagind)
    for l=0,norigmagind-1 do origmag[l]=orig1.(origmagind[l])
    gdorig = where(finalmag lt 50 and origmag lt 50,ngdorig)
    if ngdorig gt 0 then omagdiff = mean(abs(origmag[gdorig]-finalmag[gdorig])) else omagdiff=99.99
    ofinaldiff = sqrt(odist^2 + omagdiff^2)

    ; SYNTH
    MATCH,aind1,finaldbl1,ind2,/sort
    synthdbl[k] = ind2  ; index into AIND1/2
    synth1 = synth[aind2[ind2]]
    adist = sqrt( (final1.x-synth1.x)^2 + (final1.y-synth1.y)^2 )
    synthmag = fltarr(nsynthmagind)
    for l=0,nsynthmagind-1 do synthmag[l]=synth1.(synthmagind[l])
    gdsynth = where(finalmag lt 50 and synthmag lt 50,ngdsynth)
    if ngdsynth gt 0 then amagdiff = mean(abs(synthmag[gdsynth]-finalmag[gdsynth])) else amagdiff=99.99
    afinaldiff = sqrt(adist^2 + amagdiff^2)

    ; Which one is the match
    if ofinaldiff lt afinaldiff then begin
      com='  REAL'
      flag[k] = 0
    endif else begin
      com='  AST'
      flag[k] = 1
    endelse
    if keyword_set(verbose) then $
      printlog,logfile,'  '+strtrim(i+1,2),odist,omagdiff,ofinaldiff,' ',adist,amagdiff,afinaldiff,com
    ;if amagdiff lt omagdiff and adist gt odist then stop
  endfor  ; duplicates loop

  ; Remove ASTs from the "ORIG" list
  bdomatch = where(flag eq 1,nbdomatch,ncomp=ngdomatch)  
  print,strtrim(nbdomatch,2),' are ASTs and ',strtrim(ngdomatch,2),' are REAL sources'
  if nbdomatch gt 0 then begin
    if nbdomatch lt nomatch then begin
      bdorigdbl = origdbl[bdomatch]
      REMOVE,bdorigdbl,oind1,oind2
      nomatch = n_elements(oind1)
    endif else begin
      undefine,oind1,oind2
      nomatch = 0
    endelse
  endif
endif ; duplicates

; "Prune" the real sources from the list of recovered sources
left = final
leftind = lindgen(n_elements(final))
if nomatch gt 0 then remove,oind1,left,leftind
; Now rematch SYNTH to the leftover sources
SRCMATCH,left.x,left.y,synth.x,synth.y,2,aind1,aind2,count=nastmatch
print,strtrim(nastmatch,2),' ASTs recovered'

; Final indices
find = leftind[aind1]
sind = aind2
nmatch = nastmatch

stop

end
