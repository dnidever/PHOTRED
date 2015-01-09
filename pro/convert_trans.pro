pro convert_trans,infile,outfile,mag,color,stp=stp

;+
;
; CONVERT_TRANS
;
; This converts the SKAWDPHOT transformation equation output
; into something more useful, similar to what STDRED_TRANSPHOT.PRO
; spits out.
;
; The "infile" should be a SKAWDPHOT output file with transformation
; equations of this form.  The last set of equations in the file
; will be used.
;
;   Inversion completed successfully
;   Error of Solution, SIGMA =  6.576E-03
;   Coefficients and errors:
;   a( 1) =     2.8101 sigma =  9.142E-03  Night  1 zero point
;   a( 2) =     2.8534 sigma =  9.231E-03  Night  2 zero point
;   a( 3) =     2.8584 sigma =  9.522E-03  Night  3 zero point
;   a( 4) =     2.8593 sigma =  1.114E-02  Night  4 zero point
;   a( 5) =     2.8144 sigma =  9.119E-03  Night  5 zero point
;   a( 6) =     0.1344 sigma =  5.712E-03  airmass
;   a( 7) =    -0.0387 sigma =  1.858E-03  color
;
;
; INPUTS:
;  infile   The skawdphot output file (e.g. M.data.out).
;  outfile  The output filename (e.g. M.trans)
;  mag      The name of the filter passband (e.g. 'M')
;  color    The name of the color (e.g. 'M-T')
;  /stp     Stop at the end of the program.
; 
; OUTPUTS:
;  The output file of converted transformation equations.  They
;  will be in the PHOTRED/STDRED format.
;
;   T  M-T   3.2027    0.0268   -0.0338   0.0000   0.0000
;            0.0127    0.0072    0.0037   0.0000   0.0000
;
;   Magnitude and color, then the zeropoint, airmass and color
;   terms.  On the second line are the uncertainties in those
;   quantities.
;
; USAGE:
;  IDL>convert_trans,infile,outfile,mag,color,stp=stp
;
; By D. Nidever  May 2008
;-

ninfile = n_elements(infile)
if ninfile eq 0 then begin
  print,'Syntax - convert_trans,infile,outfile,mag,color,stp=stp'
  return
endif

readline,infile,lines,count=count
if count eq 0 then begin
  print,'Problem with ',infile
  return
endif
lines = strtrim(lines,2)
nlines = n_elements(lines)

; Where does the last trans. eqns output start
ind = where(stregex(lines,'Coefficients and errors:',/boolean) eq 1,nind)
lo = max(ind)+1

; Blank lines
blankind = where(lines eq '',nblankind)
if nblankind gt 0 then begin
  gd = where(blankind gt lo,ngd)
  if ngd gt 0 then hi = blankind[gd[0]]-1 else hi=nlines-1
endif else begin
  hi = nlines-1
endelse


; Starting the final output lines array
undefine,flines
push,flines,'SKAWDPHOT OUTPUT'
push,flines,''
push,flines,lines[lo-3:hi]
push,flines,''


; Getting the transformation equations
tlines = lines[lo:hi]

nightind = where(stregex(tlines,'Night',/boolean) eq 1,nnightind)

tlines2 = strmid(tlines,8,100)
tlines2 = REPSTR(tlines2,'sigma =')

tarr = strsplitter(tlines2,' ',/extract)
par = float(reform(tarr[0,*]))
sigma = float(reform(tarr[1,*]))

airind = where(stregex(tlines2,'airmass',/boolean) eq 1,nairind)
colind = where(stregex(tlines2,'color',/boolean) eq 1,ncolind)

airterm = par[airind[0]]
airsig = sigma[airind[0]]
colterm = par[colind[0]]
colsig = sigma[colind[0]]

nightnum = strtrim(reform(tarr[3,nightind]),2)
nightterm = par[nightind]
nightsig = sigma[nightind] 

; OUTPUT

if n_elements(mag) eq 0 then magname='MAG' else magname=mag
if n_elements(color) eq 0 then colname='COL' else colname=color
magname = strtrim(strupcase(magname),2)
colname = strtrim(strupcase(colname),2)

; Loop through the nights
fmt2 = '(F7.4)'
for i=0,nnightind-1 do begin
  push,flines,'Night '+nightnum[i]+' Transformation Equations'
  add='  '+magname+'  '+colname+'  '
  nadd = strlen(add)
  magline = add+string(nightterm[i],format=fmt2)+'   '+string(airterm,format=fmt2)+$
       '   '+string(colterm,format=fmt2)+'   0.0000   0.0000'
  push,flines,magline
  spaces = string(byte(lonarr(nadd)+32))
  errline = spaces+string(nightsig[i],format=fmt2)+'   '+string(airsig,format=fmt2)+$
       '   '+string(colsig>0.0001,format=fmt2)+'   0.0000   0.0000'
  push,flines,errline
  if i ne nnightind-1 then push,flines,''
end

; Print everything to the screen
printline,flines

; Output
if n_elements(outfile) gt 0 then begin
  WRITELINE,outfile,flines
  print,'Trans equations written to >>',outfile,'<<'
endif else begin
  print,'NO output file written'
endelse

if keyword_set(stp) then stop

end
