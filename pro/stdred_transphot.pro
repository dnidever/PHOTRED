pro std_trans_dummy
; This makes it so that you don't have to compile before running
FORWARD_FUNCTION std_devtrans, std_transfunc
end

;--------------------

function std_devtrans,p,y=y,mag=mag,col=col,am=am,night=night,chip=chip,$
                      mapnight=mapnight,weight=weight

model = std_transfunc(p,mag=mag,col=col,am=am,night=night,chip=chip,mapnight=mapnight)

return,(y-model)*weight

end

;---------------------------

function std_transfunc,p,mag=mag,col=col,am=am,night=night,chip=chip,mapnight=mapnight


;nnights = max(night)-min(night)+1
;zeropar = P[0:nnights-1]
;ampar = P[nnights]
;colpar = P[nnights+1]
;nightzero = zeropar[night-min(night)]

npar = n_elements(p)

;zeropar = P[0:npar-3]
;ampar = P[npar-2]
;colpar = P[npar-1]
zeropar = P[0:npar-4]
ampar = P[npar-3]
colpar = P[npar-2]
amcolpar = P[npar-1]

nightzero = zeropar[mapnight[night]]

;   V = mV - v1 - v2 * XV - v3 * (B-V) - v4 * XV*(B-V) - v5 * (B-V) * (B-V)
;         +(aperture correction) + (time correction)
;model = mag + nightzero + col*colpar + am*ampar
model = mag - nightzero - col*colpar - am*ampar - col*am*amcolpar

return,model

end 

;---------------------

pro stdred_transphot,input,stp=stp,arr=arr,plotresid=plotresid,$
                  yrange=yrange,fixam=fixam,fixcol=fixcol,trans=trans,$
                  tlines=tlines,rms=rms,errlim=errlim,nooutput=nooutput,$
                  inparr=inparr,sepchip=sepchip,fixac=fixac,fitac=fitac,silent=silent

;+
;
; STDRED_TRANSPHOT
;
; This is very similar to FIT_TRANSPHOT.PRO, but slightly modified
; for the STDRED pipeline.
;
; This program finds the photometric transformation equations
; with iterative outlier rejection given the proper input file
; (same as SKAWDPHOT.PRO).
;
; INPUTS:
;  input    Filename of input data (i.e. M.data).  Same as for
;            SKAWDPHOT.PRO.
;  =fixam   Fix the airmass term to this value.
;  =fixcol  Fix the color term to this value.
;  =fixac   Fix the airmass*color term to this value.  The default is 0.0.
;  /fitac   Fit the airmass*color term.  The default is keep it
;             fixed at 0.0.
;  =errlim  Only use stars with observed errors lower than this value.
;  =inparr  Input the structure instead of the filename
;  /plotresid  Plot the residuals.
;  =yrange     Yrange for residual plotting.
;  /nooutput  Don't print out any of the transformation files.
;  /sepchip  Each chip separately.
;  /stp     Stop at the end of the program.
;
; OUTPUTS:
;  The transformation equations are written to a file called
;  input+'.trans'.  Separate transformation equations for each
;  night are also put in the file.  Also, each night's transformation
;  equations are written to their own files called 'n#MAG.trans'
;  (i.e. n1M.trans').
;
; =arr      The data structure.
; =tlines   An array of the transformation equations for each night.
; =trans    The transformation structure, one element for each night.  
; =rms      The final rms of the fit.
;
; USAGE:
;  IDL>fit_transphot,'M.data',
;
; By D. Nidever   Jan. 2008
;-

; Not enough inputs
if n_elements(input) eq 0 and n_elements(inparr) eq 0 then begin
  print,'Syntax - stdred_transphot,input,inparr=inparr,stp=stp,plotresid=plotresid,arr=arr'
  print,'                          fixam=fixam,fixcol=fixcol,fixac=fixac,fitac=fitac,sepchip=sepchip'
  return
endif

; Loading the data
if n_elements(input) gt 0 then if input ne '' then begin
  arr = importascii(input[0],/header,/noprint)

; Using structure input
endif else begin
  arr = inparr
endelse
narr = n_elements(arr)
orig = arr

; Each chip separately
;  This is a quick KLUDGE
if keyword_set(sepchip) and tag_exist(arr,'CHIP') then begin
  print,'Solving all CHIPS separately'

  print,'' & print,'Initial solution of each chip separately'
  ui = uniq(arr.chip,sort(arr.chip))
  chips = arr[ui].chip
  nchips = n_elements(chips)
  undefine,alltrans1
  for i=0,nchips-1 do begin
    print,'' & print,'Chip = ',strtrim(chips[i],2)
    ind = where(arr.chip eq chips[i],nind)
    inparr = arr[ind]
    stdred_transphot,'',stp=stp,plotresid=plotresid,$
                 yrange=yrange,fixam=fixam,fixcol=fixcol,fixac=fixac,$
                 fitac=fitac,trans=trans,$
                 tlines=tlines,rms=rms,errlim=errlim,$
                 inparr=inparr,/nooutput ;,/silent
    add_tag,trans,'chip',chips[i],trans
    PUSH,alltrans1,trans
  endfor

  ; Map from chip number to index in alltrans
  mapchip = lonarr(max(chips)+1)-1
  mapchip[chips] = indgen(nchips)

  ; Remove zero-point offsets and fit color-term, airmass and airmass*color terms
  print,'' & print,'Removing zero-point offsets and solving for color, airmass and airmass*color terms'
  zpterm = alltrans1[mapchip[arr.chip]].zpterm
  tarr = arr
  tarr.mag -= zpterm
  stdred_transphot,'',stp=stp,plotresid=plotresid,$
                 yrange=yrange,fixam=fixam,fixcol=fixcol,fixac=fixac,fitac=fitac,$
                 trans=trans0,tlines=tlines,rms=rms,errlim=errlim,$
                 inparr=tarr,/nooutput  ;,/silent
  print,'Airmass term = ',trans0.amterm
  print,'Color term = ',trans0.colterm
  print,'Airmass*color term = ',trans0.amcolterm
  
  ; Now fix color-term, airmass and airmass-color terms and refit chip zeropoints
  print,'' & print,'Refitting zero-point offsets and keeping color, airmass and airmass*color terms fixed'
  fixam = trans0[0].amterm
  fixcol = trans0[0].colterm
  fixac = trans0[0].amcolterm
  undefine,alltrans
  for i=0,nchips-1 do begin
    print,'' & print,'Chip = ',strtrim(chips[i],2)
    ind = where(arr.chip eq chips[i],nind)
    inparr = arr[ind]
    stdred_transphot,'',stp=stp,plotresid=plotresid,$
                 yrange=yrange,fixam=fixam,fixcol=fixcol,fixac=fixac,$
                 trans=trans,tlines=tlines,rms=rms,errlim=errlim,$
                 arr=arr1,inparr=inparr,/nooutput  ;,/silent
    push,allarr,arr1
    add_tag,trans,'chip',chips[i],trans
    PUSH,alltrans,trans
  endfor
  ; Calculate final RMS
  gd = where(allarr.rejected eq 0,npts)
  npts = n_elements(tarr.mag)
  totwt = total(allarr[gd].weight)
  rms = sqrt( total(allarr[gd].weight*allarr[gd].resid^2.)*npts/((npts-1.)*totwt) )

  ; Write to file

  ; Print out the transformation equations
  ;---------------------------------------
  fmt2 = '(F7.4)'
  undefine,lines
  push,lines,'#'
  push,lines,'# FINAL TRANSFORMATION COEFFICIENTS:'
  push,lines,'# Final RMS      = '+string(rms,format='(F9.6)')
  push,lines,'#'
  ;for i=0,nnights-1 do begin
  ;  push,lines,'Night '+strtrim(unights[i],2)+' zpoint = '+string(fpar[i],format=fmt2)+$
  ;             ' ('+string(sigpar[i],format=fmt2)+')'
  ;end
  push,lines,'# Airmass term   = '+string(alltrans[0].amterm,format=fmt2)+$
             '#  ('+string(alltrans[0].amerr,format=fmt2)+')'
  push,lines,'# Color term     = '+string(alltrans[0].colterm,format=fmt2)+$
             '#  ('+string(alltrans[0].colerr,format=fmt2)+')'
  push,lines,'# Am*Color term  = '+string(alltrans[0].amcolterm,format=fmt2)+$
             '#  ('+string(alltrans[0].amcolerr,format=fmt2)+')'
  if not keyword_set(silent) then printline,lines

  ; Figure out what the magnitude and color names are
  tags = TAG_NAMES(arr)
  magname = tags[6]
  colname = tags[8]
  colname = REPSTR(colname,'_','-')


  ; Print transformation equations to file
  ;---------------------------------------
  ;    M    M-T  -0.9990    0.1402     -0.1345    0.0000   0.0000
  ;              1.094E-02  5.037E-03  2.010E-03  0.0000   0.0000

  ; Final transformation equations
  if not keyword_set(nooutput) then $
    WRITELINE,magname+'.trans',lines

  ; Transformation equations for each night and chip
  ui = uniq(alltrans.night,sort(alltrans.night))
  unights = alltrans[ui].night
  nnights = n_elements(unights)
  for i=0,nnights-1 do begin
    print,'Night ',strtrim(unights[i],2)
    undefine,lines1,all
    tlines = strarr(nchips,2)
    for j=0,nchips-1 do begin
      tind = where(alltrans.night eq unights[i] and alltrans.chip eq chips[j],ntind)
      itrans = alltrans[tind[0]]

      undefine,lines1
      add=string(itrans.chip,format='(I3)')+'  '+magname+'  '+colname+'  '
      nadd = strlen(add)
      magline = add+string(itrans.zpterm,format=fmt2)+'   '+string(itrans.amterm,format=fmt2)+$
           '   '+string(itrans.colterm,format=fmt2)+'   '+string(itrans.amcolterm,format=fmt2)+'   0.0000'
      push,lines1,magline
      spaces = string(byte(lonarr(nadd)+32))
      errline = spaces+string(itrans.zperr,format=fmt2)+'   '+string(itrans.amerr,format=fmt2)+$
           '   '+string(itrans.colerr>0.0001,format=fmt2)+'   '+string(itrans.amcolerr>0.0001,format=fmt2)+'   0.0000'
      push,lines1,errline

      ; Add to trans lines array
      tlines[i,*] = lines1

      push,all,' '
      push,all,lines1
    endfor
    printline,all
    print,''
    WRITELINE,'n'+strtrim(unights[i],2)+magname+'.trans',all

    ; Add to MAG.trans file
    WRITELINE,magname+'.trans',['','Night '+strtrim(unights[i],2)],/append
    WRITELINE,magname+'.trans',all,/append
  endfor


  ; Output the stars used for deriving the transformation
  ; with the REJECTED tag
  arr = allarr
  
  return

endif  ; each chip separately


; Add CMAG, COL to arr
ADD_TAG,arr,'CMAG',0.0,arr
arr.cmag = arr.(6)
ADD_TAG,arr,'COL',0.0,arr
arr.col = arr.(8)

; Add a rejection tag
ADD_TAG,arr,'rejected',0,arr
arr.rejected = 0

; Select only stars with good instrumental and calibrated photometry
gdphot = where(arr.mag gt 0 and arr.mag lt 50 and $
               arr.cmag gt 0 and arr.cmag lt 50 and $
               arr.col gt -5 and arr.col lt 20,ngdphot)
if ngdphot eq 0 then stop,'NO stars'
arr = arr[gdphot]


; Error limit
if n_elements(errlim) gt 0 then if errlim gt 0.0001 then begin
  if not keyword_set(silent) then $
    print,'Imposing Error limit <= ',strtrim(errlim,2)
  ;gderr = where(arr.err le errlim,ngderr)
  gderr = where(arr.err le errlim and arr.(7) le errlim,ngderr) ; instr and cal errors
  if ngderr eq 0 then gderr = where(arr.err le 2*errlim,ngderr)
  if ngderr eq 0 then stop,'no stars'
  if not keyword_set(silent) then $
    print,strtrim(ngderr,2),'/',strtrim(narr,2),' observations retained'
  arr = arr[gderr] 
endif


; Figure out what the magnitude and color names are
tags = TAG_NAMES(arr)
magname = tags[6]
colname = tags[8]
colname = REPSTR(colname,'_','-')


; LO and HI indices for the unique stars
si = sort(arr.id)
arr = arr[si]
idlo = where(shift(arr.id,1) ne arr.id)
idhi = where(shift(arr.id,-1) ne arr.id)


;; Get a number for each unique "real" star
; Use the calibrated color and magnitude to do this.
ADD_TAG,arr,'REALSTAR',0L,arr
ui = uniq(arr.id,sort(arr.id))
uniqid = arr[ui].id
nuniqid = n_elements(ui)
for i=0,nuniqid-1 do begin
  gd = where(arr.id eq uniqid[i],ngd)
  arr[gd].realstar=i+1
end


; Nights
; mapnight[night] gives the parameter zero-point
; index for the "night"
; zeronight = par[mapnight[night]]
ui = uniq(arr.night,sort(arr.night))
unights = arr[ui].night
nnights = n_elements(unights)
mapnight = lonarr(max(unights)+1)-1
mapnight[unights] = indgen(nnights)

; Chips
; mapchip[chip] gives the parameter zero-point
; index for the "chip"
; zerochip = par[mapchip[chip]]
ui = uniq(arr.chip,sort(arr.chip))
uchips = arr[ui].chip
nchips = n_elements(uchips)
mapchip = lonarr(max(uchips)+1)-1
mapchip[uchips] = indgen(nchips)

; Outlier rejection loop
flag=0
count=0
nnewrej=0
nrej=0
;nnights = max(arr.night)-min(arr.night)+1
;par = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
;par = replicate(0.0,nnights+2)
par = replicate(0.0,nnights+3)
if not keyword_set(silent) then begin
  print,'-------------------------------------'
  print,' ITER  NPTS   RMS       SIG     NREJ  '
  print,'====================================='
endif
WHILE (flag ne 1) do begin

  ; Get non-rejected stars
  gd = where(arr.rejected eq 0,ngd)
  tarr = arr[gd]

  ; Fitting with MPFIT.PRO
  func = 'std_devtrans'
  ;fa = {y:arr.mag, mag:arr.cmag, col:arr.col, am:arr.airmass, night:arr.night, weight:arr.weight}
  if tag_exist(arr,'CHIP') then chip=tarr.chip else chip=lonarr(ngd)+1
  fa = {y:tarr.cmag, mag:tarr.mag, col:tarr.col, am:tarr.airmass,$
        night:tarr.night, chip:chip, mapnight:mapnight, weight:tarr.weight}

  npar = n_elements(par)
  parinfo = replicate({limited:[0,0],limits:[0.0,0.0],fixed:0},npar)
  ; zeropoint
  parinfo[npar-4].limited=1 & parinfo[npar-4].limits=[-5,5]
  ; airmass
  parinfo[npar-3].limited=1 & parinfo[npar-3].limits=[-5,5]
  if n_elements(fixam) gt 0 then begin
    par[npar-3] = fixam
    parinfo[npar-3].fixed = 1
  endif
  ; color
  parinfo[npar-2].limited=1 & parinfo[npar-2].limits=[-5,5]
  if n_elements(fixcol) gt 0 then begin
    par[npar-2] = fixcol
    parinfo[npar-2].fixed = 1
  endif
  parinfo[npar-1].fixed = 1  ; fix am*col by default
  if keyword_set(fitac) then parinfo[npar-1].fixed=0  ; let am*col float
  if n_elements(fixac) gt 0 then begin
    par[npar-1] = fixac
    parinfo[npar-1].fixed = 1
  endif  

  par[0:nnights-1] = median([tarr.mag-tarr.cmag])

  fpar = MPFIT(func,par,functargs=fa,status=status,perror=perror,bestnorm=chisq,$
               parinfo=parinfo,dof=dof,autoderivative=1,ftol=ftol,/quiet)

  sigpar = perror * sqrt(chisq/dof)

  ; total resid
  model = std_transfunc(fpar,mag=tarr.mag,col=tarr.col,am=tarr.airmass,night=tarr.night,mapnight=mapnight)
  ;resid = arr.mag-model
  resid = model-tarr.cmag
  ;rms = sqrt(total(resid^2)/(npts-n_elements(fpar)-1))
  npts = n_elements(tarr.mag)
  totwt = total(tarr.weight)
  rms = sqrt( total(tarr.weight*resid^2.)*npts/((npts-1.)*totwt) )
  sig = mad(resid)
  
  model2 = std_transfunc(fpar,mag=arr.mag,col=arr.col,am=arr.airmass,night=arr.night,mapnight=mapnight)
  resid2 = model2-arr.cmag


  ; New Rejecting stars with bad resids
  thresh = 3.0*sig > 0.015 
  bad1 = where(abs(resid2-median(resid2,/even)) gt thresh and arr.rejected eq 0,nbad1)
  if nbad1 gt 0 then begin
    arr[bad1].rejected = 1
    ;remove,bd,arr
    par = fpar
  endif

  ; Chuck any stars that are REALLY off
  nbad2 = 0
  com=''
  nstar_rej = 0
  for i=0,nuniqid-1 do begin
    nind = idhi[i]-idlo[i]+1
    ind = lindgen(nind)+idlo[i]
    ;ind = where(arr.id eq uniqid[i],nind)
    medresid = median(resid2[ind],/even)
    if (abs(medresid) gt 2.0*sig and min(arr[ind].rejected) eq 0 and nind gt 2) then begin
      newbad = where(arr[ind].rejected eq 0,nbad2)
      arr[ind].rejected = 1
      
      if com eq '' then com='Star(s) ' else com=com+', '
      ;com = com+strtrim(uniqid[i],2)

      nstar_rej++
    endif
  end
  if nstar_rej gt 0 then com=strtrim(nstar_rej,2)+' stars rejected'
  ;if com ne '' then com=com+' rejected'

  count++


  ; How many new rejects
  nnewrej = nbad1+nbad2
  ; Ending?
  if nnewrej eq 0 then flag=1

  fmt = '(I3,I6,F10.5,F10.5,I5,A5,A-45)'
  ;print,format=fmt,count,npts,rms,sig,nbd
  if not keyword_set(silent) then $
    print,format=fmt,count,npts,rms,sig,nnewrej,'',com

  ;stop

end
if not keyword_set(silent) then $
  print,'-------------------------------------'

;nnights = max(arr.night)-min(arr.night)+1


; Put the residual information into the structure
ADD_TAG,arr,'resid',0.0,arr
model = std_transfunc(fpar,mag=arr.mag,col=arr.col,am=arr.airmass,night=arr.night,mapnight=mapnight)
resid = model-arr.cmag
arr.resid = resid


; Print out the number of "good" observations per night
;------------------------------------------------------
if not keyword_set(silent) then begin
  print,''
  print,'Number of "good" observations per night'
  print,'---------------------------------------'
endif
for i=0,nnights-1 do begin
  inight = unights[i]
  dum = where(arr.night eq inight and arr.rejected eq 0,ngd_inight)
  if not keyword_set(silent) then $
    print,'Night '+strtrim(inight,2)+' Nobs = '+strtrim(ngd_inight,2)
end

; Print out RMS per night
;------------------------
if not keyword_set(silent) then begin
  print,''
  print,'RMS per night'
  print,'-------------'
endif
for i=0,nnights-1 do begin
  inight = unights[i]
  igd = where(arr.night eq inight and arr.rejected eq 0,ngd_inight)
  irms = sqrt(median([arr[igd].resid]^2.0))
  if ngd_inight eq 1 then irms = 0.0
  if not keyword_set(silent) then $
    print,'Night '+strtrim(inight,2)+' RMS = '+string(irms,format='(F10.5)')
end
if not keyword_set(silent) then print,''

; Plot the residuals
;---------------------------------
if keyword_set(plotresid) then begin

  ;window,15,xsize=600,ysize=400
  !p.multi = [0,3,2]
  charsize = 2.0
  bd = where(arr.rejected eq 1,nbd)

  dy = range(arr.resid)*0.1
  if keyword_set(yrange) then if n_elements(yrange) eq 2 then yr=yrange
  if n_elements(yr) eq 0 then yr = [min(arr.resid)-dy,max(arr.resid)+dy]

  ; Resid vs. Airmass
  ;------------------
  dx = range(arr.airmass)*0.1
  xr = [min(arr.airmass)-dx,max(arr.airmass)+dx]
  plot,arr.airmass,arr.resid,ps=1,xtit='Airmass',ytit='Residuals',tit='Airmass',$
       xr=xr,yr=yr,xs=1,ys=1,charsize=charsize
  if nbd gt 0 then oplot,[arr[bd].airmass],[arr[bd].resid],ps=1,co=250

  ; Resid vs. Color
  ;----------------
  dx = range(arr.col)*0.1
  xr = [min(arr.col)-dx,max(arr.col)+dx]
  plot,arr.col,arr.resid,ps=1,xtit='Color',ytit='Residuals',tit='Color',$
       xr=xr,yr=yr,xs=1,ys=1,charsize=charsize
  if nbd gt 0 then oplot,[arr[bd].col],[arr[bd].resid],ps=1,co=250

  ; Resid vs. Magnitude
  ;--------------------
  dx = range(arr.cmag)*0.1
  xr = [min(arr.cmag)-dx,max(arr.cmag)+dx]
  plot,arr.cmag,arr.resid,ps=1,xtit='Magnitude',ytit='Residuals',tit='Magnitude',$
       xr=xr,yr=yr,xs=1,ys=1,charsize=charsize
  if nbd gt 0 then oplot,[arr[bd].cmag],[arr[bd].resid],ps=1,co=250

  ; Resid vs. night
  ;----------------
  xr = [min(arr.night)-1,max(arr.night)+1]
  plot,arr.night,arr.resid,ps=1,xtit='Night',ytit='Residuals',tit='Night',$
       xr=xr,yr=yr,xs=1,ys=1,charsize=charsize,xminor=1
  if nbd gt 0 then oplot,[arr[bd].night],[arr[bd].resid],ps=1,co=250

  ; Resid vs. frame
  ;----------------
  ;frame = strarr(narr)
  ;for i=0,narr-1 do frame[i] = first_el(strsplit(arr[i].id,'-',/extract))
  ui = uniq(arr.frame,sort(arr.frame))
  uframes = arr[ui].frame
  nframes = n_elements(uframes)

  if nframes le 60 then begin

    if not keyword_set(silent) then begin
      print,' NUM      FRAME'
      for i=0,nframes-1 do print,format='(I3,A12)',i+1,uframes[i]
    endif

    xtickv = indgen(nframes+1)
    xtickname = [' ',strtrim(indgen(nframes)+1,2),' ']

    xr = [0,nframes+1]
    plot,[0],[0],/nodata,ps=1,xtit='Frame',ytit='Residuals',tit='Frame',$
         xr=xr,yr=yr,xs=1,ys=1,charsize=charsize,xtickinterval=1,$
         xtickv=xtickv,xtickname=xtickname,xminor=1
    for i=0,nframes-1 do begin
      frameind = where(arr.frame eq uframes[i],nframeind)
      oplot,[lonarr(nframeind)+i+1],[arr[frameind].resid],ps=1
      bdframeind = where(arr.frame eq uframes[i] and arr.rejected eq 1,nbdframeind)
      if nbdframeind gt 0 then oplot,[lonarr(nbdframeind)+i+1],[arr[bdframeind].resid],ps=1,co=250
    end

  endif else begin
    print,'Too many frames to plot RESID vs. FRAME'
  endelse

  ; Resid vs. realstar
  ;-------------------
  xr = [0,max(arr.realstar)+1]
  plot,arr.realstar,arr.resid,ps=1,xtit='Star',ytit='Residuals',tit='Star',$
       xr=xr,yr=yr,xs=1,ys=1,charsize=charsize
  if nbd gt 0 then oplot,[arr[bd].realstar],[arr[bd].resid],ps=1,co=250

  ui = uniq(arr.realstar,sort(arr.realstar))
  urealstar = arr[ui].realstar
  nrealstar = n_elements(urealstar)

  if not keyword_set(silent) then begin
    print,''
    print,' NUM  STAR ID'
  endif
  for i=0,nrealstar-1 do begin
    starind = where(arr.realstar eq urealstar[i],nstarind)
    if not keyword_set(silent) then $
      print,format='(I3,'+strtrim(nstarind,2)+'A12)',i+1,arr[starind[0]].id
  end

  !p.multi=0

endif



; Print out the transformation equations
;---------------------------------------
fmt2 = '(F7.4)'
undefine,lines
push,lines,'# '
push,lines,'# FINAL TRANSFORMATION COEFFICIENTS:'
push,lines,'# Final RMS      = '+string(rms,format='(F9.6)')
;push,lines,'Final RMS      = '+string(rms,format=fmt2)
push,lines,'# '
for i=0,nnights-1 do begin
  ;push,lines,'Night '+strtrim(i+min(arr.night),2)+' zpoint = '+string(fpar[i],format=fmt2)+$
  ;           ' ('+string(sigpar[i],format=fmt2)+')'
  push,lines,'# Night '+strtrim(unights[i],2)+' zpoint = '+string(fpar[i],format=fmt2)+$
             ' ('+string(sigpar[i],format=fmt2)+')'
end
push,lines,'# Airmass term   = '+string(fpar[nnights],format=fmt2)+' ('+string(sigpar[nnights],format=fmt2)+')'
push,lines,'# Color term     = '+string(fpar[nnights+1],format=fmt2)+' ('+string(sigpar[nnights+1],format=fmt2)+')'
push,lines,'# Am*Color term  = '+string(fpar[nnights+2],format=fmt2)+' ('+string(sigpar[nnights+2],format=fmt2)+')'
if not keyword_set(silent) then printline,lines



; Print transformation equations to file
;---------------------------------------
;    M    M-T  -0.9990    0.1402     -0.1345    0.0000   0.0000
;              1.094E-02  5.037E-03  2.010E-03  0.0000   0.0000

; Final transformation equations
if not keyword_set(nooutput) then $
  WRITELINE,magname+'.trans',lines


; Start the trans structure
;  one element for each night
transdum = {magname:magname,colname:colname,night:0L,rms:rms,zpterm:0.0,zperr:0.0,$
            amterm:fpar[nnights],amerr:sigpar[nnights],colterm:fpar[nnights+1],$
            colerr:sigpar[nnights+1],amcolterm:fpar[nnights+2],amcolerr:sigpar[nnights+2],$
            magline:'',errline:''}
trans = REPLICATE(transdum,nnights)


; Transformation equations for each night
undefine,lines1,all
tlines = strarr(nnights,2)
for i=0,nnights-1 do begin

  undefine,lines1
  ;push,lines1,' '
  ;push,lines1,'Night '+strtrim(i+1,2)+' Transformation Equations'
  add='  '+magname+'  '+colname+'  '
  nadd = strlen(add)
  magline = add+string(fpar[i],format=fmt2)+'   '+string(fpar[nnights],format=fmt2)+$
       '   '+string(fpar[nnights+1],format=fmt2)+'   '+string(fpar[nnights+2],format=fmt2)+'   0.0000'
  push,lines1,magline
  spaces = string(byte(lonarr(nadd)+32))
  errline = spaces+string(sigpar[i],format=fmt2)+'   '+string(sigpar[nnights],format=fmt2)+$
       '   '+string(sigpar[nnights+1]>0.0001,format=fmt2)+'   '+$
       string(sigpar[nnights+2]>0.0001,format=fmt2)+'   0.0000'
  push,lines1,errline

  ; Writing transformation equations for this night only
  ;WRITELINE,'n'+strtrim(i+1,2)+magname+'.trans',lines1
  if not keyword_set(nooutput) then $
    WRITELINE,'n'+strtrim(unights[i],2)+magname+'.trans',lines1

  ; Add to trans structure
  trans[i].night = unights[i]
  trans[i].zpterm = fpar[i]
  trans[i].zperr = sigpar[i]
  trans[i].magline = magline
  trans[i].errline = errline
  ; Add to trans lines array
  tlines[i,*] = lines1

  push,all,' '
  ;push,all,'Night '+strtrim(i+1,2)+' Transformation Equations'
  push,all,'Night '+strtrim(unights[i],2)+' Transformation Equations'
  push,all,lines1
endfor

; Write everything to a log file
if not keyword_set(nooutput) then begin
  WRITELINE,magname+'.trans',all,/append

  if not keyword_set(silent) then begin
    print,''
    print,'Transformation equations written to: ',magname+'.trans'
    print,''
  endif
endif


if keyword_set(stp) then stop

end
