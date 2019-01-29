;+
;
; MATCHSTARS
;
; This program matches two star lists of X/Y coordinates
; using cross-correlation of two "fake" images
;
; INPUTS:
;  xinp1   The X-coordinates of the first set of stars
;  yinp1   The Y-coordinates of the first set of stars
;  xinp2   The X-coordinates of the second set of stars
;  yinp2   The Y-coordinates of the second set of stars
;  /plotresid   Plot the residuals at each stag
;  /stp    Stop at the end of the program
;
; OUTPUTS:
;  ind1    The matched indices for the first set
;  ind2    The matched indices for the second set
;  trans   The transformation equation to go from (2) to (1).
;            trans = [A, B, C, D, E, F]          
;            x(1) = A + C*x(2) + E*y(2)
;            y(1) = B + D*x(2) + F*y(2)
;  =count  The number of matches.  This is set to -1 if there
;            was an error.
;
; USAGE:
;  IDL>matchstars,xinp1,yinp1,xinp2,yinp2,ind1,ind2,trans
;
; By D.Nidever  Feb. 2008
;-

pro matchstars_dummy
FORWARD_FUNCTION trans_coo, trans_coo_dev
end


;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function trans_coo,xin,yin,par

; Apply the transformation to X/Y

A = par[0]
B = par[1]
C = par[2]
D = par[3]
E = par[4]
F = par[5]

; from ccdpck.txt
;              x(1) = A + C*x(n) + E*y(n)
;              y(1) = B + D*x(n) + F*y(n)

; Apply transformation
xout = A + C*xin + E*yin
yout = B + D*xin + F*yin


nstars = n_elements(xin)
out = dblarr(2,nstars)
out[0,*] = xout
out[1,*] = yout

return,out

end


;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function trans_coo_dev,par,x1=x1,y1=y1,x2=x2,y2=y2

; Rotate coordinates(2) to coordinate system 1
; and return deviates

out = trans_coo(x2,y2,par)
nx2 = reform(out[0,*])
ny2 = reform(out[1,*])

diff = sqrt( (x1-nx2)^2.0d0 + (y1-ny2)^2.0d0 )

; Do robust outlier rejection
std = mad(diff)
med = median([diff],/even)
;rms = sqrt(mean(diff^2.0))
bd = where(diff gt (med+3.0*std),nbd)
if nbd gt 0 then diff[bd] = 0.0

return,diff

end


;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


pro matchstars_xcorr,xx1,yy1,xx2,yy2,xshift,yshift,angle,bestcorr,xcorr,gradient=gradient,$
                     smooth=smooth0,xyscale=xyscale0,fwhm=fwhm0,xoff=xoffs,yoff=yoffs,$
                     fft1=fft1,fft2=fft2,fftp=fftp,matchnum=matchnum,nsig=nsig,stp=stp,$
                     maxshift=maxshift,extra=extra,silent=silent

COMMON xcorr, factorsum, factornum

; This program cross-correlates two coordinate lists
;
; The IDL Astro Users's library has
; CORREL_IMAGES, CORREL_OPTIMIZE, and CORRMAT_ANALYZE functions,
; but these DO NOT use FFTs.  It's very slow!


; Error Handling
;------------------
; Establish error handler. When errors occur, the index of the  
; error is returned in the variable Error_status:  
CATCH, Error_status 

;This statement begins the error handler:  
if (Error_status ne 0) then begin 
   xshift = 999999.
   yshift = 999999.
   angle = 999999.
   xcorr = -1
   bestcorr = -1
   matchnum = -1
   nsig = -1
   PHOTRED_ERRORMSG
   CATCH, /CANCEL 
   return
endif

nxx1 = n_elements(xx1)
nyy1 = n_elements(yy1)
nxx2 = n_elements(xx2)
nyy2 = n_elements(yy2)

xshift = 999999.  ; initalize output values
yshift = 999999.
angle = 999999.
xcorr = -1
bestcorr = -1
matchnum = -1
nsig = -1

; Not enough inputs
if (nxx1 eq 0 or nyy1 eq 0 or nxx2 eq 0 or nyy2 eq 0) then begin
  print,'Syntax - matchstars_xcorr,xx1,yy1,xx2,yy2,xshift,yshift,angle,bestcorr,xcorr,gradient=gradient,'
  print,'                     smooth=smooth0,xyscale=xyscale0,fwhm=fwhm0,xoff=xoffs,yoff=yoffs,'
  print,'                     fft1=fft1,fft2=fft2,fftp=fftp,matchnum=matchnum,nsig=nsig,stp=stp,'
  print,'                     maxshift=maxshift,extra=extra'
  return
endif

; xx1/yy1 must have same number of elements (same for xx2/yy2)
if nxx1 ne nyy1 or nxx2 ne nyy2 then begin
  if nxx1 ne nyy1 then $
    print,'XX1/YY1 must have the same number of elements'
  if nxx2 ne nyy2 then $
    print,'XX1/YY1 must have the same number of elements'
  return
endif

; Need at least two stars
if nxx1 lt 2 or nxx2 lt 2 then begin
  if not keyword_set(silent) then print,'Need 2 or more points for matching.'
  return  
endif


; Making all of the coordinates positive
minx = MIN([xx1,xx2])
miny = MIN([yy1,yy2])
xx1s = xx1 - minx
yy1s = yy1 - miny
xx2s = xx2 - minx
yy2s = yy2 - miny

; Scaling the coordinates to make the "fake" images smaller
; and the FFTs faster
if n_elements(xyscale0) gt 0 then xyscale=xyscale0  else xyscale=4
;xyscale = 4
xx1s = xx1s/xyscale
yy1s = yy1s/xyscale
xx2s = xx2s/xyscale
yy2s = yy2s/xyscale
;xx1s = xx1/xyscale
;yy1s = yy1/xyscale
;xx2s = xx2/xyscale
;yy2s = yy2/xyscale

maxxs = max([xx1s,xx2s])
maxys = max([yy1s,yy2s])
nxs = ceil(maxxs) + 1
nys = ceil(maxys) + 1
nxorig = nxs
nyorig = nys

; The run time of the FFT is roughly proportional to the number of
; total elements in the array times the sum of its prime factors
; use FACTOR to factor numbers
; test ~10 numbers above/below the nx/ny numbers and sum the factors
; use the ones with the smallest sums   
; ALLOWING THE ARRAY TO BE LARGER THAN THE NUMBER, AND ALLOWING THERE
; TO BE BLANK SPACE MESSES UP THE CROSS-CORRELATION
if n_elements(fft1) eq 0 and n_elements(fft2) eq 0 and n_elements(fftp) eq 0 then begin

  idldir = file_search('~/idl/',/fully_qualify_path)
  test = file_test(idldir+'/factorsum.dat')
  if test eq 1 and n_elements(factornum) eq 0 then $
    restore,idldir+'/factorsum.dat'

  if n_elements(factornum) gt 0 and n_elements(factorsum) gt 0 then begin
    lo1 = round(nxs*0.80)
    hi1 = nxs 
    ;lo1 = ( round(nxs*0.80) > 100 )
    ;hi1 = round(nxs*1.20) 
    bestind1 = first_el(minloc(factorsum[lo1:hi1])) + lo1
    nxs = factornum[bestind1]
    lo2 = round(nys*0.80)
    hi2 = nys
    ;lo2 = ( round(nys*0.80) > 100 )
    ;hi2 = round(nys*1.20) 
    bestind2 = first_el(minloc(factorsum[lo2:hi2])) + lo2
    nys = factornum[bestind2]

  ; factorsum file not found
  endif else begin

    if not keyword_set(silent) then $
      print,'factorsum.dat NOT FOUND.  PLEASE copy factorsum.dat to ~/idl/'

    ntest = (100 < round(0.2*nxs))    ; 51 ;21 
    ;ntest = 201 ; 51 ;21 
    sumarr1 = fltarr(ntest)
    testarr1 = nxs - findgen(ntest)
    ;testarr1 = findgen(ntest)-(ntest-1)*0.5+nxs
    for i=0,ntest-1 do begin
      factor,testarr1[i],p,nn,/quiet
      sumarr1[i] = total(p*nn)
    end
    best1 = first_el(minloc(sumarr1))
    nxs = testarr1[best1]

    ntest = (100 < round(0.2*nys))    ; 51 ;21 
    sumarr2 = fltarr(ntest)
    testarr2 = nys - findgen(ntest)
    ;testarr2 = findgen(ntest)-(ntest-1)*0.5+nys
    for i=0,ntest-1 do begin
      factor,testarr2[i],p,nn,/quiet
      sumarr2[i] = total(p*nn)
    end
    best2 = first_el(minloc(sumarr2))
    nys = testarr2[best2]
  end

; Getting NX/NY from input FFT arrays
endif else begin

  if n_elements(fft1) gt 0 then begin
    sz = size(fft1)
    nxs = sz[1]
    nys = sz[2]
  endif
  if n_elements(fft2) gt 0 then begin
    sz = size(fft2)
    nxs = sz[1]
    nys = sz[2]
  endif
  if n_elements(fftp) gt 0 then begin
    sz = size(fftp)
    nxs = sz[1]
    nys = sz[2]
  endif

endelse

; Getting stars inside this box
gd1b = where(xx1s le (nxs-1) and yy1s le (nys-1),ngd1b)
gd2b = where(xx2s le (nxs-1) and yy2s le (nys-1),ngd2b)
if ngd1b eq 0 or ngd2b eq 0 then begin
  if not keyword_set(silent) then $
    print,'MATCHSTARS_XCORR error'
  xshift = 999999.
  yshift = 999999.
  angle = 999999.
  xcorr = -1
  bestcorr = -1
  matchnum = -1
  nsig = -1
  return
endif
xx1s = xx1s[gd1b]
yy1s = yy1s[gd1b]
xx2s = xx2s[gd2b]
yy2s = yy2s[gd2b]


; THINGS TO DO:
; -TAKE STARS AROUND THE CENTER!!!!
; -Is fwhm=5 good enough for the scaled image?
; -Is osamp=4 a good number?
; -Is smooth=5 a good number?
; -Make Xpeak analysis more robust, 
;   better estimate of size for subxcorr image
; -Is there a bias in the rotation angle because
;   of the fwhm and smoothing?  Should fwhm and
;   smooth be changed to improve the rotation estimate?
;print,'THINGS TO DO!'
;stop


;#################################
; CROSS-CORRELATION
;#################################

; Making the "fake" images of stars
if n_elements(fft1) eq 0 then begin
  newim1 = fltarr(nxs,nys)
  ; Put in a gradient
  if keyword_set(gradient) then $
  newim1[round(xx1s),round(yy1s)] = scale(yy1s,[0,max(yy1s)],[0.5,1.0]) else $
  newim1[round(xx1s),round(yy1s)]++
  ;newim1[round(xx1s),round(yy1s)] = 1.0
endif
if n_elements(fft2) eq 0 then begin
  newim2 = fltarr(nxs,nys)
  newim2[round(xx2s),round(yy2s)]++
  ;newim2[round(xx2s),round(yy2s)] = 1.0
endif
numpix = nxs*nys

; NOTE on Normalization
; IDL divides the Forward transform by Numpix and multiplies the 
; backward transform by Numpix.  Need to remove this to normalization
; so we get the right numbers for the cross-correlation
;  forwardfft = FFT(im,-1) * numpix    ; correct for normalization
;  backwardfft = FFT(im,1) / numpix    ; correct for normalization

; Convolving images with PSF in fourier space
if n_elements(fftp) eq 0 then begin
  if n_elements(fwhm0) eq 0 then fwhm=5 else fwhm=fwhm0
  save_except = !EXCEPT & !EXCEPT=0

  psf = psf_gaussian(npixel=[nxs,nys],fwhm=[1,1]*fwhm,/norm)  ; make a PSF image
  fftp = FFT(psf,-1) * numpix    ; correct for normalization

  dum = check_math()
  !EXCEPT = save_except
endif

; Fourier transforming the images
if n_elements(fft1) eq 0 then begin
  fft1 = FFT(newim1,-1) * numpix  ; correct for normalization
  fft1 = fft1 * fftp
endif
if n_elements(fft2) eq 0 then begin
  fft2 = FFT(newim2,-1) * numpix         ; correct for normalization
  fft2 = fft2 * fftp
endif

; Do the Cross-correlation
xcorr = REAL_PART( FFT( fft1 * CONJ(fft2) ,1) ) / numpix  ; correct for normalization
xcorr_orig = xcorr

; Shift to center and smooth
xhalfs = nxs/2
yhalfs = nys/2
;xhalfs = long(0.5*nxs)
;if odd(nxs) eq 1 then xhalfs = floor(0.5*nxs)
;yhalfs = long(0.5*nys)
;if odd(nys) eq 1 then yhalfs = floor(0.5*nys)

; Shift the Xcorr array
xcorr = shift(xcorr,xhalfs,yhalfs)
; Smooth the Xcorr array
if keyword_set(smooth0) then begin
  smlen = smooth0
  if smlen eq 1 then smlen=5
  xcorr = smooth(xcorr,smlen,/edge_truncate)
endif


; Remove a smoothed background
; This removed bad structure in the background
bsmooth = nxs*0.10 > nys*0.10
backgd = SMOOTH(xcorr,bsmooth,/edge_truncate)
xcorr = xcorr - backgd

;; Maybe cross-correlation VERY smoothed versions
;; of newim1 and newim2 and then subtract that
;; from the "normal" cross-correlation image
;; Does not seem to do as good of a job as
;; removing a smooth background from xcorr.
; bsmooth = nxs*0.10 > nys*0.10
; smim1 = SMOOTH(newim1,bsmooth,/edge)
; smim2 = SMOOTH(newim2,bsmooth,/edge)
; crosscorr,smim1,smim2,xsh,ysh,bestc,smxcorr
; xcorr2 = xcorr - smxcorr

; Finding the best shift
;------------------------
; Constraining the shift
if keyword_set(maxshift) then begin
  xcorr2 = xcorr
  maxsh = round(abs(maxshift)/xyscale)
  szx = size(xcorr)
  if maxsh lt xhalfs then begin
    xcorr2[0:xhalfs-maxsh-1,*]=-1000
    xcorr2[xhalfs+maxsh:szx[1]-1,*]=-1000
  endif
  if maxsh lt yhalfs then begin
    xcorr2[*,0:yhalfs-maxsh-1]=-1000
    xcorr2[*,yhalfs+maxsh:szx[2]-1]=-1000
  endif
  bestcorr = max(xcorr2)
  bestind = first_el(maxloc(xcorr2)) 
  bestind2 = ARRAY_INDICES(xcorr2,bestind)
  xoffs = bestind2[0]
  yoffs = bestind2[1]
  xshift = xoffs-xhalfs
  yshift = yoffs-yhalfs
endif else begin
  bestcorr = max(xcorr)
  bestind = first_el(maxloc(xcorr)) 
  bestind2 = ARRAY_INDICES(xcorr,bestind)
  xoffs = bestind2[0]
  yoffs = bestind2[1]
  xshift = xoffs-xhalfs
  yshift = yoffs-yhalfs
endelse
xshift_orig = xshift
yshift_orig = yshift

;##############################################################
; Find rotation angle and XSHIFT/YSHIFT as center of XCORR peak
;##############################################################
; Use partial sums to find the widths of the xcorr peak
; this tells you the rotation angle
;
; if the heights of the peaks have a gradient in y
; and there is rotation then one side of xcorr will be
; higher than the other (in x).  Compare the mean of
; the left half to the mean of the right half.

if keyword_set(extra) then begin

  ; Get narrow profiles across the peak
  size = round(0.1*nxs) > 40.
  ; Narrow X profile
  ;-----------------
  xlo1 = (xoffs-size) > 0
  xhi1 = (xoffs+size) < (nxs-1)
  ylo1 = (yoffs-1) > 0
  yhi1 = (yoffs+1) < (nys-1)
  nysum = (yhi1-ylo1)+1
  xprofile = total(xcorr[xlo1:xhi1,ylo1:yhi1],2)/nysum
  xprofile = xprofile-min(xprofile)       ; normalize
  xprofile = xprofile/max(xprofile)
  nxprofile = n_elements(xprofile)

  ; Where does the slope become negative left of the peak
  ; and positive right of the peak
  dxprofile = SLOPE(xprofile,/acc)
  maxxind = first_el(maxloc(xprofile))
  xlftind = first_el(where(dxprofile[0:(maxxind-1)>0] le 0.0),/last)
  if xlftind eq -1 then xlftind = 0
  xrtend = (maxxind+1) < (nxprofile-1)
  xrtind = first_el(where(dxprofile[xrtend:*] ge 0.0)) + xrtend
  if xrtind eq -1 then xrtind = n_elements(xprofile)-1
  xwidth1 = (xrtind-xlftind)+1    ; this is going to be larger than the width

  ; Narrow Y profile
  ;-----------------
  xlo1 = (xoffs-1) > 0
  xhi1 = (xoffs+1) < (nxs-1)
  ylo1 = (yoffs-size) > 0
  yhi1 = (yoffs+size) < (nys-1)
  nxsum = (yhi1-ylo1)+1
  yprofile = total(xcorr[xlo1:xhi1,ylo1:yhi1],1)/nxsum
  yprofile = yprofile-min(yprofile)       ; normalize
  yprofile = yprofile/max(yprofile)
  nyprofile = n_elements(yprofile)

  ; Where does the slope become negative left of the peak
  ; and positive right of the peak
  dyprofile = SLOPE(yprofile,/acc)
  maxyind = first_el(maxloc(yprofile))
  ylftind = first_el(where(dyprofile[0:(maxyind-1)>0] le 0.0),/last)
  if ylftind eq -1 then ylftind = 0
  yrtend = (maxyind+1) < (nyprofile-1)
  yrtind = first_el(where(dyprofile[yrtend:*] ge 0.0)) + yrtend
  if yrtind eq -1 then yrtind = n_elements(yprofile)-1
  ywidth1 = (yrtind-ylftind)+1    ; this is going to be larger than the width

  ; Getting an image of JUST the PEAK
  ;----------------------------------
  ; Looking at a broader range
  ;size = round(0.1*nxs)  ; 200
  xsize = round(xwidth1*0.5) > 5
  ysize = round(ywidth1*0.5) > 5
  xlo = (xoffs - xsize) > 0
  xhi = (xoffs + xsize) < (nxs-1)
  ylo = (yoffs - ysize) > 0
  yhi = (yoffs + ysize) < (nys-1)
  ;subxcorr = smxcorr2[xlo:xhi,ylo:yhi]
  subxcorr = xcorr[xlo:xhi,ylo:yhi]
  szsub = size(subxcorr)
  nxsub = szsub[1]
  nysub = szsub[2]
  med = median([subxcorr],/even)
  std = mad(subxcorr)
  subxcorr = subxcorr - med
  highind = first_el(maxloc(subxcorr))
  highind2 = array_indices(subxcorr,highind)
  xhighind = highind2[0]
  yhighind = highind2[1]

  ; Using CONTOUR
  ;-----------------
  minsub = min(subxcorr)
  maxsub = max(subxcorr)
  rsub = maxsub-minsub
  levels = scale_vector(findgen(10), minsub+0.1*rsub, maxsub-0.1*rsub)
  contour,subxcorr,levels=levels,/path_data_coords,path_info=path_info,$
          path_xy=path_xy
  nlevels = n_elements(path_info)   ; there might be less than were requested

  ; Looping through the levels 
  xwidarr = fltarr(nlevels)-1
  ywidarr = fltarr(nlevels)-1
  radarr = fltarr(nlevels)-1
  xloarr = fltarr(nlevels)-1
  xhiarr = fltarr(nlevels)-1
  yloarr = fltarr(nlevels)-1
  yhiarr = fltarr(nlevels)-1

  for i=0,nlevels-1 do begin
    lo = path_info[i].offset
    hi = lo + path_info[i].n - 1
    x = reform(path_xy[0,lo:hi])
    y = reform(path_xy[1,lo:hi])

    xloarr[i] = min(x)
    xhiarr[i] = max(x)
    yloarr[i] = min(y)
    yhiarr[i] = max(y)

    ; Does this include the peak?
    ROI_CUT,x,y,[xhighind],[yhighind],ind,cutind,fac=1,/silent

    ; Peak inside, calculate widths
    if n_elements(cutind) eq 0 then cutind=-1
    if cutind[0] ne -1 then begin
      dist = sqrt( (x-xhighind[0])^2.0 + (y-yhighind)^2.0 )
      xwidarr[i] = range(x)
      ywidarr[i] = range(y)
      radarr[i] = median([dist],/even)
      ;stop
    endif

  end

  ; We have some good levels
  gdlev = where(xwidarr ne -1,ngdlev)
  if ngdlev gt 0 then begin
    xwidarr2 = xwidarr[gdlev]
    ywidarr2 = ywidarr[gdlev]
    radarr2 = radarr[gdlev]
    levels2 = levels[gdlev]
    xloarr2 = xloarr[gdlev]
    xhiarr2 = xhiarr[gdlev]
    yloarr2 = yloarr[gdlev]
    yhiarr2 = yhiarr[gdlev]

    ; Spline the xwidth/ywidth at 1/2 max height
    halfheight = max(subxcorr)*0.5
    xwidth = ( interpol(xwidarr2, levels2, max(subxcorr)*0.5) > 1 ) < nxsub
    ywidth = ( interpol(ywidarr2, levels2, max(subxcorr)*0.5) > 1 ) < nysub
    rwidth = ( interpol(radarr2, levels2, max(subxcorr)*0.5) > 1 ) < nxsub

    ; Spline the positions of the Half max region
    xlo_halfbox = ( ceil(interpol(xloarr2, levels2, max(subxcorr)*0.5)) > 0 ) < (nxsub-1)
    xhi_halfbox = ( floor(interpol(xhiarr2, levels2, max(subxcorr)*0.5)) > 0 ) < (nxsub-1)
    ylo_halfbox = ( ceil(interpol(yloarr2, levels2, max(subxcorr)*0.5)) > 0 ) < (nysub-1)
    yhi_halfbox = ( floor(interpol(yhiarr2, levels2, max(subxcorr)*0.5)) > 0 ) < (nysub-1)

    ; Estimate xshift/yshift from half max region
    ; This can give bad results if the peak is on the edge
    xhalfarr = (xhiarr2-xloarr2)*0.5 + xloarr2
    yhalfarr = (yhiarr2-yloarr2)*0.5 + yloarr2
    xmid_halfbox = interpol(xhalfarr, levels2, max(subxcorr)*0.5)
    ymid_halfbox = interpol(yhalfarr, levels2, max(subxcorr)*0.5)
    xshift2s = xmid_halfbox + xlo - xhalfs
    yshift2s = ymid_halfbox + ylo - yhalfs

    ;setdisp
    ;displayc,subxcorr
    ;contour,subxcorr,nlevels=20,/over
    ;print,xwidth,ywidth,rwidth*2.0
    ;print,xlo_halfbox, xhi_halfbox, ylo_halfbox, yhi_halfbox
    ;stop


  ; Contour failed, use MARGINAL SUMS
  ; the marginal sums method can be skewed if the peak has an elongated
  ; shape at an angle
  endif else begin

    ; X marginal sum
    xpsum = total(subxcorr,2)  ; sum in y
    xpsum = xpsum - min(xpsum)
    xpsum = xpsum/max(xpsum)   ; normalize
    xstd = mad(xpsum)

    ; spline   
    osamp = 10
    nx2 = n_elements(xpsum)
    x1 = findgen(nx2) + xlo
    x2 = scale_vector(findgen(nx2*osamp),min(x1),max(x1))
    xpsum2 = spline(x1,xpsum,x2)
    maxxind = first_el(maxloc(xpsum2))
    xlftind = first_el(where(xpsum2[0:maxxind] ge 0.5))
    if xlftind eq -1 then xlftind = 0
    xrtind = first_el(where(xpsum2[maxxind:*] le 0.5))
    if xrtind eq -1 then xrtind = (nx2*osamp-1) else xrtind=xrtind+maxxind
    xwidth = x2[xrtind]-x2[xlftind]
    ; get a shift from the middle of the peak, not necessarily the highest point
    xshift2s = mean(x2[[xlftind,xrtind]]) - xhalfs

    ; Y marginal sum
    ypsum = total(subxcorr,1)  ; sum in x
    ypsum = ypsum - min(ypsum)
    ypsum = ypsum/max(ypsum)   ; normalize
    ystd = mad(ypsum)

    ; spline
    ny2 = n_elements(ypsum)
    y1 = findgen(ny2) + ylo
    y2 = scale_vector(findgen(ny2*osamp),min(y1),max(y1))
    ypsum2 = spline(y1,ypsum,y2)
    maxyind = first_el(maxloc(ypsum2))
    ylftind = first_el(where(ypsum2[0:maxyind] ge 0.5)) 
    if ylftind eq -1 then ylftind = 0
    yrtind = first_el(where(ypsum2[maxyind:*]  le 0.5))
    if yrtind eq -1 then yrtind = (ny2*osamp-1) else yrtind=yrtind+maxyind
    ywidth = y2[yrtind]-y2[ylftind]
    ; get a shift from the middle of the peak, not necessarily the highest point
    yshift2s = mean(y2[[ylftind,yrtind]]) - yhalfs

    ; Positions of the half max region
    szsub = size(subxcorr)
    xlo_halfbox = ceil(xlftind/osamp) > 0
    xhi_halfbox = floor(xrtind/osamp) < (szsub[1]-1)
    ylo_halfbox = ceil(ylftind/osamp) > 0
    yhi_halfbox = floor(yrtind/osamp) < (szsub[2]-1)

  endelse  ; marginal sums


  ; Estimate rotation angle, from the width of the peak
  ;-----------------------------------------------------
  ang1 = xwidth/nys
  ang2 = ywidth/nxs
  angamp = mean([ang1,ang2])*!radeg    ; amplitude of rotation

  ; Get sign of rotation
  peakim = subxcorr[xlo_halfbox:xhi_halfbox, ylo_halfbox:yhi_halfbox]
  peakim = peakim - min(subxcorr)
  peakim = peakim/max(peakim)

  szpeakim = size(peakim)
  nxpeakim = szpeakim[1]
  xmid = round(nxpeakim*0.5)
  lftmn = mean(peakim[0:xmid,*])
  rtmn = mean(peakim[xmid:*,*])
  sign = 1.0
  if (lftmn gt rtmn) then sign=-1.0
  angle = sign*angamp
  if not keyword_set(gradient) then angle=angamp   ; no sign information


; Don't do any extra fancy work
;------------------------------
endif else begin

  xshift2s = xshift
  yshift2s = yshift
  angle = 0.0

endelse

;-----------------------------------------
; NOW BACK TO THE **ORIGINAL** SCALE
;-----------------------------------------
xshift = xshift2s*xyscale
yshift = yshift2s*xyscale
xhalf = xhalfs*xyscale
yhalf = yhalfs*xyscale


; Estimate of how many matched
if n_elements(fwhm) eq 0 then fwhm=5
if n_elements(psf) eq 0 then $
tpsf = psf_gaussian(npixel=[11,11],fwhm=[1,1]*fwhm,/norm) else $
 tpsf = psf
; Smooth the Xcorr array
if keyword_set(smooth0) then begin
  smlen = smooth0   
  if smlen eq 1 then smlen=5
  tpsf = smooth(tpsf,smlen,/edge_truncate)
endif
save_except=!EXCEPT & !EXCEPT=0
onematch = total(tpsf*tpsf)    ; this gives math errors sometimes
dum = check_math() & !EXCEPT = save_except
matchnum = (bestcorr-median([xcorr],/even))/onematch

; How significant is this peak
nsig = (bestcorr-median([xcorr],/even))/mad(xcorr)

; How significant does it need to be for the number of pixels
;http://en.wikipedia.org/wiki/Normal_distribution
;prof=erf(n/sqrt(2)), where n is number of sigma.
;this is the probability of being within +/-n sigma.
;we need the probability of being N sigma to one side.
;so we want (1.0-(erf(n/sqrt(2)))*0.5

nsigma = scale_vector(dindgen(100),3.0,8.3)
prob =0.5d0*(1.0d0-erf(nsigma/sqrt(2.0d0))) 
npix = nxs*nys
tot = npix*prob
bestind = min(where(tot lt 1.0,ngd))
if bestind gt -1 then nsiglim = nsigma[bestind[0]] else nsiglim=5.0

;print,nsiglim
;stop

; CORRECT ANGLE FOR FWHM EFFECTS
; Doesn't seem to make a big effect
;stop

if keyword_set(stp) then stop

end



;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


pro matchstars_linefit,par0,xx1m0,yy1m0,xx2m0,yy2m0,fpar,rms1,plotresid=plotresid,stp=stp,$
                       silent=silent


; This removes any leftover rotation by fitting lines to the x/y residuals
; iteratively

; Error Handling
;------------------
; Establish error handler. When errors occur, the index of the  
; error is returned in the variable Error_status:  
CATCH, Error_status 

;This statement begins the error handler:  
if (Error_status ne 0) then begin 
   fpar = 999999.
   rms1 = 999999.
   PHOTRED_ERRORMSG
   CATCH, /CANCEL 
   return
endif


xx1m = xx1m0
yy1m = yy1m0
xx2m = xx2m0
yy2m = yy2m0

out = trans_coo(xx2m,yy2m,par0)
xx2c = reform(out[0,*])
yy2c = reform(out[1,*])

x2med = median([xx2m],/even)
y2med = median([yy2m],/even)

; Nothing to fit
xrms = mad(xx1m-xx2m)
yrms = mad(yy1m-yy2m)
if xrms lt 0.01 and yrms lt 0.01 then begin
  rms1 = 0.01
  fpar = [mean(xx1m-xx2m), mean(yy1m-yy2m), 1.0, 0.0, 0.0, 1.0]
  return
endif

count = 0
rms1 = 99.0
flag = 0
;par = fltarr(6)+99999.
par = par0

; Loop until converged
WHILE (flag eq 0) do begin

  oldrms = rms1
  oldpar = par

  out = trans_coo(xx2m,yy2m,par)
  xx2c = reform(out[0,*])
  yy2c = reform(out[1,*])

  ; Do robust line fitting to the residuals
  ; deltaX vs. Y
  ; deltaY vs. X
  coef1x = ROBUST_LINEFIT(xx2m,xx1m-xx2c,yfit1,sig1,coef_sig1)
  coef1y = ROBUST_LINEFIT(yy2m,xx1m-xx2c,yfit2,sig2,coef_sig2)
  coef2x = ROBUST_LINEFIT(xx2m,yy1m-yy2c,yfit1,sig1,coef_sig1)
  coef2y = ROBUST_LINEFIT(yy2m,yy1m-yy2c,yfit2,sig2,coef_sig2)

  if n_elements(coef1x) lt 2 or n_elements(coef1y) lt 2 or $
     n_elements(coef2x) lt 2 or n_elements(coef2y) lt 2 then begin
    if not keyword_set(silent) then $
      print,'MATCHSTARS_LINEFIT problem'
    undefine,fpar
    rms1 = 999999.
    return
  endif


  ; Making transformation array
  ; par = [A, B, C, D, E, F]
  ; x(1) = A + C*x(n) + E*y(n)
  ; y(1) = B + D*x(n) + F*y(n)
  ; Adding to the linear coefficients
  par[2] = par[2] + coef1x[1]                      ; x coef for x
  par[3] = par[3] + coef2x[1]                      ; x coef for y
  par[4] = par[4] + coef1y[1]                      ; y coef for x
  par[5] = par[5] + coef2y[1]                      ; y coef for y
  ;par = [coef1[0], coef2[0], cos(coef2[1]), coef2[1], coef1[1], cos(coef2[1])]

  ; Adding constant offsets
  ; changing the linear coefficients will also produce an offset
  xoff = median([xx1m-xx2c],/even)
  yoff = median([yy1m-yy2c],/even)
  par[0] = par[0] + xoff - ( x2med*coef1x[1] + y2med*coef1y[1] )
  par[1] = par[1] + yoff - ( x2med*coef2x[1] + y2med*coef2y[1] )

  ;plot,yy1m,xx1m-xx2c,ps=3
  ;oplot,xx1m,yy1m-yy2c,ps=3,co=250
  ;stop

  out = trans_coo(xx2m,yy2m,par)
  xx2c = reform(out[0,*])
  yy2c = reform(out[1,*])

  ; Match them
  dcr = (3.0*rms1 > 5.0) < 15.
  SRCMATCH,xx1m,yy1m,xx2c,yy2c,dcr,ind1b,ind2b,count=nind1b
  if nind1b lt 2 then begin
    undefine,fpar
    if not keyword_set(silent) then $
      print,'NOT ENOUGH MATCHES'
    return
  endif
  xx1m = xx1m[ind1b]
  yy1m = yy1m[ind1b]
  xx2m = xx2m[ind2b]
  yy2m = yy2m[ind2b]

  ; Plot residuals
  ;plotresid=1
  if keyword_set(plotresid) then begin
    plot,yy1m,xx1m-xx2c[ind2b],ps=3
    oplot,xx1m,yy1m-yy2c[ind2b],ps=3,co=250
    stop
  endif

  ; Remove outliers
  xresid = xx1m-xx2c[ind2b]
  yresid = yy1m-yy2c[ind2b]
  xmed = median([xresid],/even)
  ymed = median([yresid],/even)
  xrms1 = mad(xresid)
  yrms1 = mad(yresid)
  diff = sqrt(xresid^2.0 + yresid^2.0)
  rms1 = sqrt(mean(diff^2.0))
  ;rms1 = sqrt(xrms1^2.0 + yrms1^2.0)
  xlim = (3.0*xrms1) > 0.1
  ylim = (3.0*yrms1) > 0.1
  gd = where(abs(xresid-xmed) lt xlim and abs(yresid-ymed) lt ylim,ngd)
  xx1m = xx1m[gd]
  yy1m = yy1m[gd]
  xx2m = xx2m[gd]
  yy2m = yy2m[gd]

  ; Check if we've converged
  dpar = oldpar-par
  dpar2 = abs(dpar)/abs(par)*100.
  drms = oldrms - rms1
  drms2 = drms/rms1*100.
  if (((drms2 lt 0.5) or (drms lt 0.01)) and (max(dpar2) lt 1.0)) then flag=1
  if (max(dpar2) lt 0.1) then flag=1
  if (count ge 20) then flag=1


  if keyword_set(plotresid) then begin
    print,'count = ',count
    print,'rms=',rms1
    print,'drms=',drms,' drms2=',drms2
    print,'coef1x=',coef1x
    print,'coef1y=',coef1y
    print,'coef2x=',coef2x
    print,'coef2y=',coef2y
    print,'par=',par
    print,'dpar2=',dpar2
  endif

  count++

  ;stop

ENDWHILE

fpar = par

if keyword_set(stp) then stop

end


;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


PRO matchstars,xinp1,yinp1,xinp2,yinp2,ind1,ind2,trans,plotresid=plotresid,$
               norot=norot,maxshift=maxshift,silent=silent,count=count,stp=stp,rms=rms

COMMON xcorr, factorsum, factornum

undefine,rms,trans
rms = 999999.

t0 = systime(1)
undefine,ind1,ind2,trans
count = 0

; Error Handling
;------------------
; Establish error handler. When errors occur, the index of the  
; error is returned in the variable Error_status:  
CATCH, Error_status 

;This statement begins the error handler:  
if (Error_status ne 0) then begin 
   ind1 = -1
   ind2 = -1
   trans = -1
   count = -1          ; There was an error
   PHOTRED_ERRORMSG
   CATCH, /CANCEL 
   return
endif


; There can be problems sometimes if any of the star lists
; are empty on the corners.

nxinp1 = n_elements(xinp1)
nyinp1 = n_elements(yinp1)
nxinp2 = n_elements(xinp2)
nyinp2 = n_elements(yinp2)

; Not enough inputs
if nxinp1 eq 0 or nyinp1 eq 0 or nxinp2 eq 0 or nyinp2 eq 0 then begin
  if not keyword_set(silent) then $
  print,'Syntax - matchstars,xinp1,yinp1,xinp2,yinp2,ind1,ind2,trans,plotresid=plotresid,'
  print,'                    norot=norot,maxshift=maxshift,stp=stp'
  return
endif

; Array sizes don't match
if nxinp1 ne nyinp1 then begin
  print,'xinp1 and yinp1 must have same number of elements'
  return
endif
if nxinp2 ne nyinp2 then begin
  print,'xinp2 and yinp2 must have same number of elements'
  return
endif

; Restoring factorsum for Xcorr
idldir = file_search('~/idl/',/fully_qualify_path)
test = file_test(idldir+'/factorsum.dat')
if test eq 1 then begin
  restore,idldir+'/factorsum.dat'
endif else begin
  print,'factorsum.dat NOT FOUND.  PLEASE copy factorsum.dat to ~/idl/'
endelse

; Getting positive coordinates
minx = floor(min([xinp1,xinp2]))
miny = floor(min([yinp1,yinp2]))
maxx = ceil(max([xinp1,xinp2]))
maxy = ceil(max([yinp1,yinp2]))

; Applying an offset so that all of the coordinates are positive
xx1 = xinp1 - minx
yy1 = yinp1 - miny
xx2 = xinp2 - minx
yy2 = yinp2 - miny

maxx2 = maxx-minx
maxy2 = maxy-miny

;######################################
; CROSS-CORRELATION
;######################################

;-----------------------------------------------------------------------------------
; Checking SMALL Rotations
;-----------------------------------------------------------------------------------


xflip = 0
if not keyword_set(silent) then $
  print,'CHECKING for SMALL Rotations at default orientation'

; Try a quick cross-correlation to check if it's a good match
; This can pick up good fits for rotations of up to +/-4 degrees
xyscale = 4
fwhm = 5
smooth = 5
MATCHSTARS_XCORR,xx1,yy1,xx2,yy2,xshift1,yshift1,angle1,bestcorr,xcorr,smooth=smooth,$
                 xyscale=xyscale,fwhm=fwhm,maxshift=maxshift,nsig=nsig1,matchnum=matchnum1,$
                 /extra,silent=silent
nsig = nsig1

; No rotation
if keyword_set(norot) and nsig gt 5 then begin
  xshift = xshift1
  yshift = yshift1
  angle = 0.0
  matchnum = matchnum1
endif


; Good SMALL Rotation Fit
; GETTING the sign of the angle by comparing XSHIFT from NORTHERN and SOUTHERN stars
;-----------------------------------------------------------------------------------

if (nsig gt 5.0) and not keyword_set(norot) then begin

  ; Shift the stars and just get those that overlap
  xx2sh = xx2 + xshift1
  yy2sh = yy2 + yshift1
  txmin = min(xx1) > min(xx2sh)
  txmax = max(xx1) < max(xx2sh)
  tymin = min(yy1) > min(yy2sh)
  tymax = max(yy1) < max(yy2sh)
  tyhalf = mean([tymin,tymax])

  gd1 = where(xx1 ge txmin and xx1 le txmax and yy1 ge tymin and yy1 le tymax,ngd1)
  gd2 = where(xx2sh ge txmin and xx2sh le txmax and yy2sh ge tymin and yy2sh le tymax,ngd2)

  txx1 = xx1[gd1]
  tyy1 = yy1[gd1]
  txx2 = xx2sh[gd2]
  tyy2 = yy2sh[gd2]

  xyscale = 4
  smooth = 5
  fwhm = 5
              
  ; Find the NORTH shift
  ;xhalf = round(0.5*maxx2)
  ;yhalf = round(0.5*maxy2)
  Ngd1 = where(tyy1 ge tyhalf,nNgd1)   
  Ngd2 = where(tyy2 ge tyhalf,nNgd2)
  xx1s = txx1[Ngd1]
  yy1s = tyy1[Ngd1]  ;-tyhalf
  xx2s = txx2[Ngd2]
  yy2s = tyy2[Ngd2]  ;-tyhalf
  medyN = median([tyy2[Ngd2]],/even)

  MATCHSTARS_XCORR,xx1s,yy1s,xx2s,yy2s,xshiftN,yshiftN,angleN,bestcorrN,xcorrN,silent=silent,$
                   smooth=smooth,xyscale=xyscale,fwhm=fwhm,nsig=nsigN,matchnum=matchnumN,/extra

  ; Find the SOUTH shift
  Sgd1 = where(tyy1 le tyhalf,nSgd1)
  Sgd2 = where(tyy2 le tyhalf,nSgd2)
  xx1s = txx1[Sgd1]
  yy1s = tyy1[Sgd1]
  xx2s = txx2[Sgd2]
  yy2s = tyy2[Sgd2]
  medyS = median([tyy2[Sgd2]],/even)

  MATCHSTARS_XCORR,xx1s,yy1s,xx2s,yy2s,xshiftS,yshiftS,angleS,bestcorrS,xcorrS,silent=silent,$
                   smooth=smooth,xyscale=xyscale,fwhm=fwhm,nsig=nsigS,matchnum=matchnumS,/extra

  ; Good fits
  ;if nsigS gt 5.0 and nsigN gt 5.0 then begin
  if nsigS gt 3.5 and nsigN gt 3.5 then begin
    
    ; Average Xshift/Yshift
    ;xshift = mean([xshiftN,xshiftS])
    ;xshift = mean([xshift,xshift1])
    ;yshift = mean([yshiftN,yshiftS])
    ;yshift = mean([yshift,yshift1])
    xshift = xshift1
    yshift = yshift1
    ; Calculating angle from the shifts
    angle = (xshiftN-xshiftS)/(medyN-medyS)*!radeg

    ; Number of matches
    matchnum = matchnumN + matchnumS

    ; Is the angle too large
    if (abs(angle) gt 2.0) then begin
      if not keyword_set(silent) then $
        print,'Angle is too large.  Doing larger search.'
      nsig = 0.
    endif

  ; Bad fits
  endif else begin

    xshift = xshift1
    yshift = yshift1
    angle = 0.0
    matchnum = matchnum1
    nsig = 0

    ;print,'Xshift = ',strtrim(xshift,2)
    ;print,'Yshift = ',strtrim(yshift,2)
    ;print,'Angle unknown'


  endelse


  ; Printing info if we have a good fit
  if (nsig gt 5.0) and not keyword_set(silent) then begin
    print,'Rough cross-correlation shifts and rotation angle:'
    print,'Xshift = ',strtrim(xshift,2),' pixels'
    print,'Yshift = ',strtrim(yshift,2),' pixels'
    print,'Angle  = ',strtrim(angle,2),' degrees'
    print,'Nsig = ',strtrim(nsig,2)
    print,'Rough number of Xcorr matches = ',strtrim(matchnum,2)
    print,''
  endif

  ;stop

end

;stop

;----------------------------------------------------------------------
; Checking INTERMEDIATE Rotations
; Fit is bad, check rotations of up to +/-5 deg around the 8 orientations
;----------------------------------------------------------------------
xhalf = round(0.5*maxx2)
yhalf = round(0.5*maxy2)

if (nsig lt 5.) then begin

  if not keyword_set(silent) then $
    print,'CHECKING all 8 orientations with intermediate rotations'

  ; Go through the orientations
  rotarr1 = [0.0, 90.0, 180.0, 270.0, 0.0, 90.0, 180.0, 270.0]
  nsigarr1 = fltarr(8)
  bestrot1 = fltarr(8)
  xshiftarr1 = fltarr(8)
  yshiftarr1 = fltarr(8)
  nmatcharr1 = fltarr(8)

  for i=0,7 do begin

    ; Check the rotations
    nrot = 11
    rotarr2 = findgen(nrot)*1.0-5.0 + rotarr1[i] 
    nsigarr2 = fltarr(nrot)
    xshiftarr2 = fltarr(nrot)
    yshiftarr2 = fltarr(nrot)
    nmatcharr2 = fltarr(nrot)

    for j=0,nrot-1 do begin

      ; Rotate the coordinates
      rxx2 =  (xx2-xhalf)*cos(rotarr2[j]/!radeg) + (yy2-yhalf)*sin(rotarr2[j]/!radeg) + xhalf
      ryy2 = -(xx2-xhalf)*sin(rotarr2[j]/!radeg) + (yy2-yhalf)*cos(rotarr2[j]/!radeg) + yhalf
      if (i ge 4) then rxx2 = maxx2-rxx2  ; flip

      ; Only use a central box
      ;box = 0.5*(maxx < maxy)
      ;gd1 = where(xx1 ge (xhalf-box*0.5) and xx1 le (xhalf+box*0.5) and $
      ;            yy1 ge (yhalf-box*0.5) and yy1 le (yhalf+box*0.5),ngd1)
      ;gd2 = where(rxx2 ge (xhalf-box*0.5) and rxx2 le (xhalf+box*0.5) and $
      ;            ryy2 ge (yhalf-box*0.5) and ryy2 le (yhalf+box*0.5),ngd2)
      ;xx1s = xx1[gd1]-(xhalf-box*0.5)
      ;yy1s = yy1[gd1]-(yhalf-box*0.5)
      ;xx2s = rxx2[gd2]-(xhalf-box*0.5)
      ;yy2s = ryy2[gd2]-(yhalf-box*0.5)

      gd2 = where(rxx2 ge 0 and rxx2 le max(xx1) and ryy2 ge 0 and ryy2 le max(yy1),ngd2)
      xx1s = xx1
      yy1s = yy1
      xx2s = rxx2[gd2]
      yy2s = ryy2[gd2]


      MATCHSTARS_XCORR,xx1s,yy1s,xx2s,yy2s,xshift,yshift,angle,bestcorr,xcorr,silent=silent,$
                       smooth=5,xyscale=10,fwhm=3,matchnum=matchnum,nsig=nsig

      ; with xyscale=15 it is accurate to +/-10 deg
      ; only step by 5 deg
      ; stars are probably overlapping now

      nsigarr2[j] = nsig
      xshiftarr2[j] = xshift
      yshiftarr2[j] = yshift
      nmatcharr2[j] = matchnum

    endfor  ; rotations

    ; Find best rotation, spline
    rotarr2b = scale_vector(findgen(100),min(rotarr2),max(rotarr2))
    nsigarr2b = spline(rotarr2,nsigarr2,rotarr2b)
    xshiftarr2b = spline(rotarr2,xshiftarr2,rotarr2b)
    yshiftarr2b = spline(rotarr2,yshiftarr2,rotarr2b)
    nmatcharr2b = spline(rotarr2,nmatcharr2,rotarr2b)
    bestind = first_el(maxloc(nsigarr2b))

    nsigarr1[i] = nsigarr2b[bestind]
    bestrot1[i] = rotarr2b[bestind]
    ;if (i le 4) then bestrot1[i] = rotarr1[i] + rotarr2b[bestind]  ; normal
    ;if (i ge 4) then bestrot1[i] = rotarr1[i] - rotarr2b[bestind]  ; XFLIP
    xshiftarr1[i] = xshiftarr2b[bestind]
    yshiftarr1[i] = yshiftarr2b[bestind]
    nmatcharr1[i] = nmatcharr2b[bestind]

    ;stop

  endfor  ; 8 orientations

  ; Find best rotation
  bestind1 = first_el(maxloc(nsigarr1))
  nsig1 = nsigarr1[bestind1]
  if not keyword_set(silent) then $
    print,'Nsig = ',strtrim(nsig1,2) 


  ; We have a decent solution, Refining
  if (nsig1 gt 5.) then begin

    angle = bestrot1[bestind1]
    if bestind1 ge 4 then xflip=1 else xflip=0
    xshift1 = xshiftarr1[bestind1]
    yshift1 = yshiftarr1[bestind1]
    matchnum1 = nmatcharr1[bestind1]


    ; Do one more XCORR at higher resolution
    if not keyword_set(silent) then $
      print,'Refining the solution'

    ; Rotate the coordinates
    rxx2 =  (xx2-xhalf)*cos(angle/!radeg) + (yy2-yhalf)*sin(angle/!radeg) + xhalf
    ryy2 = -(xx2-xhalf)*sin(angle/!radeg) + (yy2-yhalf)*cos(angle/!radeg) + yhalf
    if xflip eq 1 then rxx2 = maxx2-rxx2  ; flip

    gd2 = where(rxx2 ge 0 and rxx2 le max(xx1) and ryy2 ge 0 and ryy2 le max(yy1),ngd2)
    xx2s = rxx2[gd2]
    yy2s = ryy2[gd2]
    
    MATCHSTARS_XCORR,xx1,yy1,xx2s,yy2s,xshift,yshift,angle0,bestcorr,xcorr,smooth=5,$
                 xyscale=4,fwhm=5,nsig=nsig,matchnum=matchnum,/extra,silent=silent
    ;             xyscale=4,fwhm=5,nsig=nsig,matchnum=matchnum,/extra,/gradient


    ; THESE angles are NOT that reliable
    ; Updating the angle
    ;angle_orig = angle
    ;if xflip eq 1 then angle = angle + angle0 ;- angle0
    ;if xflip eq 0 then angle = angle - angle0 ;+ angle0

    ; Printing info if we have a good fit
    if (nsig gt 5.0) and not keyword_set(silent) then begin
      print,'Rough cross-correlation shifts and rotation angle:'
      print,'Xshift = ',strtrim(xshift,2)
      print,'Yshift = ',strtrim(yshift,2)
      print,'Angle = ',strtrim(angle,2),' degrees'
      if xflip eq 1 then print,'XFLIP'
      print,'Nsig = ',strtrim(nsig,2)
      print,'Rough number of Xcorr matches = ',strtrim(matchnum,2)
      print,''
    endif else begin
      if not keyword_set(silent) then $
        print,'Fit not good'
    endelse

  ; No good solution
  endif else begin
    if not keyword_set(silent) then $
      print,'No good Intermediate Rotation solution'
  endelse


endif  ; not a good fit




;---------------------------------------------------------------------
; Check LARGE Rotations
;---------------------------------------------------------------------
; THIS WORKS EXCEPTIONALLY WELL NOW
; This gives *very* good transformations

if (nsig lt 5.) then begin

  if not keyword_set(silent) then $
    print,'CHECKING all possible rotations'

  dangle = 3 ;5.
  nrot = 360/dangle
  rotarr1 = findgen(nrot)*dangle
  nsigarr1 = fltarr(nrot)
  xshiftarr1 = fltarr(nrot)
  yshiftarr1 = fltarr(nrot)
  nmatcharr1 = fltarr(nrot)
  xflip = 0

  for i=0,nrot-1 do begin

    rxx2 =  (xx2-xhalf)*cos(rotarr1[i]/!radeg) + (yy2-yhalf)*sin(rotarr1[i]/!radeg) + xhalf
    ryy2 = -(xx2-xhalf)*sin(rotarr1[i]/!radeg) + (yy2-yhalf)*cos(rotarr1[i]/!radeg) + yhalf

    ;box = 0.5*(maxx < maxy)
    ;gd1 = where(xx1 ge (xhalf-box*0.5) and xx1 le (xhalf+box*0.5) and $
    ;            yy1 ge (yhalf-box*0.5) and yy1 le (yhalf+box*0.5),ngd1)
    ;gd2 = where(rxx2 ge (xhalf-box*0.5) and rxx2 le (xhalf+box*0.5) and $
    ;            ryy2 ge (yhalf-box*0.5) and ryy2 le (yhalf+box*0.5),ngd2)
    ;xx1s = xx1[gd1]-(xhalf-box*0.5)
    ;yy1s = yy1[gd1]-(yhalf-box*0.5)
    ;xx2s = rxx2[gd2]-(xhalf-box*0.5)
    ;yy2s = ryy2[gd2]-(yhalf-box*0.5)

    gd2 = where(rxx2 ge 0 and rxx2 le max(xx1) and ryy2 ge 0 and ryy2 le max(yy1),ngd2)
    xx1s = xx1
    yy1s = yy1
    xx2s = rxx2[gd2]
    yy2s = ryy2[gd2]

    MATCHSTARS_XCORR,xx1s,yy1s,xx2s,yy2s,xshift,yshift,angle,bestcorr,xcorr,silent=silent,$
                     smooth=5,xyscale=10,fwhm=3,matchnum=matchnum,nsig=nsig

    ; with xyscale=15 it is accurate to +/-10 deg
    ; only step by 5 deg
    ; stars are probably overlapping now

    nsigarr1[i] = nsig
    xshiftarr1[i] = xshift
    yshiftarr1[i] = yshift
    nmatcharr1[i] = matchnum

  endfor

  ; Find best rotation, spline
  rotarr1b = scale_vector(findgen(3600),min(rotarr1),max(rotarr1))
  nsigarr1b = spline(rotarr1,nsigarr1,rotarr1b)
  xshiftarr1b = spline(rotarr1,xshiftarr1,rotarr1b)
  yshiftarr1b = spline(rotarr1,yshiftarr1,rotarr1b)
  nmatcharr1b = spline(rotarr1,nmatcharr1,rotarr1b)
  bestind1 = first_el(maxloc(nsigarr1b))
  nsig1 = nsigarr1b[bestind1]
  angle1 = rotarr1b[bestind1]
  xshift1 = xshiftarr1b[bestind1]
  yshift1 = yshiftarr1b[bestind1]
  matchnum1 = nmatcharr1b[bestind1]

  ; We have a decent solution
  if (nsig1 gt 5.) then begin
    nsig = nsig1
    angle = angle1
    xshift = xshift1
    yshift = yshift1
    matchnum = matchnum1
    xflip = 0

    if not keyword_set(silent) then begin
      print,'Rough cross-correlation shifts and rotation angle:' 
      print,'Xshift = ',strtrim(xshift,2)
      print,'Yshift = ',strtrim(yshift,2)
      print,'Angle = ',strtrim(angle,2),' degrees'
      print,'Nsig = ',strtrim(nsig,2)
      print,'Rough number of Xcorr matches = ',strtrim(matchnum,2)
      print,''
    endif

  ; No good solution
  endif else begin
    if not keyword_set(silent) then $
      print,'No good LARGE Rotation solution'
  endelse




  ; Checking FLIPPED rotations
  ;--------------------------------------
  if (nsig1 lt 5.0) then begin

    if not keyword_set(silent) then $
      print,'Checking LARGE Rotations with one FLIPPED axis'

    ; Now try flipping the X coordinates
    dangle = 5 ; 3.  ; 5.
    nrot = 360/dangle
    rotarr2 = findgen(nrot)*dangle
    nsigarr2 = fltarr(nrot)
    xshiftarr2 = fltarr(nrot)
    yshiftarr2 = fltarr(nrot)
    nmatcharr2 = fltarr(nrot)

    for i=0,nrot-1 do begin
   
      rxx2 =  (xx2-xhalf)*cos(rotarr2[i]/!radeg) + (yy2-yhalf)*sin(rotarr2[i]/!radeg) + xhalf
      ryy2 = -(xx2-xhalf)*sin(rotarr2[i]/!radeg) + (yy2-yhalf)*cos(rotarr2[i]/!radeg) + yhalf
      rxx2 = maxx2-rxx2   ; flip x-coordinates
 
      ;box = 0.5*(maxx < maxy)
      ;gd1 = where(xx1 ge (xhalf-box*0.5) and xx1 le (xhalf+box*0.5) and $
      ;            yy1 ge (yhalf-box*0.5) and yy1 le (yhalf+box*0.5),ngd1)
      ;gd2 = where(rxx2 ge (xhalf-box*0.5) and rxx2 le (xhalf+box*0.5) and $
      ;            ryy2 ge (yhalf-box*0.5) and ryy2 le (yhalf+box*0.5),ngd2)
      ;xx1s = xx1[gd1]-(xhalf-box*0.5)
      ;yy1s = yy1[gd1]-(yhalf-box*0.5)
      ;xx2s = rxx2[gd2]-(xhalf-box*0.5)
      ;yy2s = ryy2[gd2]-(yhalf-box*0.5)

      gd2 = where(rxx2 ge 0 and rxx2 le max(xx1) and ryy2 ge 0 and ryy2 le max(yy1),ngd2)
      xx1s = xx1
      yy1s = yy1
      xx2s = rxx2[gd2]
      yy2s = ryy2[gd2]
  

      MATCHSTARS_XCORR,xx1s,yy1s,xx2s,yy2s,xshift,yshift,angle,bestcorr,xcorr,silent=silent,$
                       smooth=5,xyscale=10,fwhm=3,matchnum=matchnum,nsig=nsig


      ; with xyscale=15 it is accurate to +/-10 deg
      ; only step by 5 deg
      ; stars are probably overlapping now


      nsigarr2[i] = nsig
      xshiftarr2[i] = xshift
      yshiftarr2[i] = yshift
      nmatcharr2[i] = matchnum

    end

    ; Find best rotation, spline
    rotarr2b = scale_vector(findgen(3600),min(rotarr2),max(rotarr2))
    nsigarr2b = spline(rotarr2,nsigarr2,rotarr2b)
    xshiftarr2b = spline(rotarr2,xshiftarr2,rotarr2b)
    yshiftarr2b = spline(rotarr2,yshiftarr2,rotarr2b)
    nmatcharr2b = spline(rotarr2,nmatcharr2,rotarr2b)
    bestind2 = first_el(maxloc(nsigarr2b))
    nsig2 = nsigarr2b[bestind2]
    angle2 = rotarr2b[bestind2]
    xshift2 = xshiftarr2b[bestind2]
    yshift2 = yshiftarr2b[bestind2]
    matchnum2 = nmatcharr2b[bestind2]



    ; We have a decent solution
    if (nsig2 gt 5.) then begin
      nsig = nsig2
      angle = angle2
      xshift = xshift2
      yshift = yshift2
      matchnum = matchnum2
      xflip = 1

      if not keyword_set(silent) then begin
        print,'Rough cross-correlation shifts and rotation angle:'    
        print,'Xshift = ',strtrim(xshift,2)
        print,'Yshift = ',strtrim(yshift,2)
        print,'Angle = ',strtrim(angle,2),' degrees'
        print,'XFLIP'
        print,'Nsig = ',strtrim(nsig,2)
        print,'Rough number of Xcorr matches = ',strtrim(matchnum,2)
        print,''
      endif

    ; No good solution
    endif else begin
      if not keyword_set(silent) then $
        print,'No good LARGE Rotation solution'
    endelse
  

  endif ; flipped rotations



  ; Refining the solution, checking smaller angles intervals
  ;----------------------------------------------------------
  if (nsig gt 5.0) then begin

    if not keyword_set(silent) then $
      print,'Refining the solution'

    nrot = 11
    rotarr3 = findgen(nrot)*1.0-5.0 + angle
    nsigarr3 = fltarr(nrot)
    xshiftarr3 = fltarr(nrot)
    yshiftarr3 = fltarr(nrot)
    nmatcharr3 = fltarr(nrot)

    for i=0,nrot-1 do begin

      rxx2 =  (xx2-xhalf)*cos(rotarr3[i]/!radeg) + (yy2-yhalf)*sin(rotarr3[i]/!radeg) + xhalf
      ryy2 = -(xx2-xhalf)*sin(rotarr3[i]/!radeg) + (yy2-yhalf)*cos(rotarr3[i]/!radeg) + yhalf
      if xflip eq 1 then rxx2 = maxx2-rxx2   ; flip x-coordinates
      
      ;box = 0.5*(maxx < maxy)
      ;gd1 = where(xx1 ge (xhalf-box*0.5) and xx1 le (xhalf+box*0.5) and $
      ;            yy1 ge (yhalf-box*0.5) and yy1 le (yhalf+box*0.5),ngd1)
      ;gd2 = where(rxx2 ge (xhalf-box*0.5) and rxx2 le (xhalf+box*0.5) and $
      ;            ryy2 ge (yhalf-box*0.5) and ryy2 le (yhalf+box*0.5),ngd2)
      ;xx1s = xx1[gd1]-(xhalf-box*0.5)
      ;yy1s = yy1[gd1]-(yhalf-box*0.5)
      ;xx2s = rxx2[gd2]-(xhalf-box*0.5)
      ;yy2s = ryy2[gd2]-(yhalf-box*0.5)

      gd2 = where(rxx2 ge 0 and rxx2 le max(xx1) and ryy2 ge 0 and ryy2 le max(yy1),ngd2)
      xx1s = xx1
      yy1s = yy1
      xx2s = rxx2[gd2]
      yy2s = ryy2[gd2]      

      MATCHSTARS_XCORR,xx1s,yy1s,xx2s,yy2s,xshift,yshift,angle0,bestcorr,xcorr,silent=silent,$
                       smooth=5,xyscale=4,fwhm=5,matchnum=matchnum,nsig=nsig,/extra

      nsigarr3[i] = nsig
      xshiftarr3[i] = xshift
      yshiftarr3[i] = yshift 
      nmatcharr3[i] = matchnum
    end
       
    ; Find best rotation, spline
    rotarr3b = scale_vector(findgen(3600),min(rotarr3),max(rotarr3))
    nsigarr3b = spline(rotarr3,nsigarr3,rotarr3b)
    xshiftarr3b = spline(rotarr3,xshiftarr3,rotarr3b)
    yshiftarr3b = spline(rotarr3,yshiftarr3,rotarr3b)
    nmatcharr3b = spline(rotarr3,nmatcharr3,rotarr3b)
    bestind3 = first_el(maxloc(nsigarr3b))
    nsig = nsigarr3b[bestind3]
    angle = rotarr3b[bestind3]
    xshift = xshiftarr3b[bestind3]
    yshift = yshiftarr3b[bestind3]
    matchnum = nmatcharr3b[bestind3]

    if not keyword_set(silent) then begin
      print,'Rough cross-correlation shifts and rotation angle:'
      print,'Xshift = ',strtrim(xshift,2)
      print,'Yshift = ',strtrim(yshift,2)
      print,'Angle = ',strtrim(angle,2),' degrees'
      if xflip eq 1 then print,'XFLIP'
      print,'Nsig = ',strtrim(nsig,2)
      print,'Rough number of Xcorr matches = ',strtrim(matchnum,2)
      print,''
    endif

  endif  ; refining

end  ; still bad


; No good cross-correlation fit
if (nsig lt 5.0) then begin
  undefine,ind1,ind2,trans
  if not keyword_set(silent) then $
    print,'NO GOOD CROSS-CORRELATION'
  count = 0
  return
endif





;###################################################
; MATCH UP with initial rough XSHIFT/YSHIFT/ANGLE
;###################################################

; Transform coordinates(2) to coordinates(1)
; The xshift/yshift from Xcorr are for the center
; but for MPFIT they are for (0,0)

; The angle must be correct to within about +/-1 deg
; or we won't get enough good matches

; Take into account that we are rotating around the origin
; and not the center of the field
xoff = xhalf - xhalf*cos(angle/!radeg) - yhalf*sin(angle/!radeg)
yoff = yhalf - yhalf*cos(angle/!radeg) + xhalf*sin(angle/!radeg)
par = [xoff, yoff, cos(angle/!radeg), -sin(angle/!radeg), $   
       sin(angle/!radeg), cos(angle/!radeg)]
if xflip eq 1 then begin
  par[[0,2,4]] = -par[[0,2,4]]    ; XFLIP
  par[0] = par[0] + maxx2
endif
par[[0,1]] = par[[0,1]] + [xshift, yshift]

out = trans_coo(xx2,yy2,par)
xx2c = reform(out[0,*])
yy2c = reform(out[1,*])

; Match them
dcr = 15.  ;10.0
SRCMATCH,xx1,yy1,xx2c,yy2c,dcr,ind1a,ind2a,count=nind1a
if nind1a eq 0 then begin
  undefine,ind1,ind2,trans
  if not keyword_set(silent) then $
    print,'NO MATCHES'
  count = 0
  return
endif
xx1m = xx1[ind1a]
yy1m = yy1[ind1a]
xx2m = xx2[ind2a]
yy2m = yy2[ind2a]

dum = where(ind1a ne -1,nind1a)
if not keyword_set(silent) then begin
  print,strtrim(nind1a,2),' matches from Xcorr solution'
  print,''
endif
diff = sqrt((xx1m-xx2c[ind2a])^2.0 + (yy1m-yy2c[ind2a])^2.0)
rms = sqrt(mean(diff^2.0))

; Plot the residuals
if keyword_set(plotresid) then begin
  plot,yy1m,xx1m-xx2c[ind2a],ps=3
  oplot,xx1m,yy1m-yy2c[ind2a],ps=3,co=250
  stop
endif




;#################################################
; ITERATIVE LINE FITTING
;#################################################


MATCHSTARS_LINEFIT,par,xx1m,yy1m,xx2m,yy2m,fpar1,rms1,plotresid=plotresid,silent=silent

if n_elements(fpar1) eq 0 or rms1 gt 100 then begin
  undefine,ind1,ind2,trans
  rms = 999999.
  count=0
  return
endif
; Match them using the newer transformation
out = trans_coo(xx2,yy2,fpar1)
xx2c = reform(out[0,*])
yy2c = reform(out[1,*])
dcr = (3.0*rms1) < 10.
SRCMATCH,xx1,yy1,xx2c,yy2c,dcr,ind1b,ind2b,count=nind1b
if nind1b eq 0 then begin
  undefine,ind1,ind2,trans
  if not keyword_set(silent) then $
    print,'NO MATCHES'
  count = 0
  return
endif
xx1m = xx1[ind1b]
yy1m = yy1[ind1b]
xx2m = xx2[ind2b]
yy2m = yy2[ind2b] 
nind1b = n_elements(ind1b)

if not keyword_set(silent) then begin
  print,'Linefit finished'
  print,'RMS = ',strtrim(rms1,2)
  print,'Nmatch = ',strtrim(nind1b,2)
  print,''
endif


; Need at least 6 stars to use MPFIT
; since there are 6 parameters
if (nind1b lt 6) then begin
  if not keyword_set(silent) then $
    print,'Not enough matches for MPFIT'

  ; Return what we got so far
  if rms1 lt 100 and nind1b gt 0 then begin
    ; Refit trans with original coordinates (w/o minx/miny offsets)
    MATCHSTARS_LINEFIT,fpar1,xx1m+minx,yy1m+miny,xx2m+minx,yy2m+miny,fpar,rms
    ind1 = ind1b
    ind2 = ind2b
    trans = fpar
    count = nind1b
    rms = rms
  ; Noting useful to return
  endif else begin
    undefine,ind1,ind2,trans
    count=0
    rms = 999999.
  endelse
  return
endif



;########################################
; MPFIT - FIRST time
;########################################

; Now do a robust fit for the transformation
func='trans_coo_dev'
fa = {x1:double(xx1m), y1:double(yy1m), x2:double(xx2m), y2:double(yy2m)}
par = double(fpar1)
fpar2 = MPFIT(func,par,functargs=fa, perror=perror, niter=iter, status=status,$
             bestnorm=chisq, dof=dof, autoderivative=1, /quiet)  ;ftol=1d-10

if (status lt 1) then begin
  if not keyword_set(silent) then $
    print,'MPFIT problem'
  undefine,ind1,ind2,trans
  count=0
  return
endif

diff2 = trans_coo_dev(fpar2,x1=xx1m,y1=yy1m,x2=xx2m,y2=yy2m)
rms2 = sqrt(mean(diff2^2.))


;--------------------------------
; Match-up original coordinates
;--------------------------------
; Need to apply the MINX/MINY offset that we did at the beginning
out = trans_coo(xinp2-minx,yinp2-miny,fpar2)
txinp2 = reform(out[0,*])
tyinp2 = reform(out[1,*])

; Match again
dcr = ((rms2*3.0) < 3.0) > 0.0001
; Need to apply the MINX/MINY offset that we did at the beginning
SRCMATCH,xinp1-minx,yinp1-miny,txinp2,tyinp2,dcr,ind1c,ind2c,count=nind1c
if nind1c eq 0 then begin
  undefine,ind1,ind2,trans
  if not keyword_set(silent) then $
    print,'NO MATCHES'
  count = 0
  return
endif
xinp1m = xinp1[ind1c]
yinp1m = yinp1[ind1c]
xinp2m = xinp2[ind2c]
yinp2m = yinp2[ind2c]


; Plot residuals
if keyword_set(plotresid) then begin
  plot,yinp1m,(xinp1m-minx)-txinp2[ind2c],ps=3
  oplot,xinp1m,(yinp1m-miny)-tyinp2[ind2c],ps=3,co=250
  print,fpar2
  print,rms2
  stop
endif


; Second fit with MPFIT
func='trans_coo_dev'
fa = {x1:xinp1m, y1:yinp1m, x2:xinp2m, y2:yinp2m}
;ftol = 1e-10
par = fpar2
fpar3 = MPFIT(func,par,functargs=fa, perror=perror, niter=iter, status=status,$
             bestnorm=chisq, dof=dof, autoderivative=1, /quiet)

if (status lt 1) then begin
  if not keyword_set(silent) then $
    print,'MPFIT problem'
  undefine,ind1,ind2,trans
  count=0
  return
endif

diff3 = trans_coo_dev(fpar3,x1=xinp1m,y1=yinp1m,x2=xinp2m,y2=yinp2m)
rms3 = sqrt(mean(diff3^2.))

sigpar = perror * sqrt(chisq/dof)
chisq2 = chisq/n_elements(xinp1m)



;##############################################
; Final matchup with the ORIGINAL coordinates
;##############################################
out = trans_coo(xinp2,yinp2,fpar3)
fxinp2 = reform(out[0,*])
fyinp2 = reform(out[1,*])

; Match again
dcr = ((3.0*rms3) < 3.0) > 0.0001
SRCMATCH,xinp1,yinp1,fxinp2,fyinp2,dcr,ind1,ind2,count=nind1
if nind1 eq 0 then begin
  undefine,ind1,ind2,trans
  if not keyword_set(silent) then $
    print,'NO MATCHES'
  count = 0
  return
endif
xinp1m = xinp1[ind1]
yinp1m = yinp1[ind1]
xinp2m = xinp2[ind2]
yinp2m = yinp2[ind2]

diff4 = trans_coo_dev(fpar3,x1=xinp1m,y1=yinp1m,x2=xinp2m,y2=yinp2m)
rms4 = sqrt(mean(diff4^2.))

; Plot the residuals
if keyword_set(plotresid) then begin
  plot,yinp1m,xinp1m-fxinp2[ind2],ps=3
  oplot,xinp1m,yinp1m-fyinp2[ind2],ps=3,co=250
  print,fpar3
  print,rms3
  stop
endif

nmatch = n_elements(ind1)
count = nmatch     ; the final number of matches

; Final results
if not keyword_set(silent) then begin
  print,'Final MATCHSTARS Results:'
  ;print,strtrim(nmatch,2),' MATCHES within '+string(dcr,format='(F5.2)')+' pixels'
  print,strtrim(nmatch,2),' MATCHES within '+strtrim(dcr,2)+' pixels'
  print,'RMS = ',strtrim(rms4,2),' pixels'
endif



; Final transformation equation
; trans = [A, B, C, D, E, F]
; x(1) = A + C*x(2) + E*y(2)
; y(1) = B + D*x(2) + F*y(2)
trans = fpar3

; Final RMS
rms = rms4

; Runtime
if not keyword_set(silent) then $
  print,'Runtime = ',stringize(systime(1)-t0,ndec=2),' seconds'


if keyword_set(stp) then stop

end
