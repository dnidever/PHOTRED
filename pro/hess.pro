;+
;  NAME:
;   HESS
;
;  PURPOSE:
;   This program makes a "Hess"-like density diagram of input arrays.
;
;  INPUTS:
;   x        Array of x-values
;   y        Array of y-values
;   z        Array of z-values.  Required for /total,/avg,/mean,/std,/maximum, perc=
;   cut      A cut made in the data that can be made with the IDL
;            WHERE statement.  This should be a string such as 'X gt 5.0'.
;               Use 'X','Y','Z' to make cuts on the x/y/z input arrays
;   /color   Plot in color
;   dx       Interval in x-axis.  Either set the # of bins or interval, not both
;   dy       Interval in y-axis
;   nx       Number of bins in x-axis.  Either set the # of bins or interval, not both
;   ny       Number of bins in y-axis
;   xtit     Title for x-axis on plot
;   ytit     Title for y-axis on plot
;   tit      Overall plot title
;   /log     Plot the logarithm.
;   /interp  Interpolate between the points
;   /noplot  Don't plot anything
;   /save    Save the plot to a poscript file
;   file     Name of postscript file to save the plot to
;   xrange   X-axis plot range. If center=1 (default) then this should
;            be the range of the centers of the pixels.  The order of
;            the elements can be reversed to flip the axis (identical
;            to setting /xflip)
;   yrange   Y-axis plot range. If center=1 (default) then this should
;            be the range of the centers of the pixels.  The order of
;            the elements can be reversed to flip the axis (identical
;            to setting /yflip)
;   /stp     Stop at the end of the program
;   top      The largest value to plot
;   bot      The lowest value to plot
;   /total   Show total of z values
;   /avg     Color indicates the average of the third dimension
;   /std     The standard deviation of zarr
;   /mean    The mean of zarr
;   /maximum The maximum of zarr
;   =perc    Use the Nth percentile (perc=90, 90th percentile).
;   /colnorm Normalize each column by the number of points in the column
;   /rownorm Normalize each row by the number of points in each row
;   /nocolorbar   Don't overplot the colorbar
;   position The position to put the plot in, in normalized coordinates
;   /onlyim  Don't plot the colorbar and don't leave space for it.
;   /noerase Don't erase the screen before you plot the image
;   charsize Character size
;   minhue   The minimum hue value (between 0 and 360). For /AVG
;   maxhue   The maximum hue value (between 0 and 360). For /AVG
;   minbright The minimum brightness value (between 0 and 1). For /AVG
;   maxbright The maximum brightness value (between 0 and 1). For /AVG
;   saturation The staturation value (between 0 and 1). For /AVG
;   posim    Position of image (in normalized coordinates)
;   poscol   Positin of color bar (in normalized coordinates)
;   thick    The thickness of lines
;   charthick The thickness of the annotations
;   framecolor  The color of the box and annotations
;   background_color The color of the background
;   /center   The coordinates plotted should be the center of bins (default).
;   /force    Force the xrange and yrange to be exactly as the way they are set (default)
;   =format   The format to use for the colorbar annotation.
;   =xsize    Sets the x image size for PS output.
;   =ysize    Sets the y image size for PS output.
;   =zbot     Set the bottom of the z-scale for /AVG
;   =ztop     Set the top of the z-scale for /AVG
;   =maskvalue    The value to which bad/missing data is flagged.  Pixels
;                   with this value are set to 0B
;   =maskcolor    The color to use for "bad data".  The default is 0B.
;
;  OUTPUTS
;   im        The image that is plotted
;   =xout     The array of x-values
;   =yout     The array of y-values
;   =zout     The array of z-values
;   =xarr     The x-values of the image
;   =yarr     The y-values of the image
;   =indim    The 2D array giving the index of points used for IM
;               for =PERC and /MAXIMUM.
;
;  PROGRAMS USED:
;   DISPLAY.PRO
;   IMGSCL.PRO
;   COLORBAR.PRO
;   SCALE_VECTOR.PRO
;   CLOSEST.PRO
;   MINLOC.PRO
;   FIRST_EL.PRO
;
; Created by David Nidever August 2005
;-


pro hess_hist,im,x,y,xmin,ymin,dx,dy,colnorm=colnorm,rownorm=rownorm

; This little function does the histogram part of HESS.PRO

ng = n_elements(x)

; Number of points
;for i=0.,ng-1. do dum = ++im( floor((x(i)-xmin)/dx), floor((y(i)-ymin)/dy) )
++im( floor((x-xmin)/dx), floor((y-ymin)/dy) )

; COLNORM
if keyword_set(colnorm) then begin
  sz = size(im)
  imn = im  ; number of points
  for i=0LL,sz[1]-1LL do begin
    tot = (total(imn(i,*)) > 1.)
    im(i,*) /= tot
  endfor
endif ; /colnorm

; ROWNORM
if keyword_set(rownorm) then begin
  sz = size(im)
  imn = im  ; number of points
  for i=0LL,sz(2)-1LL do begin
    tot = (total(imn(*,i)) > 1.)
    im(*,i) /= tot
  endfor
endif ; /rownorm

end


;-------------------------------------------------------------------------

pro hess_total,im,x,y,z,xmin,ymin,dx,dy

ng = n_elements(x)

; Total
for i=0LL,ng-1LL do begin
  xind = floor((x(i)-xmin)/dx)
  yind = floor((y(i)-ymin)/dy)
  im(xind,yind) += z(i)
endfor

end

;-------------------------------------------------------------------------

pro hess_std,im,x,y,z,xmin,ymin,dx,dy,colnorm=colnorm,rownorm=rownorm

; STD (standard deviation)

ng = n_elements(x)

; Number of points
imn = im      ; number of gaussians
;for i=0.,ng-1. do dum = ++imn( floor((x(i)-xmin)/dx), floor((y(i)-ymin)/dy) )
++imn( floor((x-xmin)/dx), floor((y-ymin)/dy) )

; Finding the mean
ima = im
for i=0LL,ng-1LL do begin
  xind = floor((x(i)-xmin)/dx)
  yind = floor((y(i)-ymin)/dy)
  ima(xind,yind) = ima(xind,yind)+z(i)
endfor
ima = ima/(imn > 1.)   ; divide by number of points

; Std. Dev.
; stdev = SQRT( SUM( (x-meanx)^2. )/(N-1) )
ims = im      ; standard deviation 
for i=0LL,ng-1LL do begin
  xind = floor((x(i)-xmin)/dx)
  yind = floor((y(i)-ymin)/dy)
  ims(xind,yind) += (z(i)-ima(xind,yind))^2.
endfor
ims = ims/( (imn-1) > 1. )   ; divide by # of points
ims = sqrt(ims)              ; square-root

; COLNORM
if keyword_set(colnorm) then begin
  sz = size(im)
  for i=0LL,sz[1]-1LL do begin
    tot = (total(imn(i,*)) > 1.)
    ims(i,*) /= tot
  endfor
endif ; /colnorm

; ROWNORM
if keyword_set(rownorm) then begin
  sz = size(im)
  for i=0LL,sz(2)-1LL do begin
    tot = (total(imn(*,i)) > 1.)
    ims(*,i) /= tot
  endfor
endif ; /rownorm

; Creating the final image
im = ims

end

;-------------------------------------------------------------------------

pro hess_max,im,x,y,z,xmin,ymin,dx,dy,indim=indim

; MAXIMUM

ng = n_elements(x)

indim = im*0L-1

; Finding the maxmimum
;ima = im
for i=0LL,ng-1LL do begin
  xind = floor((x(i)-xmin)/dx)
  yind = floor((y(i)-ymin)/dy)
  if z[i] gt im[xind,yind] then indim[xind,yind]=i
  im(xind,yind) = im(xind,yind) > z(i)
endfor

end

;-------------------------------------------------------------------------

pro hess_mean,im,x,y,z,xmin,ymin,dx,dy,colnorm=colnorm,rownorm=rownorm

; MEAN

ng = n_elements(x)

; Number of points
imn = im      ; number of gaussians
;for i=0.,ng-1. do dum = ++imn( floor((x(i)-xmin)/dx), floor((y(i)-ymin)/dy) )
++imn( floor((x-xmin)/dx), floor((y-ymin)/dy) )

; Finding the mean
ima = im
for i=0LL,ng-1LL do begin
  ;xind = floor((x(i)-xmin)/dx)
  ;yind = floor((y(i)-ymin)/dy)
  ;ima(xind,yind) = ima(xind,yind)+z(i)
  ima[ floor((x[i]-xmin)/dx), floor((y[i]-ymin)/dy) ] += z[i]
endfor
ima = ima/(imn > 1.)   ; divide by number of points

; COLNORM
if keyword_set(colnorm) then begin
  sz = size(im)
  for i=0LL,sz[1]-1LL do begin
    tot = (total(imn(i,*)) > 1.)
    ima(i,*) /= tot
  endfor
endif ; /colnorm

; ROWNORM
if keyword_set(rownorm) then begin
  sz = size(im)
  for i=0LL,sz(2)-1LL do begin
    tot = (total(imn(*,i)) > 1.)
    ima(*,i) /= tot
  endfor
endif ; /rownorm

; Creating the final image
im = ima

end

;-------------------------------------------------------------------------


pro hess_avg,im,x,y,z,xmin,ymin,dx,dy,nx,ny,colnorm=colnorm,rownorm=rownorm,bot=bot,top=top,$
             zbot=zbot,ztop=ztop,log=log,minhue=minhue,maxhue=maxhue,minbright=minbright,$
             maxbright=maxbright,saturation=saturation,maskcolor=maskcolor,maskvalue=maskvalue

; AVERAGE

ng = n_elements(x)

if n_elements(z) eq 0 then begin
  print,'ZARR IS REQUIRED WITH /AVG'
  return
endif

imn = im      ; number of gaussians
++imn( floor((x-xmin)/dx), floor((y-ymin)/dy) )

; Total and Average
ima = im      ; average
for i=0LL,ng-1LL do begin
  xind = floor((x(i)-xmin)/dx)
  yind = floor((y(i)-ymin)/dy)
  ima(xind,yind) += z(i)
  ;dum = ++imn( floor((x(i)-xmin)/dx), floor((y(i)-ymin)/dy) )
endfor

im = imn              ; # of points
ima = ima/(imn > 1.)  ; average

; COLNORM
if keyword_set(colnorm) then begin
  sz = size(im)
  for i=0LL,sz[1]-1LL do begin
    tot = (total(imn(i,*)) > 1.)
    ima(i,*) /= tot
  endfor
endif ; /colnorm

; ROWNORM
if keyword_set(rownorm) then begin
  sz = size(im)
  for i=0LL,sz(2)-1LL do begin
    tot = (total(imn(*,i)) > 1.)
    ima(*,i) /= tot
  endfor
endif ; /rownorm

; Setting TOP/BOT if not input
;if (n_elements(top) eq 0) then top=max(ima)
;if (n_elements(bot) eq 0) then bot=min(ima)
if (n_elements(top) eq 0) then top=max(im)
if (n_elements(bot) eq 0) then bot=min(im)
;if (n_elements(ztop) eq 0) then ztop=max(ima)
;if (n_elements(zbot) eq 0) then zbot=min(ima)


; TOP
if keyword_set(top) then begin
  ;tbd = where(ima gt top,ntbd)
  ;if ntbd gt 0 then (ima)(tbd) = top
  tbd = where(im gt top,ntbd)
  if ntbd gt 0 then (im)(tbd) = top
endif

; BOT
if keyword_set(bot) then begin
  ;bbd = where(ima lt bot,nbbd)
  ;if nbbd gt 0 then (ima)(bbd) = bot
  bbd = where(im lt bot,nbbd)
  if nbbd gt 0 then (im)(bbd) = bot
endif

byte_im = ImgScl(im, Min=bot, Max=top, $
		 Log=log, Levels=l, MaskValue=maskvalue)
im = byte_im

; Keeping the Average within the Z range
low = where(ima lt min(z),nlow)
if nlow gt 0 then ima[low] = min(z)
hi = where(ima gt max(z),nhi)
if nhi gt 0 then ima[hi] = max(z)

; Setting ZTOP/BOT if not input
if (n_elements(ztop) eq 0) then ztop=max(ima)
if (n_elements(zbot) eq 0) then zbot=min(ima)

; ZTOP, z-values
if keyword_set(ztop) then begin
  ztbd = where(ima gt ztop,nztbd)
  if nztbd gt 0 then (ima)(ztbd) = ztop
endif

; ZBOT, z-values
if keyword_set(zbot) then begin
  zbbd = where(ima lt zbot,nzbbd)
  if nzbbd gt 0 then (ima)(zbbd) = zbot
endif


; convert image to RGB, using HLS
; hue is Average (IMA)   (0-360)
; 0-red, 120-green, 240-blue
; brightness is Total (im) (0-1)
ima2 = -ima    ; (blue-green-red)
if n_elements(minhue) eq 0 then minhue = 0.
if n_elements(maxhue) eq 0 then maxhue = 240.
if n_elements(minbright) eq 0 then minbright = 0.10
if n_elements(maxbright) eq 0 then maxbright = 0.70
if n_elements(saturation) eq 0 then saturation = 0.9

;if n_elements(zbot) eq 0 then zbot=min(ima)
;if n_elements(ztop) eq 0 then ztop=max(ima)



hue = scale(ima2,[-ztop,-zbot],[minhue,maxhue])
;hue = scale_vector(ima2,minhue,maxhue)
bright = scale_vector(im,minbright,maxbright)

color_convert, hue, bright, im*0.+saturation, r, g, b, /HLS_RGB
;color_convert, hue, im*0.+1.0, bright, r, g, b, /HSV_RGB

; BUG
; setting bottom to zero
;if keyword_set(bot) then begin
;  bbd = where(ima lt bot,nbbd)
;  if nbbd gt 0 then begin
;    (r)(bbd) = 0.
;    (g)(bbd) = 0.
;    (b)(bbd) = 0.
;  endif
;endif

;bd = where(ima eq 0. or im eq 0.,nbd)
;if nbd gt 0 then begin
;  (r)(bd) = 255
;  (g)(bd) = 255
;  (b)(bd) = 255
;endif

; Interleaved image
im2 = bytarr(3,nx,ny)
im2(0,*,*) = r
im2(1,*,*) = g
im2(2,*,*) = b

; Final image to plot
orig_im = im
im = im2

;if not keyword_set(tit) then tit = 'COLOR AS '+strupcase(zname)
;stop

end


;-------------------------------------------------------------------------


pro hess_perc,im,x,y,z,xmin,ymin,dx,dy,nx,ny,perc=perc,indim=indim

; Using the Nth percentile

ng = n_elements(x)

fim = im*0.
indim = im*0L-1

imn = im      ; number of points
++imn( floor((x-xmin)/dx), floor((y-ymin)/dy) )

; Number of points
indarr = x*0.
for i=0LL,ng-1LL do begin
  ;dum = ++imn( floor((x(i)-xmin)/dx), floor((y(i)-ymin)/dy) )

  ; 1D array index is ind = ix + iy*nx
  indarr(i) = floor((x(i)-xmin)/dx) + floor((y(i)-ymin)/dy) * nx
endfor

; 1D array of # of points in each bin
binpop1d = (imn)(*) 
nbins = n_elements(binpop1d)

; We want the points to be sorted in the bin order (all points in bin 0
; first, then all the points in bin 1; but within each bin the
; points are not sorted).
; First we create the index of where each bin will start in the large array
; It's basically the sum of all points in the bins before the current
binind = binpop1d*0.
; Add up the number of points in all bins before this one.
for i=1LL,nbins-1LL do binind(i) =  binind(i-1)+binpop1d(i-1)

; Now go through and put the data points in a new array ordered by bin
zord = z*0.                 ; the z array ordered by bin
origindarr = indarr*0.
curbinpop1d = binpop1d*0.   ; the number of points per bin so far

; Looping through the data points
for i=0LL,ng-1LL do begin
  ; indarr will tell us what bin this data point is in
  ; curbinpop1d will tell us which number it is for its bin

  ; This points bin number in 1D array
  bin = indarr(i)

  ; Index for this point.  Where the bin starts + how many points
  ; have already been processed
  ind = binind(bin) + curbinpop1d(bin)

  ; Putting in the point and its original index
  zord(ind) = z(i)
  origindarr(ind) = i

  ; Incrementing the number of points processed for this bin
  ++curbinpop1d(bin)
endfor

; Loop through each bin and sort its points
; Then pick out the Nth percentile one
for i=0L,nbins-1LL do begin
    
  ; Index of this bin in the 2D array
  ind2 = array_indices(im,i)

  ; We've got points in this bin
  if (binpop1d(i) gt 0) then begin

    ; Getting the data points for this bin in the sorted z array
    points = zord(binind(i):binind(i)+binpop1d(i)-1)
    points_origind = origindarr(binind(i):binind(i)+binpop1d(i)-1)
    npoints = n_elements(points)
    si = sort(points)
    points = points(si)
    points_origind = points_origind[si]

    ; Getting the index of the Nth percentile point
    gd = ceil(npoints*perc/100.)-1
    gd = (0 > gd) < (npoints-1)         ; making sure it's in the range

    ;stop    

    fim[ind2[0],ind2[1]] = points[gd]
    indim[ind2[0],ind2[1]] = points_origind[gd]

  ; No points in this bin
  endif else begin
   ; fim(ind2[0],ind2[1]) = 0.
  endelse

endfor

im = fim

end


;-------------------------------------------------------------------------


pro hess,xx,yy,zz,im,color=color,nx=nx0,ny=ny0,dx=dx0,dy=dy0,$
          xtit=xtit,ytit=ytit,tit=tit,log=log,interpolate=interpolate,$
          noplot=noplot,save=save,file=file,xrange=xrange0,yrange=yrange0,$
          stp=stp,cut=cut,top=top,total=total,bot=bot,$
          nocolorbar=nocolorbar,position=position,$
          onlyim=onlyim,noerase=noerase,charsize=charsize,$
          avg=avg,minhue=minhue,maxhue=maxhue,minbright=minbright,$
          maxbright=maxbright, saturation=saturation, xflip=xflip,$
          yflip=yflip,year=year,ctool=ctool,xout=xout,yout=yout,zout=zout,$
          posim=posim,poscol=poscol,thick=thick,charthick=charthick,$
          framecolor=framecolor,background_color=background_color,$
          std=std,mean=mean,colnorm=colnorm,$
          rownorm=rownorm,xarr=xarr,yarr=yarr,center=center,zname=zname,$
          maximum=maximum,force=force,perc=perc,format=format,$
          xsize=xsize,ysize=ysize,zbot=zbot,ztop=ztop,$
          xtickformat=xtickformat,ytickformat=ytickformat,xstyle=xstyle,$
          ystyle=ystyle,indim=indim,maskvalue=maskvalue,maskcolor=maskcolor



nxarr = n_elements(xx)
nyarr = n_elements(yy)

; Wrong inputs
if (n_params() eq 0) or (nxarr eq 0) or (nyarr eq 0) then begin
  print,'Syntax - hess,x,y,z,im,color=color,nx=nx,ny=ny,dx=dx,dy=dy,'
  print,'             xtit=xtit,ytit=ytit,tit=tit,log=log,interpolate=interpolate,'
  print,'             noplot=noplot,save=save,file=file,xrange=xrange,yrange=yrange,'
  print,'             stp=stp,cut=cut,top=top,total=total,bot=bot,'
  print,'             nocolorbar=nocolorbar,position=position,'
  print,'             onlyim=onlyim,noerase=noerase,charsize=charsize,'
  print,'             avg=avg,minhue=minhue,maxhue=maxhue,minbright=minbright,'
  print,'             maxbright=maxbright, saturation=saturation, xflip=xflip,'
  print,'             yflip=yflip,year=year,ctool=ctool,xout=xout,yout=yout,zout=zout,'
  print,'             posim=posim,poscol=poscol,thick=thick,charthick=charthick,'
  print,'             framecolor=framecolor,background_color=background_color,'
  print,'             std=std,mean=mean,colnorm=colnorm,indim=indim,'
  print,'             rownorm=rownorm,maximum=maximum,perc=perc,zbot=zbot,ztop=ztop'
  return
endif


; Making temporary arrays
x = xx
y = yy
if n_elements(zz) gt 0 then z = zz
if n_elements(dx0) gt 0 then dx = dx0
if n_elements(dy0) gt 0 then dy = dy0
if n_elements(nx0) gt 0 then nx = nx0
if n_elements(ny0) gt 0 then ny = ny0
if n_elements(xrange0) gt 0 then xrange = xrange0
if n_elements(yrange0) gt 0 then yrange = yrange0

; Checking NANs
if n_elements(zz) gt 0 then begin
  g = where(finite(x) eq 1 and finite(y) eq 1 and finite(z) eq 1,ng)
  if ng eq 0 then begin
    print,'NO GOOD POINTS'
    return
  endif
  x = x[g]
  y = y[g]
  z = z[g]
endif else begin
  g = where(finite(x) eq 1 and finite(y) eq 1,ng)
  if ng eq 0 then begin
    print,'NO GOOD POINTS'
    return
  endif
  x = x[g]
  y = y[g]
endelse


; SETTING PARAMETERS

; Histogram, by default
if not keyword_set(mean) and not keyword_set(avg) and not keyword_set(std) $
   and not keyword_set(total) and not keyword_set(maximum) and not keyword_set(perc) then hist=1

if n_elements(force) eq 0 then force=1      ; force xrange and yrange to be exact

; Z not set
if (keyword_set(avg) or keyword_set(mean) or keyword_set(std) or keyword_set(maximum) $
    or keyword_set(perc)) and n_elements(z) eq 0 then begin
  print,'ZARR IS REQUIRED WITH /AVG'
  return
endif


; Making the cut
if keyword_set(cut) then begin
  outcut = execute('gd=where('+cut+',ngd)')
  x = x(gd)
  y = y(gd)
  if n_elements(z) gt 0 then z = z(gd)
endif

; Setting NX/NY=200 by default if not input
if (n_elements(nx) eq 0) then nx = 200.
if (n_elements(ny) eq 0) then ny = 200.

; Min and Max's
ymin = min(y)
ymax = max(y)
xmin = min(x)
xmax = max(x)

; Temporary DX and DY if not input
; X-axis
if (n_elements(nx) gt 0) and (n_elements(dx) eq 0) then begin
  ; XRANGE set
  if (n_elements(xrange) ne 0) then begin
    dx = (xrange[1]-xrange[0])/(nx-1.)
  ; XRANGE not set
  endif else begin
    dx = (xmax-xmin)/(nx-1.)
  endelse
endif
; Y-axis
if (n_elements(ny) gt 0) and (n_elements(dy) eq 0) then begin
  ; YRANGE set
  if (n_elements(yrange) ne 0) then begin
    dy = (yrange[1]-yrange[0])/(ny-1.)
  ; YRANGE not set
  endif else begin
    dy = (ymax-ymin)/(ny-1.)
  endelse
endif


; Axes reversed
if n_elements(xrange) gt 0 then begin
  ; X axis reversed
  if xrange[1] lt xrange[0] then begin
    xrange = reverse(xrange)
    xflip = 1
  endif
endif
if n_elements(yrange) gt 0 then begin
  ; Y axis reversed
  if yrange[1] lt yrange[0] then begin
    yrange = reverse(yrange)
    yflip = 1
  endif
endif

; Cutting based on the xrange/yrange if set
;   Otherwise the image size might get too large if there are
;   points far outside the xrange/yrange
; Xrange set
if (n_elements(xrange) gt 0) then begin
  cut2 = '(x ge (xrange[0]-dx)) and (x le (xrange[1]+dx))'
  outcut2 = execute('gd2=where('+cut2+',ngd2)')

  orig_x = x
  orig_y = y

  if ngd2 gt 0 then begin
    x = x(gd2)
    y = y(gd2)
    if n_elements(z) gt 0 then z = z(gd2)
  endif else begin
    x = -999999.
    y = -999999.
    if n_elements(z) gt 0 then z=-999999.
    print,'NO POINTS LEFT'
    return
  endelse  

  ; Min and Max's
  ymin = min(y)
  ymax = max(y)
  xmin = min(x)
  xmax = max(x)

  ; Setting DX/DY if not input
  if (n_elements(dx0) eq 0) then dx = (xmax-xmin)/(nx-1.)
  if (n_elements(dy0) eq 0) then dy = (ymax-ymin)/(ny-1.)

endif
; Yrange set
if (n_elements(yrange) gt 0) then begin
  cut3 = '(y ge yrange(0)-dy) and (y le yrange[1]+dy)'
  outcut3 = execute('gd3=where('+cut3+',ngd3)')

  if ngd3 gt 0 then begin
    x = x(gd3)
    y = y(gd3)
    if n_elements(z) gt 0 then z = z(gd3)
  endif else begin
    x = -999999.
    y = -999999.
    if n_elements(z) gt 0 then z=-999999.
    print,'NO POINTS LEFT'
    return
  endelse  


  ; Min and Max's
  ymin = min(y)
  ymax = max(y)
  xmin = min(x)
  xmax = max(x)

  ; Setting DX/DY if not input
  if (n_elements(dx0) eq 0) then dx = (xmax-xmin)/(nx-1.)
  if (n_elements(dy0) eq 0) then dy = (ymax-ymin)/(ny-1.)
endif

; Step sizes must be positive
dx = abs(dx)
dy = abs(dy)

; Center the coordinates, normally the coordinates refer to the bottom-left corner
if (n_elements(center) eq 0) then center=1

; Setting the range
if (n_elements(xrange) gt 0) then begin

  ; Default is that the xrange/yrange are for the centers
  ; Move xmin/xmax accordingly
  if keyword_set(center) then off=0.5*dx else off=0.0

  ; xmin and xmax must be xrange+/-integer*dx
  if keyword_set(force) then begin
    ; If xmin < xrange[0] then move it down an integer number of dx's
    if (xmin lt (xrange[0]-off)) then begin
      diff = (xrange[0]-off)-xmin
      xmin = (xrange[0]-off)-ceil(diff/dx)*dx
    endif else xmin = xrange[0]-off
    ; If xmax > xrange[1] then move it up an integer number of dx's
    if (xmax gt (xrange[1]+off)) then begin
      diff = xmax-(xrange[1]+off)
      xmax = (xrange[1]+off)+ceil(diff/dx)*dx
    endif else xmax = xrange[1]+off

  ; Don't force the xrange
  endif else begin
    if ((xrange[0]-off) lt xmin) or ((xrange[1]+off) gt xmax) then begin
      xmin = (xrange[0]-off) < xmin
      xmax = (xrange[1]+off) > xmax
    endif
  endelse
endif
if (n_elements(yrange) gt 0) then begin

  ; Default is that the xrange/yrange are for the centers
  ; Move ymin/ymax accordingly
  if keyword_set(center) then off=0.5*dy else off=0.0

  ; ymin and ymax must be yrange+/-integer*dy
  if keyword_set(force) then begin
    ; If ymin < yrange[0] then move it down and integer number of dy's
    if (ymin lt (yrange[0]-off)) then begin
      diff = (yrange[0]-off)-ymin
      ymin = (yrange[0]-off)-ceil(diff/dy)*dy
    endif else ymin = yrange[0]-off
    if (ymax gt (yrange[1]+off)) then begin
      diff = ymax-(yrange[1]+off)
      ymax = (yrange[1]+off)+ceil(diff/dy)*dy
    endif else ymax = yrange[1]+off

  ; Don't force the yrange
  endif else begin
    if ((yrange[0]-off) lt ymin) or ((yrange[1]+off) gt ymax) then begin
      ymin = (yrange[0]-off) < ymin
      ymax = (yrange[1]+off) > ymax
    endif
  endelse
endif


; Setting final DX/DY and NX/NY
; X-axis
if (n_elements(dx0) eq 0) then begin
  dx = (xmax-xmin)/(nx-1.)
endif else begin
  nx = floor((xmax-xmin)/dx)+ 1.  ; only want bins fully within the range
endelse
; Y-axis
if (n_elements(dy0) eq 0) then begin
  dy = (ymax-ymin)/(ny-1.)
endif else begin
  ny = floor((ymax-ymin)/dy)+ 1.  ; only want bins fully within the range
endelse


; Starting image array
im = fltarr(nx,ny)
ng = n_elements(x)


; CREATING THE IMAGE

; Array of x and y values for image
yarr = dindgen(ny)*double(dy)+double(ymin)
xarr = dindgen(nx)*double(dx)+double(xmin)

; Center the coordinates, normally the coordinates refer to the bottom-left corner
;if n_elements(center) eq 0 then center=1
if keyword_set(center) then begin
 xarr = xarr+0.5*dx
 yarr = yarr+0.5*dy
end

; *** DOING THE CALCULATIONS ***
; These sub-programs are at the top of this file

; HISTOGRAM, Number of Points
if keyword_set(hist) then $
  hess_hist,im,x,y,xmin,ymin,dx,dy,colnorm=colnorm,rownorm=rownorm

; TOTAL
if keyword_set(total) then $
  hess_total,im,x,y,z,xmin,ymin,dx,dy

; STANDARD DEVIATION
if keyword_set(std) then $
  hess_std,im,x,y,z,xmin,ymin,dx,dy,colnorm=colnorm,rownorm=rownorm

; MAXIMUM
if keyword_set(maximum) then $
  hess_max,im,x,y,z,xmin,ymin,dx,dy,indim=indim

; MEAN
if keyword_set(mean) then $
  hess_mean,im,x,y,z,xmin,ymin,dx,dy,colnorm=colnorm,rownorm=rownorm

; AVERAGE
if keyword_set(avg) then $
  hess_avg,im,x,y,z,xmin,ymin,dx,dy,nx,ny,colnorm=colnorm,rownorm=rownorm,bot=bot,top=top,$
           zbot=zbot,ztop=ztop,log=log,minhue=minhue,maxhue=maxhue,minbright=minbright,$
           maxbright=maxbright,saturation=saturation

; Nth PERCENTILE
if keyword_set(perc) then $
  hess_perc,im,x,y,z,xmin,ymin,dx,dy,nx,ny,perc=perc,indim=indim



; TOP and BOT (not /avg, has its own)
if not keyword_set(avg) then begin

  ; top
  if (n_elements(top) gt 0) then begin
    bd = where(im ge top,nbd)
    if nbd gt 0 then im(bd) = top
  endif

  ; bot
  if (n_elements(bot) gt 0) then begin
    bd = where(im le bot,nbd)
    if nbd gt 0 then im(bd) = bot
  endif
endif


; *** PLOTTING ***
;if not keyword_set(noplot) then begin
  if not keyword_set(file) then file='ghess'

  if not keyword_set(noplot) then begin
    if keyword_set(save) then begin
      ps_open,file,color=color
      ;if keyword_set(color) then loadct,39,/silent
      ;notop=1

      ; Assigning bottom values to the top
      ; most color (white)
      ;bd = where(im eq min(im),nbd)
      ;gd = where(im ne min(im),ngd)
      ;tmin = min(im(gd))
      ;tmax = max(im(gd))
      ; ncolor = 256.
      ;dim = (tmax-tmin)/ncolor
      ;im(bd) = tmax+dim

    endif else begin
      if (!d.name ne 'PS' and !d.name ne 'Z') then begin
        if keyword_set(color) then begin
          ;loadct,39,/silent
          ;device,decomposed=0
        endif ;else device,decomposed=1
      endif else begin
        ;loadct,39,/silent
      endelse
    endelse
  endif ; not noplot

  ; Using only a certain range, xrange, yrange
  if keyword_set(xrange) then begin
    ;if keyword_set(center) then off=0.5*dx else off=0.    ; center adds 0.5*dx to xarr
    ;dum = closest(xrange[0],xarr-off,ind=xlo)
    ;dum = closest(xrange[1],xarr-off,ind=xhi)
    dum = closest(xrange[0],xarr,ind=xlo)
    dum = closest(xrange[1],xarr,ind=xhi)
    imsz = size(im)
    if imsz[0] eq 3 then im = im(*,xlo:xhi,*) else im=im(xlo:xhi,*)
    xarr = xarr(xlo:xhi)
  endif
  if keyword_set(yrange) then begin
    ;if keyword_set(center) then off=0.5*dy else off=0.    ; center adds 0.5*dx to yarr
    ;dum = closest(yrange[0],yarr-off,ind=ylo)
    ;dum = closest(yrange[1],yarr-off,ind=yhi)
    dum = closest(yrange[0],yarr,ind=ylo)
    dum = closest(yrange[1],yarr,ind=yhi)
    imsz = size(im)
    if imsz[0] eq 3 then im = im(*,*,ylo:yhi) else im=im(*,ylo:yhi)
    yarr = yarr(ylo:yhi)
  endif


  ;Axis titles
  ndec = 1
  if keyword_set(ra) or keyword_set(dec) then if long(year)-float(year) eq 0.0 then ndec=0
  if keyword_set(year) then stryear = stringize(year,/nocomma,ndec=ndec)
  if not keyword_set(xtit) then if keyword_set(xtit) then xtit=strupcase(xtit) else xtit='X'
  ;if strlowcase(xtit) eq 'ra' or strlowcase(xtit) eq 'dec' then xtit=strupcase(xtit)+' (J'+stryear+')'
  if not keyword_set(ytit) then if keyword_set(ytit) then ytit=strupcase(ytit) else ytit='Y'
  ;if strlowcase(ytit) eq 'ra' or strlowcase(ytit) eq 'dec' then ytit=strupcase(ytit)+' (J'+stryear+')'
  if not keyword_set(tit) then begin
    ;if keyword_set(xtit) and keyword_set(ytit) then tit=strlowcase(ytit)+' vs. '+strlowcase(xtit) $
    if keyword_set(xtit) and keyword_set(ytit) then tit=ytit+' vs. '+xtit $
       else tit='Y vs. X'
    if keyword_set(hist) then tit='Number of Points ('+tit+')'
    if keyword_set(total) then if keyword_set(zname) then tit=tit+' (Color as total '+zname+')' $
       else tit=tit+' (Color as total Z)'
    if keyword_set(avg) or keyword_set(mean) then if keyword_set(zname) then tit=tit+' (Color as average '+zname+')' $
       else tit=tit+' (Color as average Z)'
    if keyword_set(std) then if keyword_set(zname) then tit=tit+' (Color as std. dev. of '+zname+')' $
       else tit=tit+' (Color as std. dev. of Z)'
  endif

  ; Plotting it
  if not keyword_set(position) then position = [0.,0.,1.,1.]
  dx1 = position[2]-position[0]
  dy1 = position[3]-position[1]
  x0 = position[0]
  y0 = position[1]
  if not keyword_set(onlyim) then $
  pos = [0.08*dx1+x0,0.10*dy1+y0,0.95*dx1+x0,0.85*dy1+y0]
  if keyword_set(onlyim) then $
  pos = [0.08*dx1+x0,0.10*dy1+y0,0.95*dx1+x0,0.95*dy1+y0]
  if keyword_set(posim) then pos = posim

  ;pos = [0.094*dx1+x0,0.078*dy1+y0,0.97*dx1+x0,0.96*dy1+y0]
  ;position = [ 0.0937550, 0.078130, 0.973443, 0.962896]
  ;pos = [0.08,0.10,0.95,0.85]

  if not keyword_set(noplot) then $
  display,im,xarr,yarr,interpolate=interpolate,log=log,xtit=xtit,$
          ytit=ytit,tit=tit,pos=pos,notop=notop,noerase=noerase,$
          charsize=charsize,xflip=xflip,yflip=yflip,thick=thick,$
          charthick=charthick,framecolor=framecolor,$
          background_color=background_color,color=color,max=top,min=bot,$
          xsize=xsize,ysize=ysize,xtickformat=xtickformat,ytickformat=ytickformat,$
          xstyle=xstyle,ystyle=ystyle,maskvalue=maskvalue,maskcolor=maskcolor

  colpos = [0.08*dx1+x0,0.92*dy1+y0,0.95*dx1+x0,0.95*dy1+y0]
  ;colpos = [0.08,0.92,0.95,0.95]
  if keyword_set(poscol) then colpos = poscol

  ; Overplotting the colorbar
  if not keyword_set(noplot) then begin
  if not keyword_set(nocolorbar) and not keyword_set(onlyim) then begin
    ;if keyword_set(color) then loadct,39,/silent

    if (n_elements(charsize) eq 0) then charsize=1.0
    minim = min(im)
    maxim = max(im)
    if keyword_set(avg) then begin
      minim = min(z)
      maxim = max(z)
      if n_elements(zbot) ne 0 then minim=zbot
      if n_elements(ztop) ne 0 then maxim=ztop
    endif

    ; Setting the format of the colorbar annotation
    if n_elements(format) eq 0 then begin
      len = strlen(strtrim(long(top),2))
      ;len = strlen(strtrim(long(maxim),2))
      if minim lt 0. then len = len+1
      if maxim lt 100 then form = '(G'+strtrim(len+3,2)+'.2)'
      if maxim lt 0.1 then form = '(G8.2)'                  ; room for two sig figs, sign and exponent
      if maxim lt 1.0 then form = '(G7.2)'                  ; room for two sig figs, sign and exponent
      if maxim ge 100 then form = '(I'+strtrim(len,2)+')'
      if maxim gt 1e5 then form = '(G8.2)'                  ; room for two sig figs, sign and exponent
      if keyword_set(log) and minim lt 1.0 then form = '(G'+strtrim(len+4,2)+'.2)'
      div = abs(minim)
      if div eq 0 then div=1
      npow = alog10(abs(maxim)/div)
      if (npow gt 4) then form = '(G8.2)'                   ; for large ranges
    endif else form=format

    ;if not keyword_set(color) then loadct,0,/silent
    if keyword_set(log) and not keyword_set(avg) then begin
      xlog = 1
      minor = 9
      divisions = 0
    endif
    ;if minim eq 0. then minrange = 1. else minrange = minim
    minrange = minim
    maxrange = maxim
    if (n_elements(top) gt 0) then maxrange = top
    if (n_elements(bot) gt 0) then minrange = bot
    ;if keyword_set(mean) or keyword_set(avg) then begin
    if keyword_set(avg) then begin
      minrange = min(z)
      maxrange = max(z)
      if n_elements(zbot) ne 0 then minrange=zbot
      if n_elements(ztop) ne 0 then maxrange=ztop
    endif
    if minim eq 0. and keyword_set(log) then minrange = 1.

    bottom = 10
    ncolors = 245
    if keyword_set(avg) then begin
      bottom=50
      ncolors=200
    endif

    ; Don't want any tickmarks on y-axis
    yorig = !y
    !y.minor = 0

    df_colorbar,position=colpos,minrange=minrange,maxrange=maxrange,$
             charsize=charsize,format=form,bottom=bottom,ncolors=ncolors,$
             xlog=xlog,minor=minor,divisions=divisions,thick=thick,$
             charthick=charthick,color=framecolor,xthick=thick,$
             ythick=thick

    !y = yorig

  endif ; not /nocolorbar
  endif ; not /noplot

if keyword_set(save) then ps_close


; The final cut x/y arrays
xout = x
yout = y
if n_elements(z) gt 0 then zout = z

; Returning originals
if n_elements(orig_xrange) gt 0 then xrange = orig_xrange
if n_elements(orig_yrange) gt 0 then yrange = orig_yrange

; Stopping
if keyword_set(stp) then stop

end
