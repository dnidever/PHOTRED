pro displayc,im,xarr,yarr,color=color,nx=nx,ny=ny,$
          xtit=xtit,ytit=ytit,tit=tit,log=log,interpolate=interpolate,$
          noplot=noplot,save=save,file=file,xrange=xrange,yrange=yrange,$
          stp=stp,maxim=maxim,minim=minim,$
          nocolorbar=nocolorbar,position=position,$
          onlyim=onlyim,noerase=noerase,charsize=charsize,$
          xflip=xflip,yflip=yflip,$
          posim=posim,poscol=poscol,thick=thick,charthick=charthick,$
          framecolor=framecolor,background_color=background_color,$
          colcharsize=colcharsize,xstyle=xstyle,ystyle=ystyle,$
          invert=invert,center=center,xtickformat=xtickformat,$
          ytickformat=ytickformat,xticks=xticks,yticks=yticks,$
          out_posim=out_posim,out_poscol=out_poscol,xsize=xsize,$
          ysize=ysize,aitoff=aitoff,rev=rev,grid=grid,coldivisions=coldivisions,$
          colformat=colformat,zscale=zscale,contrast=contrast,$
          xtickv=xtickv,xtickname=xtickname,ytickv=ytickv,ytickname=ytickname,$
          xminor=xminor,yminor=yminor,xticklen=xticklen,yticklen=yticklen,$
          maskvalue=maskvalue,coltickvalue=coltickvalue,colticknames=colticknames,$
          colxticklen=colxticklen,colright=colright,colvertical=colvertical,$
          colyticklen=colyticklen,maskcolor=maskcolor,notop=notop,top=top,$
          coltitle=coltitle

;+
;
; DISPLAYC.PRO
;
; PURPOSE:
;   To display an image to the screen with a colorbar.
;
;  INPUTS:
;   im            2D image
;   x             Array of x-values
;   y             Array of y-values
;   /color        Plot in color
;   =xtit         Title for x-axis on plot
;   =ytit         Title for y-axis on plot
;   =tit          Overall plot title
;   /log          Plot the logarithm
;   /interp       Interpolate between the points
;   /noplot       Don't plot anything
;   /save         Save the plot to a poscript file
;   =file         Name of postscript file to save the plot to
;   =xrange       X-axis plot range
;   =yrange       Y-axis plot range
;   /stp          Stop at the end of the program
;   =maxim        The largest value to plot
;   =minim        The lowest value to plot
;   /nocolorbar   Don't overplot the colorbar
;   =position     The position to put the plot in, in normalized coordinates
;   /onlyim       Don't plot the colorbar and don't leave space for it.
;   /noerase      Don't erase the screen before you plot the image
;   =charsize     Character size
;   =colcharsize  Character size for the colorbar alone
;   =posim        Position of image (in normalized coordinates)
;   =poscol       Positin of color bar (in normalized coordinates)
;   =thick        The thickness of lines
;   =charthick    The thickness of the annotations
;   =framecolor   The color of the box and annotations
;   =background_color The color of the background
;   /invert       Black on white instead of white on black.
;   /center       The coordinates plotted should be the center of bins (default).
;   /zscale       Use the IRAF zscaling algorithm
;   =contrast     Contrast for the zscaling.  By default zscale=0.25
;   =maskvalue    The value to which bad/missing data is flagged.  Pixels
;                   with this value are set to 0B
;   =maskcolor    The color to use for "bad data".  The default is 0B.
;   =colformat    Format string to use for the colorbar annotations.
;   =coldivisions The number of divisions to use for the colorbar.
;                   The number of colobar annotations are
;                   coldivisions+1
;   =coltickvalue An array of values at which to put the annotations
;                   on the colorbar.
;   =colticknames A string array of names to use for the colorbar annotations
;   =coltitle     Title for colorbar.
;
; OUTPUTS:
;   Image to the screen
;
; USAGE:
;   IDL>displayc,im 
;
; Created by David Nidever October 2005
;-

nim = n_elements(im)

;if (n_params() eq 0) or (nim eq 0) then begin
if (n_params() eq 0) and (n_elements(file) eq 0) then begin
  print,'Syntax displayc,im,x,y,color=color,'
  print,'                xtit=xtit,ytit=ytit,tit=tit,log=log,interpolate=interpolate'
  print,'                noplot=noplot,save=save,file=file,xrange=xrange,yrange=yrange'
  print,'                stp=stp'
  return
endif

; Loading FITS file
if keyword_set(file) then begin
  fil = file_search(file)

  ; File doesn't exist
  if fil(0) eq '' then begin
    print,file,' NOT FOUND'
    return
  endif

  ; Importing the image
  message=''
  FITS_READ,file,im,head,/no_abort,message=message

  ; Fits_read error
  if (message ne '') then begin
    print,file,' ERROR: ',message
    return
  endif

end ; importing FITS file


sz = size(im)
if sz(0) lt 2 or sz(0) gt 3 then begin
  print,'The image is not the right dimension'
  return
endif

if sz(0) eq 2 then begin
  if n_elements(xarr) eq 0 then xarr=findgen(sz(1))
  if n_elements(yarr) eq 0 then yarr=findgen(sz(2))
endif else begin
  if n_elements(xarr) eq 0 then xarr=findgen(sz(2))
  if n_elements(yarr) eq 0 then yarr=findgen(sz(3))
endelse

orig_im = im
orig_xarr = xarr
orig_yarr = yarr

ymin = min(yarr)
ymax = max(yarr)
xmin = min(xarr)
xmax = max(xarr)

if keyword_set(xrange) then begin
  if (xrange(0) lt xmin) or (xrange(1) gt xmax) then begin
    xmin = xrange(0) < xmin
    xmax = xrange(1) > xmax
  endif
endif
if keyword_set(yrange) then begin
  if (yrange(0) lt ymin) or (yrange(1) gt ymax) then begin
    ymin = yrange(0) < ymin
    ymax = yrange(1) > ymax
  endif
endif

;; MAXIM
;if keyword_set(maxim) then begin
;  tbd = where(im gt maxim,ntbd)
;  if ntbd gt 0 then (im)(tbd) = maxim
;endif
;
;; MINIM
;if keyword_set(minim) then begin
;  bbd = where(im lt minim,nbbd)
;  if nbbd gt 0 then (im)(bbd) = minim
;endif


; Plotting
if not keyword_set(noplot) then begin
  if not keyword_set(file) then file='displayc'

  if keyword_set(save) then begin
    ps_open,file,color=color
    ;if keyword_set(color) then loadct,39

  endif else begin
    if (!d.name ne 'PS' and !d.name ne 'Z') then begin
      if keyword_set(color) then begin
        ;loadct,39
        ;device,decomposed=0
      endif ;else device,decomposed=1
    endif else begin
      ;loadct,39
    endelse
  endelse

  ; Using only a certain range
  if keyword_set(xrange) then begin
    dum = closest(xrange(0),xarr,ind=xlo)
    dum = closest(xrange(1),xarr,ind=xhi)
    imsz = size(im)
    if imsz(0) eq 3 then im = im(*,xlo:xhi,*) else im=im(xlo:xhi,*)
    xarr = xarr(xlo:xhi)
  endif
  if keyword_set(yrange) then begin
    dum = closest(yrange(0),yarr,ind=ylo)
    dum = closest(yrange(1),yarr,ind=yhi)
    imsz = size(im)
    if imsz(0) eq 3 then im = im(*,*,ylo:yhi) else im=im(*,ylo:yhi)
    yarr = yarr(ylo:yhi)
  endif

  ; Plotting it
  if not keyword_set(position) then position = [0.,0.,1.,1.]
  dx1 = position(2)-position(0)
  dy1 = position(3)-position(1)
  x0 = position(0)
  y0 = position(1)
  if not keyword_set(onlyim) then $
  pos = [0.08*dx1+x0,0.10*dy1+y0,0.95*dx1+x0,0.85*dy1+y0]
  if keyword_set(onlyim) then $
  pos = [0.08*dx1+x0,0.10*dy1+y0,0.95*dx1+x0,0.95*dy1+y0]
  if keyword_set(posim) then pos = posim

  pim = im
  if keyword_set(invert) then pim=-pim

  display,pim,xarr,yarr,interpolate=interpolate,log=log,xtit=xtit,$
          ytit=ytit,tit=tit,pos=pos,notop=notop,noerase=noerase,$
          charsize=charsize,xflip=xflip,yflip=yflip,thick=thick,$
          charthick=charthick,framecolor=framecolor,$
          background_color=background_color,color=color,$
          xstyle=xstyle,ystyle=ystyle,xtickformat=xtickformat,$
          ytickformat=ytickformat,xticks=xticks,yticks=yticks,$
          max=maxim,min=minim,xsize=xsize,ysize=ysize,$
          aitoff=aitoff,rev=rev,grid=grid,zscale=zscale,contrast=contrast,$
          xtickv=xtickv,xtickname=xtickname,ytickv=ytickv,ytickname=ytickname,$
          xminor=xminor,yminor=yminor,xticklen=xticklen,yticklen=yticklen,$
          maskvalue=maskvalue,maskcolor=maskcolor,top=top

  colpos = [0.08*dx1+x0,0.92*dy1+y0,0.95*dx1+x0,0.95*dy1+y0]
  ;colpos = [0.08,0.92,0.95,0.95]
  if keyword_set(poscol) then colpos = poscol

  ; The actual positions used
  out_posim = pos
  out_poscol = colpos

  ; Overplotting the colorbar
  if not keyword_set(nocolorbar) and not keyword_set(onlyim) then begin
    ;if keyword_set(color) then loadct,39

    ; Log scaling
    ylog = 0
    if keyword_set(log) and not keyword_set(avg) then begin

      xlog = 1
      minor = 9
      divisions = 0
      if keyword_set(colvertical) then begin
        xlog = 0
        ylog = 1
      endif

      len = strlen(strtrim(long(minim),2)) > strlen(strtrim(long(maxim),2))
      if minim lt 1.0 then form = '(G'+strtrim(len+4,2)+'.2)'

      minrange = minim
      if minim eq 0. and keyword_set(log) then minrange = 1.   ; check imgscl.pro, log=0.01
      maxrange= maxim


    ; Linear scaling
    endif else begin

      len = strlen(strtrim(long(minim),2)) > strlen(strtrim(long(maxim),2))
      if minim lt 0.0 then len=len+1    ; need another space for negative sign      

      scal = abs(maxim) > abs(minim)
      frac = abs(maxim-minim)/6.      ; 6 divisions
      pow = alog10(frac)
      npow = round(pow-1)        ; want an extra decimal space
      if npow ge 0 then form = '(I'+strtrim(len,2)+')'
      if npow lt 0 then form = '(F'+strtrim(len+abs(npow)+1,2)+'.'+strtrim(abs(npow),2)+')'

      ; Really large numbers
      if maxim gt 1e5 then form = '(G8.2)'                  ; room for two sig figs, sign and exponent

      minrange = minim
      maxrange = maxim

    endelse


    ; Scale of colors
    if keyword_set(avg) then begin
      bottom=50
      ncolors=200
    endif else begin
      ; display/imgscl reserves 0 and 255 for annotations
      bottom=1       ;10
      ncolors=254    ;245
      ;bottom=0       ;10
      ;ncolors=255    ;245
    endelse

    if keyword_set(colformat) then form=colformat
    if keyword_set(coldivisions) then divisions=coldivisions
    if not keyword_set(charsize) then charsize=1.0
    if not keyword_set(colcharsize) then colcharsize=charsize

    ; Displaying the colorbar
    df_colorbar,position=colpos,minrange=minrange,maxrange=maxrange,$
               charsize=colcharsize,format=form,bottom=bottom,ncolors=ncolors,$
               xlog=xlog,ylog=ylog,minor=minor,divisions=divisions,thick=thick,$
               charthick=charthick,color=framecolor,xthick=thick,$
               ythick=thick,invertcolor=invert,xtickv=coltickvalue,$
               ticknames=colticknames,xticklen=colxticklen,yticklen=colyticklen,$
               vertical=colvertical,right=colright,title=coltitle

  endif ; not /nocolorbar

  if keyword_set(save) then ps_close

endif

if keyword_set(stp) then stop

im = orig_im
xarr = orig_xarr
yarr = orig_yarr

end
