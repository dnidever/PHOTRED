;+
; NAME:
;	DISPLAY
;
; PURPOSE:
;	This procedure will display an image with the TV command that fills
;	the plotting window.  It handles scale, annotations, X and PostScript
;	devices, aspect ratios, logarithmic scaling, and interpolation.  The
;	first colormap entry is reserved for the background (pixels flagged
;	with the MASKVALUE value are mapped to this color) and the last entry
;	is reserved for user defined colored annotations.  The annotation
;	plotted by this procedure are in the color !P.Color.
;
; CATEGORY:
;	Image display.
;
; CALLING SEQUENCE:
;	DISPLAY, Image, XS, YS
;
; INPUTS:
;	Image:	Two-dimensional array to be displayed.
;
; OPTIONAL INPUTS:
;	XS:	Vector of x-axis values.  The length must equal the number of
;		rows in <Image>
;
;	YS:	Vector of y-axis values.  The length must equal the number of
;		columns in <Image>
;
; KEYWORD PARAMETERS:
;	TITLE=	Set this keyword to a string containing the title annotation
;		to be used by PLOT.
;
;	XTITLE=	Set this keyword to a string containing the x-axis annotation
;		to be used by PLOT.
;
;	YTITLE=	Set this keyword to a string containing the y-axis annotation
;		to be used by PLOT.
;
;	ASPECT=	Set this keyword to the aspect ratio (width/height) of the
;		pixels.  /ASPECT is the same as ASPECT=1 and produces square
;		pixels.
;
;	/INTERPOLATE:
;		Set this switch to enable bilinear interpolation for pixels
;		in the expanded image.  See /PS_FINE for information
;		on using this switch on a PostScript device.
;
;	/HIST	Set this for histogram equalization colormap scaling
;
;	MASKVALUE=
;		Set this keyword to the value that pixels with bad data or
;		no data have been flagged with.  These will be mapped to 0B.
;
;       MASKCOLOR=
;               The color to use for "bad data".  The default is 0B.  
;
;	MIN=	The minimum value of <Image> to be considered.  If MIN is not
;		provided, <Image> is searched for its minimum value.  All
;		values less than or equal to MIN are set to 1 in the Result.
;
;	MAX=	The maximum value of <Image> to be considered.  If MAX is not
;		provided, <Image> is searched for its maximum value.  All
;		values greater than or equal to MAX are set to TOP in the
;		Result.
;
;	TOP=	The maximum value of the scaled result.  If TOP is not
;		specified, 255 is used. Note that the minimum value of the
;		scaled result is always 1 (NOT 0 as in BYTSCL).
;
;	LEVELS=	Set this keyword to a vector of data value boundaries between
;		which all elements of <Image> have the same scaled byte
;		value.  e.g. LEVELS=[0,1,2,5] maps all values below 0 and
;		above 5 to 0B, map values between 0 and 1 to 1B, map values
;		between 1 and 2 to 128B, and map values between 2 and 5 to
;		255B.  This does not plot contours.
;
;	/LOG:	Set this switch to cause a logarithmic mapping.  This is
;		overridden by the LEVELS keyword.
;
;	/PS_FINE:
;		Set the switch to enable higher resolution images on a
;		PostScript device.  This is only useful with /INTERPOLATE and
;		will increase the size of the PostScript file.
;
;	/NOERASE:
;		Set the switch to prevent output device from being erased
;		before the image, scales, and annotations are displayed.
;
;	/NO_EXPAND:
;		Set this switch to prevent the image from being expanded
;		to fill the plotting window.  Scaling to byte type is still
;		performed.
;
;       POSITION=
;               This sets the position of the plot window in
;               normalized coordinates.  Lower-left and upper-right coordinates.  
;               POSITION = [x_lower_left, y_lower,left, x_upper_right,y_upper_right]
;               e.g. POSITION = [0.1, 0.1, 0.9, 0.9]
;               The default is set by IDL but normally is (with charsize=1.0)
;               position = [ 0.0937550, 0.078130, 0.973443, 0.962896]
;
;       /DECOMP
;               Changes the device composed setting from 1 to 0 for
;               the TV display.  Causes it to be in color on certain windows.
;
;       THICK =  
;              The thickness of the lines that make up the box.
; 
;       CHARTHICK=
;              The thickness of the annotations.
;              
;       BACKGROUND_COLOR=
;              The color of the background (0 for black, !p.color for white).
;              Default: PS,white; other,black.
;
;       CENTER=
;              When interpolating the pixels are assumed to lie at the midpoint of their
;              coordinates rather than at their lower-left corner.  By default center=1
;
;       /COLOR:
;              Set this keyword when making color PS plots of Truecolor images.
;
;       /XFLIP: 
;              Reverses the image and axes in the x-direction.  The same thing
;              can be accomplished by reversing the order in xrange (decreasing
;              instead of increasing).
;
;       /YFLIP: 
;              Reverses the image and axes in the y-direction.  The same thing
;              can be accomplished by reversing the order in xrange (decreasing
;              instead of increasing).
;
;      /SQUARE:
;              Create a square window.
;
;      /GRID:
;              Plot a grid over the aitoff or orthographic projection.
;
;      XSIZE =
;              Set the number of pixels (in the x-dimension) in the output image.
;              Only for PS output.
;
;      YSIZE =
;              Set the number of pixels (in the y-dimension) in the output image
;              Only for PS output.
;      FILE =
;              Read the image from this FITS file.
;      /RESCALE
;              Rescale before taking the log.  This is the way it was
;              originally done (in imgscl.pro).
;
;      /ZSCALE
;              Use IRAF zscaling.
;      CONTRAST =
;              Sets the zscaling contrast.  The default is
;              contrast=0.25
;
;	TOP=
;              The maximum value of the scaled result.  If TOP is not
;              specified, 255 is used.
;
; SIDE EFFECTS:
;	TV display is altered.
;
; RESTRICTIONS:
;	This routine may work for other devices, but it has only been tested
;	on 'X' and 'PS'.
;
; PROCEDURE:
;	Straight forward.  :-)
;
; EXAMPLE:
;	LoadCT, 3
;	image = SHIFT(DIST(20, 20), 10, 10)
;	scale = FINDGEN(20) - 10.0
;	DISPLAY, image, scale, scale, /INTERPOLATE, TITLE='!6Smooth Slope', $
;		/ASPECT
;	;Use CONTOUR with /OVERPLOT to overlay contours.
;	CONTOUR, image, scale, scale, LEVELS=1.0+FINDGEN(4)*2.0, /OVERPLOT
;
;	DISPLAY		;prints out a "Usage:" line
;
; MODIFICATION HISTORY:
; 	Written by:	Fen Tamanaha, July 10, 1993.  Release 3.1
;	July 13, 1993	Fen: (3.2) Fixed /No_Expand
;	July 16, 1993	Fen: (3.3) Really fixed /No_Expand
;       March 13, 2005  D. Nidever, added POSITION keyword
;
; ** NOTE: **
;       When using the interpolate option the scale is exact, while
;       the scale is extended 0.5 pixels in every direction if the
;       interpolate option is not set.  
;-


pro	Display, image, $
		xs, ys, $
		Title=t, XTitle=xt, YTitle=yt, $
		MIN=minval, MAX=maxval, $
		LOG=log_scaling, $
		LEVELS=l, HIST=hist,$
		ASPECT=aspect, $
		INTERPOLATE=interp, $
		MASKVALUE=maskvalue, $
                maskcolor=maskcolor, $
		PSFINE=psfine, $
		NO_EXPAND=no_expand, $
		NOERASE=noerase, $
		HELP=help,$
                charsize=charsize,$
                noframe=noframe,$
                xstyle=xstyle,$
                ystyle=ystyle,$
                xticks=xticks,$
                yticks=yticks,$
                framecolor=framecolor,$
                position=position,$
                decomp=decomp,$
                xrange=xrange,yrange=yrange,$
                notop=notop,xflip=xflip,yflip=yflip,$
                stp=stp,thick=thick,charthick=charthick,$
                background_color=background_color,$
                color=color,center=center,aitoff=aitoff,$
                iso=iso,orthographic=orthographic,$
                p0lat=p0lat,p0lon=p0lon,rotation=rotation,$
                rev=rev,xlog=xlog,ylog=ylog,square=square,$
                xminor=xminor,yminor=yminor,xticklen=xticklen,$
                yticklen=yticklen,xtickformat=xtickformat,$
                ytickformat=ytickformat,grid=grid,$
                xsize=xsize,ysize=ysize,file=file,$
                xtickname=xtickname,xtickv=xtickv,ytickname=ytickname,$
                ytickv=ytickv,imout=imout,rescale=rescale,$
                zscale=zscale,contrast=contrast,top=top


    On_Error, 2

if !d.name eq 'X' and keyword_set(square) then begin
 device,window_state = opnd                    ;Get state of each window
 free = where(opnd eq 0,nfree)
 efree = min(where(free mod 4 eq 0))
 iplane = free(efree)
 window,iplane,xsize=640,ysize=600
end

if n_elements(framecolor) eq 0 then framecolor=!p.color

; Loading FITS file
if keyword_set(file) then begin

  LOADINPUT,file,inpfiles,count=ninpfiles

  ; Multiple FITS files input
  if ninpfiles gt 1 then begin

    for f=0,ninpfiles-1 do begin

      print,'Displaying ',inpfiles[f]

      Display, Title=inpfiles[f], XTitle=xt, YTitle=yt, $
		MIN=minval, MAX=maxval, $
		LOG=log_scaling, $
		LEVELS=l, HIST=hist,$
		ASPECT=aspect, $
		INTERPOLATE=interp, $
		MASKVALUE=maskvalue, $
		maskcolor=maskcolor, $
		PSFINE=psfine, $
		NO_EXPAND=no_expand, $
		NOERASE=noerase, $
		HELP=help,$
                charsize=charsize,$
                noframe=noframe,$
                xstyle=xstyle,$
                ystyle=ystyle,$
                xticks=xticks,$
                yticks=yticks,$
                framecolor=framecolor,$
                position=position,$
                decomp=decomp,$
                xrange=xrange,yrange=yrange,$
                notop=notop,xflip=xflip,yflip=yflip,$
                stp=stp,thick=thick,charthick=charthick,$
                background_color=background_color,$
                color=color,center=center,aitoff=aitoff,$
                iso=iso,orthographic=orthographic,$
                p0lat=p0lat,p0lon=p0lon,rotation=rotation,$
                rev=rev,xlog=xlog,ylog=ylog,square=square,$
                xminor=xminor,yminor=yminor,xticklen=xticklen,$
                yticklen=yticklen,xtickformat=xtickformat,$
                ytickformat=ytickformat,grid=grid,$
                xsize=xsize,ysize=ysize,file=inpfiles[f],$
                xtickname=xtickname,xtickv=xtickv,ytickname=ytickname,$
                ytickv=ytickv,imout=imout,rescale=rescale,$
                zscale=zscale,contrast=contrast

      line=''
      read,'Hit Enter to go on, q to quit: ',line
      if strmid(strlowcase(strtrim(line,2)),0,1) eq 'q' then return

    end ; file loop

    return

  endif  ; multiple files input

  ; Only one file input
  fil = file_search(file)

  ; File doesn't exist
  if fil(0) eq '' then begin
    print,file,' NOT FOUND'
    return
  endif

  ; Importing the image
  message=''
  FITS_READ,file,image,head,/no_abort,message=message

  ; Fits_read error
  if (message ne '') then begin
    print,file,' ERROR: ',message
    return
  endif

end ; importing FITS file


;
; Validate arguments and get size of the image
;

    nim = n_elements(image)
    if nim eq 0 then goto,USAGE
    sz = Size(image)

    ; NOT three dimensional image
    if (sz[0] ne 3) then begin

      nparms = N_Params()
      If ( Keyword_Set(help) ) Then nparms = 0	;force a "Usage:" line
      Case ( nparms ) Of
          1: Begin
              sz = Size(image)
              If ( sz[0] NE 2 ) Then Begin
                  Message, '<image> must be an array.'
              EndIf
              xs = FIndGen(sz[1])
              ys = FIndGen(sz[2])
          End
          2: Begin
              sz = Size(image)
              If ( sz[0] NE 2 ) Then Begin
                  Message, '<image> must be an array.'
              EndIf
              If ( N_Elements(xs) NE sz[1] ) Then Begin
                  Message, '<xs> does not match <image> dimensions.'
              EndIf
              ys = FIndGen(sz[2])
          End
          3: Begin
              sz = Size(image)
              If ( sz[0] NE 2 ) Then Begin
                  Message, '<image> must be an array.'
              EndIf
              If ( N_Elements(xs) NE sz[1] ) Then Begin
                  Message, '<xs> does not match <image> dimensions.'
              EndIf
              If ( N_Elements(ys) NE sz[2] ) Then Begin
                  Message, '<ys> does not match <image> dimensions.'
              EndIf
          End
          Else: Begin
          ;    USAGE:
          ;     Message, 'Usage: DISPLAY, image [,xs [,ys]] [,TITLE=] [,XTITLE=] [,YTITLE=]', /Info
          ;    Message, '           [,MIN=] [,MAX=] [,/LOG] [,LEVELS=] [/HIST]', /Info
          ;    Message, '           [,ASPECT=] [,/INTERPOLATE] [MASKVALUE=]', /Info
	  ;    Message, '           [,/NO_EXPAND] [,/NOERASE] [,/PSFINE]', /Info
          ;    Return
          End
        Endcase

        ; Not enough inputs
        if n_params() eq 0 and n_elements(file) eq 0 then begin
          USAGE:
          Message, 'Usage: DISPLAY, image [,xs [,ys]] [,TITLE=] [,XTITLE=] [,YTITLE=]', /Info
          Message, '           [,MIN=] [,MAX=] [,/LOG] [,LEVELS=] [/HIST]', /Info
          Message, '           [,ASPECT=] [,/INTERPOLATE] [MASKVALUE=]', /Info
          Message, '           [,/NO_EXPAND] [,/NOERASE] [,/PSFINE]', /Info
          Return
        endif

    ; Three dimensional image
    endif else begin

      if sz[1] eq 3 then true=1
      if sz[2] eq 3 then true=2
      if sz[3] eq 3 then true=3
      if not keyword_set(true) then begin
        print,'If IMAGE is 3D then one dim must have 3 elements'
        return
      endif

    endelse      

     ; SETTING THE AXIS RANGES
     If (n_elements(xrange) ne 0) then old_xrange = xrange
     If (n_elements(yrange) ne 0) then old_yrange = yrange
     If (not keyword_set(xrange) or not keyword_set(yrange)) then begin

       nx = sz[1]
       ny = sz[2]

       if sz[0] eq 3 then begin
         case true of
           1: begin
                nx = sz[2] & ny = sz[3]
              end
           2: begin
                nx = sz[1] & ny = sz[3]
              end
           3: begin
                nx = sz[1] & ny = sz[2]
              end
         endcase
       endif

       ; x and y arrays NOT input
       xrange = [-0.5, Float(nx)-0.5]    ;plus 0.5 on both ends
       yrange = [-0.5, Float(ny)-0.5]

       ; x and y array input
       ; NO interpolation, 1/2 pixel offset
       if not keyword_set(interp) then begin
         if keyword_set(xs) then begin
           dx = xs[1]-xs[0]
           xrange = [min(xs)-dx*0.5,max(xs)+dx*0.5]
         endif
         if keyword_set(ys) then begin
           dy = ys[1]-ys[0]
           yrange = [min(ys)-dy*0.5,max(ys)+dy*0.5]
         endif
       endif

       ; interpolation, NO 1/2 pixel offset
       if keyword_set(interp) then begin
         xrange = [0, Float(nx)-1] 
         yrange = [0, Float(ny)-1]

         if keyword_set(xs) then xrange = [min(xs),max(xs)]
         if keyword_set(ys) then yrange = [min(ys),max(ys)]
       endif
     Endif   ; xrange not input
     if (n_elements(old_xrange) ne 0) then xrange = old_xrange
     if (n_elements(old_yrange) ne 0) then yrange = old_yrange


; Using ZSCALE
if keyword_set(zscale) and size(image,/n_dim) eq 2 then begin
  if n_elements(contrast) eq 0 then contrast=0.25
  zscale,image,z1,z2,contrast=contrast
  minval = z1
  maxval = z2
endif
if keyword_set(zscale) and size(image,/n_dim) ne 2 then $
  print,'ZSCALE only works with 2D images'


;
; The plotting device must be erased to reset the system variables so that
;	IMGEXP will get the default values.  The /NOERASE keyword should
;	be used to prevent this.  One typical situation is when DISPLAY
;	is called after a !P.MULTI change.  An ERASE at this point would
;	destroy the previous plots.
;
    If ( Not Keyword_Set(noerase) ) Then Begin
        if not keyword_set(background_color) then $
          if (!d.name eq 'PS') then background_color=!p.color else background_color=0
	Erase,background_color
    EndIf

;
; If /PSFINE is set then up the intermediate interpolated image width.
;	This only has an effect on PostScript output.
;
    If (Keyword_Set(psfine) ) Then Begin
	psis = 512.0
    Endif


    ; Use the default IDL plot window settings, from !x.window, !y.window
    if not keyword_set(position) then begin
      ;position = [ 0.0937550, 0.078130, 0.973443, 0.962896]

      old_set = !d.name
      set_plot,'X'

      xr = xrange
      yr = yrange
      if keyword_set(xflip) then xr = reverse(xr)
      if keyword_set(yflip) then yr = reverse(yr)

      ; Make a dummy plot
      Plot, [0,1], /NoErase, /NoData, XStyle=4, YStyle=4, $
                  XRange=xr, YRange=yr, $
                  Title=title, XTitle=xtitle, YTitle=ytitle,$
                  charsize=charsize,xticks=xticks,yticks=yticks,$
                  color=0,xthick=thick,ythick=thick,$
                  charthick=charthick,xlog=xlog,ylog=ylog

      set_plot,old_set

      position = [!x.window(0),!y.window(0),!x.window(1),!y.window(1)]

    endif


 ;  ;Using the position commands
 ;  if keyword_set(position) then begin

     scalable = (!D.Flags And 1) NE 0   ; PS has scalable pixels

     ; For postscript 
     if keyword_set(scalable) then begin
       if (sz[0] eq 2) then begin
         if not keyword_set(xsize) then $
           xsize = 20*Float(sz[1]) < 600.   ; it can't be too large
         if not keyword_set(ysize) then $
           ysize = 20*Float(sz[2]) < 600.
       endif else begin
         if not keyword_set(xsize) then $
           xsize = 20*Float(sz[2]) < 600.   ; it can't be too large
         if not keyword_set(ysize) then $
           ysize = 20*Float(sz[3]) < 600.
       endelse
       
       dev_pos = position

       device = 0
       normal = 1              ; use normal coordinates (dev_pos)

     ; Normal X-windows
     endif else begin

       ; for X-Windows
       xsize = (position[2] - position[0]) * !D.X_VSIZE
       ysize = (position[3] - position[1]) * !D.Y_VSIZE
       xstart = position[0] * !D.X_VSIZE
       ystart = position[1] * !D.Y_VSIZE

       dev_pos = [xstart,ystart,xstart+xsize,ystart+ysize]

       normal = 0
       device = 1                 ; use device coordinates (dev_pos)

     endelse


     ; The points are at the center of the pixel, not the lower-left corner
     if not keyword_set(center) then center=1


     ; INTERPOLATING THE IMAGE
     if sz[0] eq 2 then $
       im = congrid(image,xsize,ysize,interp=interp,center=center)  

     if sz[0] eq 3 then begin
       case true of
         1: begin
              im = bytarr(3,xsize,ysize)
              for i=0,2 do im[i,*,*] = congrid(reform(image[i,*,*]),xsize,ysize,interp=interp,center=center)  
            end
         2: begin
              im = bytarr(xsize,3,ysize)
              for i=0,2 do im[*,i,*] = congrid(reform(image[*,i,*]),xsize,ysize,interp=interp,center=center)  
            end
         3: begin
              im = bytarr(xsize,ysize,3)
              for i=0,2 do im[*,*,i] = congrid(reform(image[*,*,i]),xsize,ysize,interp=interp,center=center)  
            end
       endcase
     endif

  ; ; normal setup
  ; endif else begin
  ;
  ;   im = ImgExp(image, xs, ys, xscale, yscale, xrange, yrange, $
  ;		Aspect=aspect, Interpolate=Keyword_Set(interp), $
  ;		MaskValue=maskvalue, Position=dev_pos, PS_Interp_Size=psis, $
  ;		No_Expand=Keyword_Set(no_expand))
  ;
  ;   normal = 0
  ;   device = 1                  ; use device coordinates (dev_pos)
  ; endelse



   sz2 = Size(im)
   im_x_width = Float(sz2[1])                   ;image width
   im_y_width = Float(sz2[2])                   ;image height

   if sz(0) eq 3 then begin
      case true of
        1: begin
             im_x_width = Float(sz2[2])
             im_y_width = Float(sz2[3])
           end
        2: begin
             im_x_width = Float(sz2[1])
             im_y_width = Float(sz2[3])
           end
        3: begin
             im_x_width = Float(sz2[1])
             im_y_width = Float(sz2[2])
           end
      endcase
   endif

;
; Determine the device coordinates of the plotting regions.
;
    dev_x_width = dev_pos[2] - dev_pos[0] + 1
    dev_y_width = dev_pos[3] - dev_pos[1] + 1
;    If ( (im_x_width GT dev_x_width) Or (im_y_width GT dev_y_width) ) Then Begin
;	Message, 'Error: Scaled image is larger than plotting window.'
;    EndIf

;
; Convert a non-byte type image to byte with IMGSCL.  The bottom entry
;	of the color table is reserved for the background/NODATA color
;	by IMGSCL.  The top color table entry will also be reserved
;	here for annotation color.
;
    If ( sz2[sz2[0]+1] GT 1 ) Then Begin
    	if keyword_set(hist) then byte_im=hist_equal(im) else $
        if n_elements(top) eq 0 then top=!D.Table_Size-2
        if keyword_set(notop) then top=!D.Table_Size-1

  	byte_im = ImgScl(im, Min=minval, Max=maxval, Top=top, $
                         Log=log_scaling, Levels=l, MaskValue=maskvalue, maskcolor=maskcolor, rescale=rescale)
        ;if not keyword_set(notop) then begin
  	;  byte_im = ImgScl(im, Min=minval, Max=maxval, Top=!D.Table_Size-2, $
	;	  	  Log=log_scaling, Levels=l, MaskValue=maskvalue, maskcolor=maskcolor, rescale=rescale)
        ;endif else begin
  	;  byte_im = ImgScl(im, Min=minval, Max=maxval, Top=!D.Table_Size-1, $
	;	  	  Log=log_scaling, Levels=l, MaskValue=maskvalue, maskcolor=maskcolor, rescale=rescale)   
        ;endelse
    EndIf Else Begin
	Message, '<Image> is already byte type. No scaling done.', /Info
	byte_im = im
    Endelse

if (n_elements(xflip) eq 0) then if (xrange(1) lt xrange(0)) then xflip=1 else xflip=0
if (n_elements(yflip) eq 0) then if (yrange(1) lt yrange(0)) then yflip=1 else yflip=0

if (size(byte_im))(0) eq 2 then begin
  if keyword_set(xflip) then byte_im = reverse(byte_im,1)
  if keyword_set(yflip) then byte_im = reverse(byte_im,2)
endif
if (size(byte_im))(0) eq 3 then begin
  if keyword_set(xflip) then byte_im = reverse(byte_im,2)
  if keyword_set(yflip) then byte_im = reverse(byte_im,3)
endif


; Different Projection
if keyword_set(orthographic) or keyword_set(aitoff) then begin
   if not keyword_set(P0lat) then P0lat=0
   if not keyword_set(P0lon) then P0lon=180      ; lon=180 in the middle
   if not keyword_set(rot) then rot=0

   ; These are the limits to use, l=[0,360], lat=[-90,+90]
   limit=[-90,0,90,360]

   map_set,P0lat,P0lon,rot,aitoff=aitoff,orthographic=orthographic,$
           isotropic=isotropic,reverse=rev,limit=limit,lonlab=lonlab,$
           latlab=latlab,position=position

   if !d.name eq 'PS' then scale=0.03 else scale=0.5

   ; 3D image
   if sz2[0] eq 3 then begin
     ;new1 = map_image(reform(byte_im(0,*,*)),sx,sy,xsiz,ysiz,/bilin,compress=1,scale=0.5)
     ;new2 = map_image(reform(byte_im(1,*,*)),sx,sy,xsiz,ysiz,/bilin,compress=1,scale=0.5)
     ;new3 = map_image(reform(byte_im(2,*,*)),sx,sy,xsiz,ysiz,/bilin,compress=1,scale=0.5)
     new1 = map_image(reform(byte_im[0,*,*]),sx,sy,xsiz,ysiz,/bilin,compress=1,scale=scale,latmin=-90,lonmin=0,$
                     latmax=90,lonmax=360)
     new2 = map_image(reform(byte_im[1,*,*]),sx,sy,xsiz,ysiz,/bilin,compress=1,scale=scale,latmin=-90,lonmin=0,$
                     latmax=90,lonmax=360)
     new3 = map_image(reform(byte_im[2,*,*]),sx,sy,xsiz,ysiz,/bilin,compress=1,scale=scale,latmin=-90,lonmin=0,$
                     latmax=90,lonmax=360)
     sz = size(new1)
     new = bytarr(3,sz[1],sz[2])
     new[0,*,*] = new1
     new[1,*,*] = new2
     new[2,*,*] = new3
     byte_im = new

   ; 2D image
   endif else begin
     ;new = map_image(byte_im,sx,sy,xsiz,ysiz,/bilin,compress=1,scale=0.5)
     new = map_image(byte_im,sx,sy,xsiz,ysiz,/bilin,compress=1,scale=scale,latmin=-90,lonmin=0,$
                     latmax=90,lonmax=360)
     byte_im = new
   endelse

   ;stop

   dev_pos = [sx,sy,xsiz+sx,ysiz+sy]
   normal=0
   device=1
   noframe=1
endif


;
; Put the image on the TV.
;
;    TV, byte_im, /Device, dev_pos(0), dev_pos(1), $
;		XSize=dev_pos(2)-dev_pos(0), YSize=dev_pos(3)-dev_pos(1)
    if keyword_set(decomp) then device,decomposed=0
    if (!d.name ne 'PS') or ((size(byte_im))(0) lt 3) then begin
      TV, byte_im, device=device, normal=normal, dev_pos[0], dev_pos[1], $
	  	XSize=dev_pos[2]-dev_pos[0], YSize=dev_pos[3]-dev_pos[1],true=true
    endif else begin
      TV_24, byte_im, device=device, normal=normal, dev_pos[0], dev_pos[1], $
             XSize=dev_pos[2]-dev_pos[0], YSize=dev_pos[3]-dev_pos[1]
    endelse
    if keyword_set(decomp) then device,decomposed=1

imout = byte_im

;
; Manage the title and axis labels.
;
    If ( Keyword_Set(t) ) Then Begin
        title = String(t)
    EndIf Else Begin
        title = ' '
    EndElse
 
    If ( Keyword_Set(xt) ) Then Begin
        xtitle = String(xt)
    EndIf Else Begin
        xtitle = ' '
    EndElse
 
    If ( Keyword_Set(yt) ) Then Begin
        ytitle = String(yt)
    EndIf Else Begin
        ytitle = ' '
    Endelse


; Overplotting the grid
if (keyword_set(orthographic) or keyword_set(aitoff)) and $
   keyword_set(grid) then begin
   ;if (!d.name eq 'PS') and keyword_set(color) then loadct,39

   map_grid,/label,co=framecolor,/horizon

   xx = dev_pos[0] + 0.5 * (dev_pos[2]-dev_pos[0])
   yy = dev_pos[3] + 0.01 * (dev_pos[3]-dev_pos[1])

   ; Overplotting the title
   if !d.name eq 'PS' then titco = 0. else titco=framecolor
   xyouts,xx,yy,title,/device,charsize=charsize,charthick=charthick,align=0.5,co=titco

   ;noframe = 1  ; don't plot the frame
end


; Overplot annotations.
;
    if not keyword_set(noframe) then begin

      if not keyword_set(xstyle) then xstyle = 1
      if not keyword_set(ystyle) then ystyle = 1

      ;if (!d.name eq 'PS') and keyword_set(color) then loadct,39

      xr = xrange
      yr = yrange
      if keyword_set(xflip) then xr = reverse(xr)
      if keyword_set(yflip) then yr = reverse(yr)

      ;  Plot, [0,1], /NoErase, /NoData, XStyle=1, YStyle=1, $
      Plot, [0,1], /NoErase, /NoData, XStyle=xstyle, YStyle=ystyle, $
                  device=device, normal=normal, Position=dev_pos, $
  ;                  /Device, Position=dev_pos, $
                  XRange=xr, YRange=yr, $
                  Title=title, XTitle=xtitle, YTitle=ytitle,$
                  charsize=charsize,xticks=xticks,yticks=yticks,$
                  color=framecolor,xthick=thick,ythick=thick,$
                  charthick=charthick,xlog=xlog,ylog=ylog,xminor=xminor,$
                  yminor=yminor,yticklen=yticklen,xticklen=xticklen,$
                  xtickformat=xtickformat,ytickformat=ytickformat,$
                  xtickname=xtickname,xtickv=xtickv,ytickname=ytickname,$
                  ytickv=ytickv

    endif ; not noframe

    if keyword_set(stp) then stop

    Return
End


