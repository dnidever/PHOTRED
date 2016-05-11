pro srcpush,arr1,arr2,count=count

count=0

narr1 = n_elements(arr1)
narr2 = n_elements(arr2)

; ARR1 already exists
if (narr1 gt 0) then begin
  if narr2 gt 0 then $
    arr1 = [temporary(arr1),arr2]

; ARR1 does NOT exist
endif else begin
  if narr2 gt 0 then arr1 = arr2
endelse

; The number of final elements
count = n_elements(narr1)

end

;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro undef,varname
 n = n_elements(varname)
 if (n gt 0) then tempvar = SIZE(TEMPORARY(varname))
end

;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

;+
; NAME:
;       HIST_ND
;
; PURPOSE:
;
;       Perform an N-dimensional histogram, also known as the joint
;       density function of N variables, ala HIST_2D.
;
; CALLING SEQUENCE:
;       hist=HIST_ND(V,BINSIZE,MIN=,MAX=,NBINS=,REVERSE_INDICES=)
;
; INPUTS:
;
;       V: A NxP array representing P data points in N dimensions.  
;
;       BINSIZE: The size of the bin to use. Either an N point vector
;         specifying a separate size for each dimension, or a scalar,
;         which will be used for all dimensions.  If BINSIZE is not
;         passed, NBINS must be.
;
; OPTIONAL INPUTS: 
;
;       MIN: The minimum value for the histogram.  Either a P point
;         vector specifying a separate minimum for each dimension, or
;         a scalar, which will be used for all dimensions.  If
;         omitted, the natural minimum within the dataset will be
;         used.
;
;       MAX: The maximum value for the histogram.  Either a P point
;         vector specifying a separate maximmum for each dimension, or
;         a scalar, which will be used for all dimensions. If omitted,
;         the natural maximum within the dataset will be used.
;
;       NBINS: Rather than specifying the binsize, you can pass NBINS,
;         the number of bins in each dimension, which can be a P point
;         vector, or a scalar.  If BINSIZE it also passed, NBINS will
;         be ignored, otherwise BINSIZE will then be calculated as
;         binsize=(max-min)/nbins.  Note that *unlike* RSI's version
;         of histogram as of IDL 5.4, this keyword actually works as
;         advertised, giving you NBINS bins over the range min to max.
;
; KEYWORD PARAMETERS:
;       
;       MIN,MAX,NBINS: See above
;       
;       REVERSE_INDICES: Set to a named variable to receive the
;         reverse indices, for mapping which points occurred in a
;         given bin.  Note that this is a 1-dimensional reverse index
;         vector (see HISTOGRAM).  E.g., to find the indices of points
;         which fell in a histogram bin [i,j,k], look up:
;
;             ind=[i+nx*(j+ny*k)]
;             ri[ri[ind]:ri[ind+1]-1]
;
;         See also ARRAY_INDICES for converting in the other
;         direction.
;
; OUTPUTS:
;
;       hist: The N-Dimensional histogram, of size N1xN2xN3x...xND
;         where the Ni's are the number of bins implied by the data,
;         and/or optional inputs min, max and binsize.
;
; OPTIONAL OUTPUTS:
;
;       The reverse indices.
;
; EXAMPLE:
;       
;       v=randomu(sd,3,100)
;       h=hist_nd(v,.25,MIN=0,MAX=1,REVERSE_INDICES=ri)
;
; SEE ALSO:
;
;       HISTOGRAM, HIST_2D
;
; MODIFICATION HISTORY:
;
;       Mon Mar 5 09:45:53 2007, J.D. Smith <jdsmith@as.arizona.edu>
;
;               Correctly trim out of range elements from the
;               histogram, when MIN/MAX are specified. Requires IDL
;               v6.1 or later.
;
;       Tue Aug 19 09:13:43 2003, J.D. Smith <jdsmith@as.arizona.edu>
;
;               Slight update to BINSIZE logic to provide consistency
;               with HIST_2D.
;
;       Fri Oct 11 10:10:01 2002, J.D. Smith <jdsmith@as.arizona.edu>
;
;               Updated to use new DIMENSION keyword to MAX/MIN.
;
;       Fri Apr 20 12:57:34 2001, JD Smith <jdsmith@astro.cornell.edu>
;
;               Slight update to NBINS logic.  More aggressive keyword
;               checking.
;
;       Wed Mar 28 19:41:10 2001, JD Smith <jdsmith@astro.cornell.edu>
;
;               Written, based on HIST_2D, and suggestions of CM.
;
;-
;##############################################################################
;
; LICENSE
;
;  Copyright (C) 2001-2003, 2004, 2007 J.D. Smith
;
;  This file is free software; you can redistribute it and/or modify
;  it under the terms of the GNU General Public License as published
;  by the Free Software Foundation; either version 2, or (at your
;  option) any later version.
;
;  This file is distributed in the hope that it will be useful, but
;  WITHOUT ANY WARRANTY; without even the implied warranty of
;  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
;  General Public License for more details.
;
;  You should have received a copy of the GNU General Public License
;  along with this file; see the file COPYING.  If not, write to the
;  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
;  Boston, MA 02110-1301, USA.
;
;##############################################################################

pro srcfilter_hist_nd,V,bs,ret,MIN=mn,MAX=mx,NBINS=nbins,REVERSE_INDICES=ri
  s=size(V,/DIMENSIONS)
  if n_elements(s) ne 2 then message,'Input must be N (dimensions) x P (points)'
  if s[0] gt 8 then message, 'Only up to 8 dimensions allowed'
  
  imx=max(V,DIMENSION=2,MIN=imn)
  
  if n_elements(mx) eq 0 then mx=imx
  if n_elements(mn) eq 0 then mn=imn
  
  if s[0] gt 1 then begin 
     if n_elements(mn)    eq 1 then mn=replicate(mn,s[0])
     if n_elements(mx)    eq 1 then mx=replicate(mx,s[0])
     if n_elements(bs)    eq 1 then bs=replicate(bs,s[0])
     if n_elements(nbins) eq 1 then nbins=replicate(nbins,s[0])
  endif 
  
  if ~array_equal(mn le mx,1b) then $
     message,'Min must be less than or equal to max.'
  
  if n_elements(bs) eq 0 then begin 
     if n_elements(nbins) ne 0 then begin 
        nbins=long(nbins)       ;No fractional bins, please
        bs=float(mx-mn)/nbins   ;a correct formulation
     endif else message,'Must pass either binsize or NBINS'
  endif else nbins=long((mx-mn)/bs+1) 
  
  total_bins=product(nbins,/PRESERVE_TYPE) ;Total number of bins
  h=long((V[s[0]-1,*]-mn[s[0]-1])/bs[s[0]-1])
  ;; The scaled indices, s[n]+N[n-1]*(s[n-1]+N[n-2]*(s[n-2]+...
  for i=s[0]-2,0,-1 do h=nbins[i]*h + long((V[i,*]-mn[i])/bs[i])
  
  out_of_range=[~array_equal(mn le imn,1b),~array_equal(mx ge imx,1b)]
  if ~array_equal(out_of_range,0b) then begin 
     in_range=1
     if out_of_range[0] then $  ;out of range low
        in_range=total(V ge rebin(mn,s,/SAMP),1,/PRESERVE_TYPE) eq s[0]
     if out_of_range[1] then $  ;out of range high
        in_range AND= total(V le rebin(mx,s,/SAMP),1,/PRESERVE_TYPE) eq s[0]
     h=(temporary(h) + 1L)*temporary(in_range) - 1L
  endif 

  ret=make_array(TYPE=3,DIMENSION=nbins,/NOZERO)
  if arg_present(ri) then $
     ret[0]=histogram(h,MIN=0L,MAX=total_bins-1L,REVERSE_INDICES=ri) $
  else $
     ret[0]=histogram(h,MIN=0L,MAX=total_bins-1L)
;  return,ret
end


;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro srcfilter_getind,nx,rev,xind,yind,allind,count=count

undef,allind
count = 0

ind = xind+nx*yind
if rev[ind] ne rev[ind+1] then $
  allind = rev[rev[ind]:rev[ind+1]-1]
count = n_elements(allind)

end


;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


pro srcfilter,list1,list2,outlist

; Filter a source list.  Only keep sources that are in list1.
; LIST1 must be a DAOPHOT aperture photometry .AP file
; 

nlist1 = n_elements(list1)
nlist2 = n_elements(list2)
noutlist = n_elements(outlist)

; Not enough inputs
if nlist1 eq 0 or nlist2 eq 0 or noutlist eq 0 then begin
  print,'Syntax - srcfilter,list1,list2,outlist'
  return
endif

; Check that the files exist
test1 = FILE_TEST(list1)
if (test1 eq 0) then begin
  print,list1,' NOT FOUND'
  return
endif
test2 = FILE_TEST(list2)
if (test2 eq 0) then begin
  print,list2,' NOT FOUND'
  return
endif


; Read in LIST1
; This must be a DAOPHOT .ap aperture photometry file
;----------------------------------------------------
nlines1 = FILE_LINES(list1)
openr,unit1,/get_lun,list1

; Read in the header, 4 lines
line = ''
header1 = strarr(4)
for i=0,3 do begin
  readf,unit1,line
  header1[i] = line
end

nap = floor((nlines1-4)/3)
x1 = fltarr(nap)
y1 = fltarr(nap)
lines1a = strarr(nap)
lines1b = strarr(nap)
line1='' & line2='' & blank=''

; Read the rest of the lines
for i=0,nap-1 do begin
  readf,unit1,line1
  readf,unit1,line2
  readf,unit1,blank
  lines1a[i] = line1
  lines1b[i] = line2
  reads,line1,id,x,y
  x1[i] = x
  y1[i] = y
end
close,unit1
free_lun,unit1


; Read in LIST2
;--------------
; This should be in the standard DAOPHOT COO/ALS format
; Three header lines, and the following lines have ID  X Y ....

nlines2 = FILE_LINES(list2)
openr,unit2,/get_lun,list2
x2 = fltarr(nlines2-3)
y2 = fltarr(nlines2-3)
lines2 = strarr(nlines2)
line=''
for i=0,nlines2-1 do begin
  readf,unit2,line
  lines2[i] = line
  if (i ge 3) then begin
    reads,line,id,x,y
    x2[i-3] = x
    y2[i-3] = y
  endif
end
close,unit2
free_lun,unit2


; Start the outlines
; Use the header lines from LIST1
outlines = header1


; Make an array of Y values from Ymin to Ymax
; then make an array that contains the X2/Y2 indices
; that closest matches those.
; maybe HISTOGRAM gives this.
;ymin = floor(min([y1,y2]))
;ymax = ceil(max([y1,y2]))
;hist2 = histogram(y2,bin=1,min=ymin,max=ymax,locations=yarr,reverse_indices=revind)
; locations are the left hand side of the bin.

; Use a 2D histogram and "reverse" indices to quickly get
; nearby sources
xmin = floor(min([x1,x2]))
xmax = ceil(max([x1,x2]))
ymin = floor(min([y1,y2]))
ymax = ceil(max([y1,y2]))
mnarr = [xmin,ymin]
mxarr = [xmax,ymax]
n2 = n_elements(x2)
V = fltarr(2,n2)
V[0,*] = x2
V[1,*] = y2

SRCFILTER_HIST_ND,V,1,hist2d,MIN=mnarr,MAX=mxarr,REVERSE_INDICES=rev
nx = n_elements(hist2d[*,0])
ny = n_elements(hist2d[0,*])

;         E.g., to find the indices of points
;         which fell in a histogram bin [i,j,k], look up:
;
;             ind=[i+nx*(j+ny*k)]
;             ri[ri[ind]:ri[ind+1]-1]
;
; The subscripts of the original array elements falling the ith bin,
; are given by the expression: R(R[i]:R[i+1]-1).  If R[i]==R[i+1] then
; no elements are present in the ith bin.



; Check that an entry in list2 is closer than MINDIST
nrepeat = 0
mindist = 1.0  ; Minimum matching distance
COMPARE:
bestdist = fltarr(nap)+999999.
; Loop through list1
FOR i=0,nap-1 do begin

  ; The left and bottom edge are included in the bin
  ; if xmin = 1 and ymin=1
  ; then [1.0,1.0] will fall into the first bin
  ; but [2.0,1.0] will fall into the second bin.

  ; Histogram index of this point
  xind = floor(x1[i]-mnarr[0])
  yind = floor(y1[i]-mnarr[1])
  ; Histogram indices of neighbors
  xlo = (xind-1) > 0
  xhi = (xind+1) < (nx-1)
  ylo = (yind-1) > 0
  yhi = (yind+1) < (ny-1)

  ; Get the indices of neighboring points
  undef,allind
  for j=long64(xlo),xhi do begin
    for k=long64(ylo),yhi do begin
      ind=j+nx*k
      ;nel = rev[rev[ind+1]-1] - rev[rev[ind]]
      ;nel = hist2d[j,k]
      if hist2d[j,k] gt 0 then begin
        if n_elements(allind) eq 0 then begin
          allind = rev[rev[ind]:rev[ind+1]-1]
        endif else begin
          allind = [allind, rev[rev[ind]:rev[ind+1]-1] ]
        endelse
      endif
    end
  end 
  nallind = n_elements(allind)

  ; Some neighbors found
  if (nallind gt 0) then begin
    dist = sqrt( (x1[i]-x2[allind])^2.0 + (y1[i]-y2[allind])^2.0 )
    bestdist[i] = min(dist)
  endif
END

ind1 = where(bestdist le mindist,count)
; Some matches
outlines = strarr(nlines1)
outlines[0:3] = header1
if (count gt 0) then begin
  for i=0L,count-1 do begin
    outlines[i*3+4] = lines1a[ind1[i]]
    outlines[i*3+5] = lines1b[ind1[i]]
    ;outlines[i*3+4+2] = ''
  end
endif

; Trim outlines
outlines = outlines[0:count*3+4-1]


; Make sure that there are at least 6 sources
if n_elements(outlines) lt 9 then begin

  ; Lower the minimum matching distance
  if (nrepeat lt 5) then begin
    ;outlines = header1       ; restart the output array
    mindist = mindist*1.5
    nrepeat++
    goto,compare

  ; Just add the middle 10 sources
  endif else begin
    lo = 3 < round(nap*0.5+3) < (nap-1)
    hi = (lo+9) < (nap-1)
    for i=lo,hi do outlines = [outlines,lines1a[i],lines1b[i],'']
  endelse

endif

; Remove last blank line
noutlines = n_elements(outlines)
if outlines[noutlines-1] eq '' then outlines=outlines[0:noutlines-2]

; Write the output file
openw,unit,/get_lun,outlist
noutlines = n_elements(outlines)
for i=0L,noutlines-1 do printf,unit,outlines[i]
close,unit
free_lun,unit

end
