; for usemap=0, the mapping is just the inverse of array_indices.
; for usemap=1, uses the index in map that corresponds to the above
; (-1 if it doesn't occur)
function wsrad_map, ngrid, a0, a1, a2, map, usemap
  index = a0 + ngrid*(a1 + ngrid*a2)
  if ~usemap then return, index
  result = value_locate(map, index)
  missing = where(map[result] ne index, nmissing)
  if nmissing gt 0 then result[missing]=-1
  return, result
end


;+
; NAME:
;    MATCHALL_SPH
;
; PURPOSE:
;    Determines which of a set of angular coordinates on the sky (or on a
;    sphere) are within a given angular distance from each of a vector of
;    points. New optimized version that uses histograms based on 3D
;    locations on the unit sphere and borrows heavily from JD's MATCH_2D.
;
; CATEGORY:
;    Astro
;
; CALLING SEQUENCE:
;    Result = MATCHALL_SPH(Ra1, Dec1, Ra2, Dec2, Sphrad, Nwithin)
;
; INPUTS:
;    Ra1:     Vector of longitude coordinates, in degrees.
;
;    Dec1:    Vector of latitude coordinates, in degrees.
;
;    Ra2:     Vector of longitude coordinates, in degrees.
;
;    Dec2:    Vector of latitude coordinates, in degrees.
;
;    Sphrad:  Maximum angular distance, in degrees.
;
; OUTPUTS:
;    The function returns the list of indices of Ra2, Dec2 that lie within
;    Sphrad of each point Ra1,Dec1. The format of the returned array is
;    similar to the REVERSE_INDICES array from HISTOGRAM: the indices
;    into Ra2,Dec2 that are close enough to element i of Ra1,Dec1 are
;    contained in Result[Result[i]:Result[i+1]-1] (note, however, that
;    these indices are not guaranteed to be sorted). If there are no matches,
;    then Result[i] eq Result[i+1].
;
; OPTIONAL OUTPUTS:
;    Nwithin: A vector containing the number of matches for each of Ra1,Dec1.
;
; KEYWORD PARAMETERS:
;    DISTANCE:  Optional output containing the distances between each pair.
;               The distances are stored in the same order as the Result
;               array but starting at 0, i.e. if j is match number k to
;               element i then
;                   j = Result[Result[i]+k]
;               and the distance between points i and j is
;                   DISTANCE[Result[i]+k-Result[0]]
;
; EXAMPLE:
;    Note that the routine is similar to finding
;      WHERE(SPHDIST(Ra1,Dec1,Ra2,Dec2,/DEGREE) LE Sphrad, Nwithin)
;    for each element of Ra1 and Dec1, but is much more efficient.
;
;    Shows which random points are within 10 degrees of various coordinates:
;
;    seed=43
;    nrandcoords = 5000l
;    ra_randcoords = 360. * RANDOMU(seed, nrandcoords)
;    dec_randcoords = ASIN( 2*RANDOMU(seed, nrandcoords)-1 ) * !RADEG
;    ra_centers = 60. * FINDGEN(5)
;    dec_centers = [0., 45., 0., -45., 90.]
;    matches = MATCHALL_SPH(ra_centers, dec_centers, ra_randcoords, $
;      dec_randcoords, 10., nmatches)
;    plot, /iso, psym=3, ra_randcoords, dec_randcoords
;    oplot, psym=1, color=fsc_color('blue'), ra_centers, dec_centers
;    oplot, psym=3, color=fsc_color('red'), ra_randcoords[matches[6:*]], $
;      dec_randcoords[matches[6:*]]
;
; MODIFICATION HISTORY:
;    Written by:    Jeremy Bailin
;    10 June 2008   Public release in JBIU as WITHINSPHRAD
;    24 April 2009  Vectorized as WITHINSPHRAD_VEC
;    25 April 2009  Polished to improve memory use
;    9 May 2009     Radical efficiency re-write as WITHINSPHRAD_VEC3D borrowing
;                   heavily from JD Smith's MATCH_2D
;    13 May 2009    Removed * from LHS index in final remapping for speed
;    8 Sept 2010    Renamed MATCHALL_SPH and added DISTANCE keyword
;    9 Aug 2011     Bug fix: incorrect bin size caused occasional matches
;                   to be missed (thanks to J. Donley for reporting)
;-
function matchall_sph, ra1, dec1, ra2, dec2, sphrad, nwithin, $
  distance=distance

if n_elements(ra2) ne n_elements(dec2) then $
  message, 'RA2 and DEC2 must have the same number of elements.'
if n_elements(ra1) ne n_elements(dec1) then $
  message, 'RA1 and DEC1 must have the same number of elements.'
if n_elements(sphrad) ne 1 then $
  message, 'SPHRAD must contain one element.'

n1 = n_elements(ra1)
n2 = n_elements(ra2)

; figure out grid spacing
; simple approximation gives a maximum theta of 2 sqrt(2) / N for
; N cells. It can really be higher, so use a 50% buffer.
; Relevant theta is 2xsphrad.
theta = 2. * sphrad / !radeg
ngrid = (floor(2. * sqrt(2) / (1.5 * theta)) + 1) > 2
gridlen = 2. / (ngrid-1)
; shift edge of box
minbox = -1. - 0.5*gridlen

; create a mapping so we don't need to histogram the entire
; 3d space
; we only need to map the points that contain either an ra1,dec1
; pair or an ra2,dec2 pair
xoff = (cos(ra1/!radeg)*cos(dec1/!radeg)-minbox)/gridlen
yoff = (sin(ra1/!radeg)*cos(dec1/!radeg)-minbox)/gridlen
zoff = (sin(dec1/!radeg)-minbox)/gridlen
xbin = floor(xoff)
ybin = floor(yoff)
zbin = floor(zoff)
x2bin = floor((cos(ra2/!radeg)*cos(dec2/!radeg)-minbox)/gridlen)
y2bin = floor((sin(ra2/!radeg)*cos(dec2/!radeg)-minbox)/gridlen)
z2bin = floor((sin(dec2/!radeg)-minbox)/gridlen)
indices = wsrad_map(ngrid,[xbin,x2bin],[ybin,y2bin],[zbin,z2bin],surfacemap,0)
indexsort = sort(indices)
surfacemap = indices[indexsort[uniq(indices[indexsort])]]
undefine, indexsort
undefine, indices

; histogram points 2
; note the extra 0 out front - used so that when we look for bin "-1"
; (meaning irrelevant), we know there are 0 entries there.
h=[0,histogram(wsrad_map(ngrid,x2bin,y2bin,z2bin,surfacemap,1), omin=hmin, $
  reverse_indices=ri)]

undefine, x2bin
undefine, y2bin
undefine, z2bin

xoff = 1 - 2*((xoff-xbin) lt 0.5)
yoff = 1 - 2*((yoff-ybin) lt 0.5)
zoff = 1 - 2*((zoff-zbin) lt 0.5)

; loop through all neighbouring cells in correct order
for xi=0,1 do begin
  for yi=0,1 do begin
    for zi=0,1 do begin
      b = wsrad_map(ngrid, xbin+xi*xoff, ybin+yi*yoff, zbin+zi*zoff, $
        surfacemap, 1)

      ; dual histogram method, loop by count in search bins (see JD's code)
      h2 = histogram(h[(b-hmin+1) > 0], omin=om, reverse_indices=ri2)

      ; loop through repeat counts
      for k=long(om eq 0), n_elements(h2)-1 do if h2[k] gt 0 then begin
        these_bins = ri2[ri2[k]:ri2[k+1]-1]

        if k+om eq 1 then begin ; single point
          these_points = ri[ri[b[these_bins]-hmin]]
        endif else begin
          targ=[h2[k],k+om]
          these_points = ri[ri[rebin(b[these_bins]-hmin,targ,/sample)]+ $
            rebin(lindgen(1,k+om),targ,/sample)]
          these_bins = rebin(temporary(these_bins),targ,/sample)
        endelse

        ; figure out which ones are really within
        these_dist = sphdist(ra1[these_bins],dec1[these_bins], $
          ra2[these_points],dec2[these_points],/degree)
        within = where(these_dist le sphrad, nwithin)
        if nwithin gt 0 then begin
          ; have there been any pairs yet?
          if n_elements(plausible) eq 0 then begin
            plausible = [[these_bins[within]],[these_points[within]],[these_dist[within]]]
          endif else begin
            ; concatenation is inefficient, but we do it at most 8 x N1 times
            plausible = [plausible,[[these_bins[within]],[these_points[within]],$
              [these_dist[within]]]]
          endelse
        endif

      endif
    endfor
  endfor
endfor

if n_elements(plausible) eq 0 then begin
  nwithin=replicate(0l,n1)
  return, replicate(-1,n1+1)
endif else begin
  ; use histogram to generate a reverse_indices array that contains
  ; the relevant entries, and then map into the appropriate elements
  ; in 2
  nwithin = histogram(plausible[*,0], min=0, max=n1-1, reverse_indices=npri)
  distance = plausible[npri[n1+1:*],2]
  npri[n1+1] = plausible[npri[n1+1:*],1]
  return, npri
endelse

end

