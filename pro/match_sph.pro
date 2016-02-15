;+
; NAME:
;    MATCH_SPH
;
; PURPOSE:
;    For each angular point in one vector, determines the closest angular match
;    from another vector. Method is to take the list returned by
;    MATCHALL_SPH and narrow it down to the closest match.
;
; CATEGORY:
;    Astro
;
; CALLING SEQUENCE:
;    Result = MATCH_SPH(Ra1, Dec1, Ra2, Dec2, Sphrad)
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
;    Sphrad:  Maximum angular distance, in degrees. Matches outside
;             of this radius are ignored.
;
; KEYWORD PARAMETERS:
;    MINDIST:  Optional output containing an array of the actual distance
;               to the closest match for each element of Ra1, Dec1.
;
;    ONE_TO_ONE:  Enforces one-to-one matching. By default, matching can
;                 be many-to-one, i.e. one entry in Ra2,Dec2 can be the
;                 closest match to several entries in Ra1,Dec1. If
;                 /ONE_TO_ONE is given, then each multiple entry in Ra2,Dec2
;                 is first assigned to its closest point in Ra1,Dec1. Then any
;                 entries in Ra1,Dec1 that have lost their match are assigned
;                 the next closest point within Sphrad. This process is
;                 iterated until all points in Ra1,Dec1 are matched to a
;                 unique point in Ra2,Dec2 or there are no more points
;                 within Sphrad
;
; OUTPUTS:
;    The function returns an array with one entry for each Ra1, Dec1
;    element containing the index in Ra2, Dec2 that is closest, or -1 if
;    there are none within Sphrad.
;
; EXAMPLE:
;
;    n1 = 25
;    n2 = 10
;    seed = 43L
;    ra1 = randomn(seed, n1)
;    dec1 = randomn(seed, n1)
;    ra2 = randomn(seed, n2)
;    dec2 = randomn(seed, n2)
;    result1 = match_sph(ra1, dec1, ra2, dec2, 1.)
;    result2 = match_sph(ra1, dec1, ra2, dec2, 1., /one_to_one)
;    !p.multi=[0,2,1]
;    plot, psym=1, ra1, dec1, xrange=[-3,3], yrange=[-3,3], title='Default'
;    oplot, psym=4, ra2, dec2
;    for i=0l,n1-1 do if result1[i] ne -1 then oplot, [ra1[i],ra2[result1[i]]], $
;      [dec1[i],dec2[result1[i]]]
;    plot, psym=1, ra1, dec1, xrange=[-3,3], yrange=[-3,3], title='/ONE_TO_ONE'
;    oplot, psym=4, ra2, dec2
;    for i=0l,n1-1 do if result2[i] ne -1 then oplot, [ra1[i],ra2[result2[i]]], $
;      [dec1[i],dec2[result2[i]]]
;    
;
; MODIFICATION HISTORY:
;    Written by:    Jeremy Bailin
;    10 June 2008   Public release in JBIU as WITHINSPHRAD
;    24 April 2009  Vectorized as WITHINSPHRAD_VEC
;    25 April 2009  Polished to improve memory use
;    9 May 2009     Radical efficiency re-write as WITHINSPHRAD_VEC3D borrowing
;                   heavily from JD Smith's MATCH_2D
;    13 May 2009    Removed * from LHS index in final remapping for speed
;    23 August 2010 Modified to only return closest match as WITHINSPHRAD_CLOSEST
;                   (importing some stuff back from MATCH_2D)
;    8 Sept 2010    Renamed MATCH_SPH. Added /ONE_TO_ONE option, and modified
;                   to explicitly call MATCHALL_SPH and then cull.
;    13 April 2011  Fixed to work around cumulative total 1-element array bug.
;-
function match_sph, ra1, dec1, ra2, dec2, sphrad, mindist=mindist, $
  one_to_one=one2onep

if n_elements(ra2) ne n_elements(dec2) then $
  message, 'RA2 and DEC2 must have the same number of elements.'
if n_elements(ra1) ne n_elements(dec1) then $
  message, 'RA1 and DEC1 must have the same number of elements.'
if n_elements(sphrad) ne 1 then $
  message, 'SPHRAD must contain one element.'

n1 = n_elements(ra1)
n2 = n_elements(ra2)

; get all matches within sphrad
matches = matchall_sph(ra1, dec1, ra2, dec2, sphrad, nwithin, distance=distance)

; if there are no matches, just quit
if total(nwithin) eq 0 then return, replicate(-1L, n1)

; sort all by distance
; create unique distance by adding the maximum distance times i1 to each
; entry corresponding to i1:
; first use histogram magic to create a list of which entry in 1 each
; match corresponds to
; note: using [] around cumulative total to deal with bug for 1-element arrays
hchunk = histogram([total(nwithin,/cumul,/int)]-1, min=0, reverse_indices=chunkri)
index1 = chunkri[0:n_elements(hchunk)-1]-chunkri[0] + (nwithin[0] eq 0)
; now create unique distance by adding max(distance) * index1:
uniqdist = distance + 1.1*max(distance)*index1
; and sort
sortdist = sort(uniqdist)
distance = distance[sortdist]
matches[n1+1] = (matches[n1+1:*])[sortdist]

matchloc = matches[0:n1-1]   ; starting point for each entry in list 1
nomatchp = nwithin eq 0      ; 1 for entries with no match in sphrad

; iterate checking for multiple matches
nincrements = lonarr(n1)
repeat begin
  nomatch = where(nomatchp, nnomatch)
  minpos = matches[matchloc]
  mindist = distance[matchloc - matches[0]]
  if nnomatch gt 0 then begin
    minpos[nomatch] = -1
    mindist[nomatch] = 0.
  endif
  ; if many-to-one matches allowed, no need to do anything else
  if ~keyword_set(one2onep) then break  

  ; histogram of entries in list 2
  h = histogram(minpos, min=0, max=n2-1, reverse_indices=ri)
  ; we're done if there are no multiple matches
  if max(h) le 1 then break

  ; find entries with multiple matches and increment the non-closest
  ; matches
; slower loop version:
;  multimatch = where(h gt 1, nmultimatch)
;  for mi=0l,nmultimatch-1 do begin
;    ; sort the distances associated with this entry
;    this2 = ri[ri[multimatch[mi]]:ri[multimatch[mi]+1]-1]
;    sortdist = sort(distance[matchloc[this2] - matches[0]])
;    notclosest = this2[sortdist[1:*]]
;    matchloc[notclosest]++
;    nincrements[notclosest]++
;  endfor
; faster histogram-of-histogram version:
  h2 = histogram(h, min=2, omax=hmax, reverse_indices=ri2)
  nh2 = n_elements(h2)
  for repci=0l,nh2-1 do if h2[repci] gt 0 then begin
    targ = [h2[repci], repci+2]
    vec_inds = ri2[ri2[repci]:ri2[repci+1]-1]  ; indices into h
    these1 = ri[rebin(ri[vec_inds], targ, /sample) + $
      rebin(transpose(lindgen(repci+2)), targ, /sample)]
    these2 = rebin(vec_inds, targ, /sample)
    ; now sort the distances for each entry in these2
    sortdist = sort_nd(distance[matchloc[these1]-matches[0]], 2)
    notclosest = these1[sortdist[*,1:*]]
    matchloc[notclosest]++
    nincrements[notclosest]++    
  endif

  ; if we've gone past the end of the points within sphrad, mark it as a nomatch
  nomatchp or= nincrements ge nwithin
endrep until 0 ne 0

return, minpos

end

