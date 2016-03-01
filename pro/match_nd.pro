;+
; NAME:
;    MATCH_ND
;
; PURPOSE:
;    For each arbitrarily-dimensioned point in one vector, determines the closest
;    point in a second vector. Method is to take the list returned by MATCHALL_ND
;    and narrow it down to the closest match.
;
; CATEGORY:
;    Astro
;
; CALLING SEQUENCE:
;    Result = MATCH_ND(P1, P2, MaxDistance)
;
; INPUTS:
;    P1:      N1 x D array of D-dimensional coordinates.
;
;    P2:      N2 x D array of D-dimensional coordinates.
;
;    MaxDistance:  Maximum D-dimensional distance.
;
; KEYWORD PARAMETERS:
;    MINDIST:  Optional output containing an array of the actual distance
;              to the closest match for each element of P1.
;
;    ONE_TO_ONE:  Enforces one-to-one matching. By default, matching can
;                 be many-to-one, i.e. one entry in P2 can be the
;                 closest match to several entries in P1. If
;                 /ONE_TO_ONE is given, then each multiple entry in P2
;                 is first assigned to its closest point in P1. Then any
;                 entries in P1 that have lost their match are assigned
;                 the next closest point within MaxDistance. This process is
;                 iterated until all points in P1 are matched to a
;                 unique point in P2 or there are no more points
;                 within MaxDistance.
;
; OUTPUTS:
;    The function returns an array with one entry for each P1
;    element containing the index in P2 that is closest, or -1 if
;    there are none within MaxDistance.
;
; EXAMPLE:
;    na = 10
;    nb = 1000
;    a = randomn(seed, na, 3)
;    b = 2. * randomu(seed, nb, 3) - 1
;    matches = match_nd(a, b, 0.3)
;    iplot, /scatter, /iso, b[*,0], b[*,1], b[*,2]
;    iplot, /overplot, /scatter, a[*,0], a[*,1], a[*,2], sym_index=6, $
;      sym_color=[255,0,0]
;    for i=0L, na-1 do if matches[i] ne -1 then $
;      iplot, /overplot, [a[i,0], b[matches[i],0]], [a[i,1], b[matches[i],1]], $
;        [a[i,2], b[matches[i],2]], sym_index=0, color=[255,0,0]
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
;    9 Feb 2011     Forked as MATCH_ND to handle Euclidean case of arbitrary dimension.
;    13 April 2011  Fixed to work around cumulative total 1-element array bug.
;-
function match_nd, p1, p2, maxdistance, mindist=mindist, one_to_one=one2onep

if (size(maxdistance))[0] ne 0 then message, 'MaxDistance must be a scalar.'
p1size = size(p1,/dimen)
p2size = size(p2,/dimen)
if n_elements(p1size) ne 2 then $
  message, 'P1 must be an N1xD dimensional array.'
if n_elements(p2size) ne 2 then $
  message, 'P2 must be an N2xD dimensional array.'
if p1size[1] ne p2size[1] then $
  message, 'P1 and P2 must have the same number of dimensions.'

ndimen = p1size[1]
n1 = p1size[0]
n2 = p2size[0]

; get all matches within maxdistance
matches = matchall_nd(p1, p2, maxdistance, nwithin, distance=distance)

; if there are no matches, just quit
if total(nwithin) eq 0 then return, replicate(-1L, n1)

; sort all by distance
; create unique distance by adding the maximum distance times i1 to each
; entry corresponding to i1:
; first use histogram magic to create a list of which entry in 1 each
; match corresponds to
; note: using [] around cumulative total to deal with 1-element array bug
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

