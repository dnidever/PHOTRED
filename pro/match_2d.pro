;+
; NAME:
;       MATCH_2D
;
; PURPOSE:
;
;       Perform a match between two sets of 2D coordinates, finding
;       the closest coordinate match to the search set, within some
;       search radius.
;
; CALLING SEQUENCE:
;       match=MATCH_2D(x1,y1,x2,y2,search_radius,MATCH_DISTANCE=)
;
; INPUTS:
;
;       x1,y1: The target list to search for matches, of length n1.
;
;       x2,y2: The search list of length n2, to be searched for
;          matches to [x1,y1].
;
;       search_radius: The search radius within which matches will be
;          found.  Only if the closest matching coordinate in [x2,y2]
;          is within the search radius will it be returned.  This is a
;          critical variable in tuning the resources required by
;          MATCH_2D.  See NOTES.
;
; KEYWORD PARAMETERS:
;
;       MATCH_DISTANCE: On output, the distances between the matches
;          and the returned coordinate from the set [x1,y1]
;          (<=search_radius), is returned.
;
; OUTPUTS:
;
;       match: A 1D vector of length n containing the indices of x2
;          and y2 for the closest match to [x1,y1], within the
;          search_radius. If no match was found within the search
;          radius, -1 is returned at that location.
;
; EXAMPLE:
;
;       n=100000
;       x1=randomu(sd,n) & y1=randomu(sd,n)
;       x2=randomu(sd,n) & y2=randomu(sd,n)
;       match=match_2d(x1,y1,x2,y2,.002,MATCH_DISTANCE=md)
;
; SEE ALSO:
;
;       HISTOGRAM, HIST_ND
;
; NOTES:
;
;       This match program uses HIST_ND to pre-bin the 2D search
;       coordinates based on position, within some canonical
;       search_radius.  See:
;
;             http://www.dfanning.com/code_tips/matchlists.html
;
;       for a discussion of its methods.  Of principle importance in
;       the efficiency and behavior of the match is the SEARCH_RADIUS
;       parameter.  This the the maximum radial distance within which
;       a match point must lie to be returned.  Here's a diagram
;       illustrating the problem and method.  For each target point,
;       four separate bins each of width 2*search_radius are searched
;       among the binned search list, depending on its location within
;       the bin.
;
;
;           +----------+----------+
;           |  t1      |          |
;           |          |          |
;           |          |          |
;           |          |          |
;           |          |          |
;           +----------+----------+   -
;           |          |          |   |
;           |          |  o       |   |
;           |          |          |   |  2*search_radius
;           |          |          |   |
;           |          |          |   |
;           +----------+----------+   -
;                         t2
;
;
;        The point `t2' is the closest to the search point `o', but it
;        is not within the search_radius, therefore it is not
;        considered.  Instead, point `t1' is found (and discarded),
;        despite the fact that `t2' is closer.
;
;        Ideally, the search radius should be set to something useful
;        in terms of the match (e.g. positional uncertainty, etc.).
;        However, if the input target coordinates ([x2,y2]) span a
;        large range (e.g. the entire sky), it may be necessary to use
;        a larger search radius to avoid an excessively large number
;        of bins.  Typically there will be an optimal search radius
;        which is fastest.  The tradoff is as follows: the larger the
;        search radius, the smaller the number of bins to search, but
;        the more search points must be considered per target point.
;        The smaller the search radius, the smaller the number of
;        search points per bin, but the greater the number of bins.
;        Something like the median inter-point separation is probably
;        close to optimal.
;
;
; MODIFICATION HISTORY:
;
;       Mon Jul 30 10:56:31 2007, J.D. Smith <jdsmith@as.arizona.edu>
;
;	        Written.
;
;-
;#############################################################################
;
; LICENSE
;
;  Copyright (C) 2007 J.D. Smith
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


function match_2d,x1,y1,x2,y2,search_radius,MATCH_DISTANCE=min_dist
  bs = 2*search_radius          ;this is the smallest binsize allowed
  h = hist_nd([1#x2,1#y2],bs,REVERSE_INDICES=ri)
  bs = bs[0]
  d = size(h,/DIMENSIONS)

  ;; Bin location of X1,Y1 in the X2,Y2 grid
  xoff = x1/bs       &  yoff = y1/bs   
  xbin = floor(xoff) &  ybin = floor(yoff) 
  bin = (xbin + d[0]*ybin)<(d[0]*d[1]-1L) ;The 1D index of the bin it's in

  ;; We must search 4 bins worth for closest match, depending on
  ;; location within bin (i.e. towards any of 4 qudrant directions
  ;;  ul, ur, ll, lr).
  xoff = 1-2*((xoff-xbin) lt 0.5)         ;add bin left or right
  yoff = 1-2*((yoff-ybin) lt 0.5)         ;add bin down or up
  
  n1=n_elements(x1) 
  min_pos = make_array(n1,VALUE=-1L)
  min_dist = fltarr(n1,/NOZERO)
  
  rad2=search_radius^2
  
  for i=0,1 do begin ;; Loop over 4 bins in the correct quadrant direction
     for j=0,1 do begin 
        ;; One of 4 search bins for all the target points
        b = 0L>(bin+i*xoff+j*yoff*d[0])<(d[0]*d[1]-1) 
         
        ;; Dual HISTOGRAM method, loop by count in the search bins
        h2 = histogram(h[b],OMIN=om,REVERSE_INDICES=ri2)
        
        ;; Process all bins with the same repeats (>= 1) at a time
        for k=long(om eq 0),n_elements(h2)-1 do begin 
           if h2[k] eq 0 then continue
           these_bins = ri2[ri2[k]:ri2[k+1]-1] ;bins with k+om search points
            
           if k+om eq 1 then begin ; single point
              these_points = ri[ri[b[these_bins]]]
           endif else begin     ; range over k+om points, (n x k+om)
              targ=[h2[k],k+om]
              these_points = ri[ri[rebin(b[these_bins],targ,/SAMPLE)]+ $
                                rebin(lindgen(1,k+om),targ,/SAMPLE)]
              these_bins = rebin(temporary(these_bins),targ,/SAMPLE)
           endelse
           
           ;; Closest distance squared within this quadrant's bin
           dist = (x2[these_points]-x1[these_bins])^2 + $
                  (y2[these_points]-y1[these_bins])^2 
           
           if k+om gt 1 then begin ;multiple points in bin: find closest
              dist = min(dist,DIMENSION=2,p)
              these_points = these_points[p] ;index of closest point in bin
              these_bins = ri2[ri2[k]:ri2[k+1]-1] ;original bin list
           endif 
            
           ;; See if a minimum is already set there
           set = where(min_pos[these_bins] ge 0, nset, $
                       COMPLEMENT=unset,NCOMPLEMENT=nunset)
           
           if nset gt 0 then begin 
              ;; Only update those where the new point is closer
              closer = where(dist[set] lt min_dist[these_bins[set]], cnt)
              if cnt gt 0 then begin 
                 set = set[closer]
                 min_pos[these_bins[set]] = these_points[set]
                 min_dist[these_bins[set]] = dist[set]
              endif 
           endif 
           
           if nunset gt 0 then begin ;; Nothing set, closest by default
              wh=where(dist[unset] lt rad2,cnt) ;demand it's within radius
              if cnt gt 0 then begin 
                 unset=unset[wh]
                 min_pos[these_bins[unset]] = these_points[unset]
                 min_dist[these_bins[unset]] = dist[unset]
              endif 
           endif 
        endfor 
     endfor 
  endfor 
  if arg_present(min_dist) then min_dist=sqrt(min_dist)
  return,min_pos
end 




































