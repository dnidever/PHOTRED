;+
;
; DOPOLYGONSOVERLAP
;
; This function checks to see if two polygons overlap by using their vertices.
;
; INPUTS:
;  xpolygon1  The array of X-values of the vertices of the first polygon.
;               The first vertex does NOT need to be repeated.
;  ypolygon1  The array of Y-values of the vertices of the first polygon.
;  xpolygon2  The array of X-values of the vertices of the second polygon.
;  ypolgyon2  The array of Y-values of the vertices of the second polygon.
;
; OUTPUTS:
;  isin    If the polygons overlap then isin=1 else isin=0
;
; USAGE:
;  IDL>isin = dopolygonsoverlap(x1,y1,x2,y2)
;
; By D. Nidever Feb 2016
;-


function isleft, x1, y1, x2, y2, x3, y3
; isLeft(): test if a point is Left|On|Right of an infinite 2D line.
;   From http://geomalgorithms.com/a01-_area.html
; Input:  three points P1, P2, and P3
; Return: >0 for P3 left of the line through P1 to P2
; =0 for P3 on the line
; <0 for P3 right of the line
return, ( (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1) )

end

;==================

function rangeoverlap,a,b
; Does the range (start1, end1) overlap with (start2, end2)
return, max(a) ge min(b) and min(a) le max(b)
end
    
;==================

function doLineSegmentsIntersect,x1, y1, x2, y2
; Do two line segments intersect.

;  Check vertical lines

; Vertical lines, but NOT same X-values
if x1[0] eq x1[1] and x2[0] eq x2[1] and x1[0] ne x2[0] then return, 0  ; no overlap

;  Vertical lines with same X values
if x1[0] eq x1[1] and x2[0] eq x2[1] and x1[0] eq x2[0] then begin
  ; Check intersection of Y ranges
  I1 = [min(y1), max(y1)]
  I2 = [min(y2), max(y2)]

  ; And we could say that Xa is included into :
  Ia = [max( min(y1), min(y2) ),$
        min( max(y1), max(y2) )]

  ;  Now, we need to check that this interval Ia exists :
  if rangeoverlap(y1,y2) eq 0 then begin
    return, 0   ; There is no mutual abcisses        
  endif else begin
    return, 1   ; There is overlap
  endelse
endif  ; vertical lines
 
; The equation of a line is:
;
; f(x) = A*x + b = y
; For a segment, it is exactly the same, except that x is included on an interval I.
; 
; If you have two segments, defined as follow:
;
; Segment1 = {(X1, Y1), (X2, Y2)}
; Segment2 = {(X3, Y3), (X4, Y4)}
; The abcisse Xa of the potential point of intersection (Xa,Ya) must be contained in both interval I1 and I2, defined as follow:
I1 = [min(x1), max(x1)]
I2 = [min(x2), max(x2)]

; And we could say that Xa is included into :
Ia = [max( [min(x1), min(x2)] ),$
      min( [max(x1), max(x2)] )]

; Now, we need to check that this interval Ia exists :
if rangeoverlap(x1,x2) eq 0 then return, 0   ; There is no mutual abcisses

; Check that the Y-ranges overlap as well
if rangeoverlap(y1,y2) eq 0 then return, 0   ; There is no mutual y-value overlap

; So, we have two line formula, and a mutual interval. Your line formulas are:
; f1(x) = m1*x + b1 = y
; f2(x) = m2*x + b2 = y

; As we got two points by segment, we are able to determine A1, A2, b1 and b2:
dx1 = x1[1]-x1[0]
if dx1 eq 0 then begin
  m1 = !values.f_nan
  b1 = 0
endif else begin
  m1 = (y1[1]-y1[0])/(x1[1]-x1[0])
  b1 = y1[0]-m1*x1[0]
endelse
dx2 = x2[1]-x2[0]
if dx2 eq 0 then begin
  m2 = !values.f_nan
  b2 = 0
endif else begin
  m2 = (y2[1]-y2[0])/(x2[1]-x2[0])
  b2 = y2[0]-m2*x2[0]
endelse
      
; If the segments are parallel, then m1 == m2:
if (m1 eq m2) and (b1 ne b2) then return, 0  ; Parallel segments

; If the segments are parallel and on top of each other, the m1==m2 and b1==b2
; we've already required that the x-ranges (abcissas) overlap
if (m1 eq m2) and (b1 eq b2) then return, 1   ; parallel segments on top of each other
    
; A point (Xa,Ya) standing on both lines must satisfy both formulas f1 and f2:
; Ya = m1 * Xa + b1
; Ya = m2 * Xa + b2
; A1 * Xa + b1 = m2 * Xa + b2

; Line segment 1 is vertical line 
if x1[0] eq x1[1] then begin
  Xa = x1[0]
  Ya = m2*Xa+b2
  if rangeoverlap(x1,[Xa]) and rangeoverlap(x2,[Xa]) and $
     rangeoverlap(y1,[Ya]) and rangeoverlap(y2,[Ya]) then begin
    return, 1
  endif else begin
    return, 0
  endelse
  
; Line semgent 2 is vertical line
endif else if x2[0] eq x2[1] then begin
  Xa = x2[0]
  Ya = m1*Xa+b1
  if rangeoverlap(x1,[Xa]) and rangeoverlap(x2,[Xa]) and $
     rangeoverlap(y1,[Ya]) and rangeoverlap(y2,[Ya]) then begin   
    return, 1
  endif else begin
    return, 0
  endelse
 
; Neither are vertical lines
endif else begin
  Xa = (b2 - b1) / (m1 - m2)
endelse
        
; The last thing to do is check that Xa is included into Ia:
if ( (Xa lt max( [min(x1), min(x2)] )) or $
     (Xa gt min( [max(x1), max(x2)] )) ) then begin
  return, 0  ; intersection is out of bound
endif else begin
  return, 1
endelse

end
    
;==================

function ispointinpolygon, xpolygon, ypolygon, xpt, ypt
; Returns boolean if a point is inside a polygon of vertices.
;    
; How to tell if a point is inside a polygon:
; Determine the change in angle made by the point and the vertices
; of the polygon.  Add up the delta(angle)'s from the first (include
; the first point again at the end).  If the point is inside the
; polygon, then the total angle will be +/-360 deg.  If the point is
; outside, then the total angle will be 0 deg.  Points on the edge will
; outside.
; This is called the Winding Algorithm
; http://geomalgorithms.com/a03-_inclusion.html

n = n_elements(xpolygon)
; Array for the angles
angle = fltarr(n)

; add first vertex to the end
xpolygon1 = [ xpolygon, xpolygon[0] ]
ypolygon1 = [ ypolygon, ypolygon[0] ]

wn = 0   ; winding number counter

; Loop through the edges of the polygon
for i=0,n-1 do begin
   ; if edge crosses upward (includes its starting endpoint, and excludes its final endpoint)
   if ypolygon1[i] le ypt and ypolygon1[i+1] gt ypt then begin
       ; if (P is  strictly left of E[i])    // Rule #4
       if isleft(xpolygon1[i], ypolygon1[i], xpolygon1[i+1], ypolygon1[i+1], xpt, ypt) gt 0 then $ 
            wn += 1   ; a valid up intersect right of P.x
   endif
       
   ; if edge crosses downward (excludes its starting endpoint, and includes its final endpoint)
   if ypolygon1[i] gt ypt and ypolygon1[i+1] le ypt then begin
       ; if (P is  strictly right of E[i])    // Rule #4
       if isleft(xpolygon1[i], ypolygon1[i], xpolygon1[i+1], ypolygon1[i+1], xpt, ypt) lt 0 then $
            wn -= 1   ; a valid up intersect right of P.x
   endif
            
endfor
            
; wn = 0 only when P is outside the polygon
if wn eq 0 then return,0 else return,1

end

;==================

function dopolygonsoverlap, xpolygon1, ypolygon1, xpolygon2, ypolygon2
; Returns True if two polygons are overlapping.
;
; How to determine if two polygons overlap.
; If a vertex of one of the polygons is inside the other polygon
; then they overlap.
    
n2 = n_elements(xPolygon2)
isin = 0

; Loop through all vertices of second polygon
for i=0,n2-1 do begin
   ; perform iterative boolean OR
   ; if any point is inside the polygon then they overlap   
   isin = isin or ispointinpolygon(xpolygon1, ypolygon1, xpolygon2[i], ypolygon2[i])
endfor

; Must do the reverse as well, they aren't the same
n1 = n_elements(xPolygon1)
for i=0,n1-1 do begin
   isin = isin or ispointinpolygon(xpolygon2, ypolygon2, xpolygon1[i], ypolygon1[i])
endfor

; Two polygons can overlap even if there are no vertices inside each other.
; Need to check if the line segments overlap
if isin eq 0 then begin
  intersect = 0
  ; Add first vertex to the end
  xp1 = [ xPolygon1, xPolygon1[0] ]
  yp1 = [ yPolygon1, yPolygon1[0] ]
  xp2 = [ xPolygon2, xPolygon2[0] ]
  yp2 = [ yPolygon2, yPolygon2[0] ]
  for i=0,3 do begin
    for j=0,3 do begin
       intersect = intersect or dolinesegmentsintersect(xp1[i:i+1],yp1[i:i+1],xp2[j:j+1],yp2[j:j+1])
       if intersect eq 1 then return, 1
       isin = isin or intersect
    endfor
  endfor
endif
       
return, isin

end





