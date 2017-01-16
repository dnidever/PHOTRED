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
   
return, isin

end





