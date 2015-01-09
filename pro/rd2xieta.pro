;+
; NAME:
;  xieta2rd
; PURPOSE: (one line)
;  ra,dec coordinates to xi, eta
; DESCRIPTION:
;  
; CATEGORY:
;  Astronomy
; CALLING SEQUENCE:
;  rd2xieta, A, D, ra, dec, xi, eta
; INPUTS:
;  A - ra of the center of a plane (radians)
;  D - dec of the center of the plane (radians)
;  ra - ra of the point (radians)
;  dec -dec of the point (radians) 
; OPTIONAL INPUT PARAMETERS:
;  none
; KEYWORD INPUT PARAMETERS:
;  none
; KEYWORD OUTPUT PARAMETERS:
;  none
; OUTPUTS:
;  xi - distance east of a point from center of plane (radians)
;  eta - distance north of a point from center of plane (radians)
; COMMON BLOCKS:
; SIDE EFFECTS:
; RESTRICTIONS:
; PROCEDURE:
; For a plane centered at A, D, find the tangent point
; coordinates xi, eta cooresponding to a target at ra, dec
;
; Smart, Text-book on spherical astronomy, 5th ed, section 161
;
; MODIFICATION HISTORY:
;  Written 2005 Jan 27, by Leslie Young, SwRI
;  Moved to $layoung/astronomy 2005 Mar 5
;-

pro rd2xieta, A, D, ra, dec, xi, eta

; cot(q) = cot(dec) cos(ra-A)     (17)
; eta = tan(q-D)                  (18)
; xi = cos(q) tan(ra-A)/cos(q-D)  (20)

q = atan(sin(dec), cos(dec) * cos(ra-A) )
eta = tan(q-D)
xi = cos(q)*tan(ra-A)/cos(q-D) 
  
; print, [xi, eta], [(ra-A)*cos(D), dec-D]

end
