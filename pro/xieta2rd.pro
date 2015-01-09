;+
; NAME:
;  xieta2rd
; PURPOSE: (one line)
;  xi,eta coordinates to ra, dec
; DESCRIPTION:
;  
; CATEGORY:
;  Astronomy
; CALLING SEQUENCE:
;  xieta2rd, A, D, xi, eta, ra, dec
; INPUTS:
;  A - ra of the center of a plane (radians)
;  D - dec of the center of the plane (radians)
;  xi - distance east of a point from center of plane (radians)
;  eta - distance north of a point from center of plane (radians)
; OPTIONAL INPUT PARAMETERS:
;  none
; KEYWORD INPUT PARAMETERS:
;  none
; KEYWORD OUTPUT PARAMETERS:
;  none
; OUTPUTS:
;  ra - ra of the point (radians)
;  dec -dec of the point (radians) 
; COMMON BLOCKS:
; SIDE EFFECTS:
; RESTRICTIONS:
; PROCEDURE:
; For a plane centered at A, D, find the ra, dec
; for the given tangent point
; coordinates xi, eta cooresponding to a target at ra, dec
;
; Green, eq. 13.13
; MODIFICATION HISTORY:
;  Written 2005 Jan 27, by Leslie Young, SwRI
;  Moved to $layoung/astronomy 2005 Mar 5
;-

pro xieta2rd, A, D, xi, eta, ra, dec 

    q = cos(D)-eta*sin(D)
    raa = atan(xi/q)
    ra = raa + A
    dec = atan( (sin(D)+eta*cos(D) ) * cos(raa), q)
    
end
