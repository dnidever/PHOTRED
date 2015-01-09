;+
; NAME:
;       SCALE_VECTOR
;
; PURPOSE:
;
;       This is a utility routine to scale the points of a vector
;       (or an array) into a given data range. The minimum value of
;       the vector (or array) is set equal to the minimum data range. And
;       the maximum value of the vector (or array) is set equal to the
;       maximum data range.
;
; AUTHOR:
;
;       FANNING SOFTWARE CONSULTING
;       David Fanning, Ph.D.
;       1645 Sheely Drive
;       Fort Collins, CO 80526 USA
;       Phone: 970-221-0438
;       E-mail: davidf@dfanning.com
;       Coyote's Guide to IDL Programming: http://www.dfanning.com
;
; CATEGORY:
;
;       Utilities
;
; CALLING SEQUENCE:
;
;       scaledVector = SCALE_VECTOR(vector, minRange, maxRange)
;
; INPUT POSITIONAL PARAMETERS:
;
;       vector:   The vector (or array) to be scaled.
;       minRange: The minimum value of the scaled vector. Set to 0 by default.
;       maxRange: The maximum value of the scaled vector. Set to 1 by default.
;
; INPUT KEYWORD PARAMETERS:
;
;       MAXVALUE: Any value in the input vector greater than this value is
;                 set to this value before scaling.
;
;       MINVALUE: Any value in the input vector less than this value is
;                 set to this value before scaling.
;
;       NAN:      Set this keyword to enable not-a-number checking. NANs
;                 in vector will be ignored.
;
; RETURN VALUE:
;
;       scaledVector: The vector (or array) values scaled into the data range.
;           This is always at least a FLOAT value.
;
; COMMON BLOCKS:
;       None.
;
; EXAMPLE:
;
;       x = [3, 5, 0, 10]
;       xscaled = SCALE_VECTOR(x, -50, 50)
;       Print, xscaled
;          -20.0000     0.000000     -50.0000      50.0000
;
;
; MODIFICATION HISTORY:
;
;       Written by:  David W. Fanning, 12 Dec 1998.
;       Added MAXVALUE and MINVALUE keywords. 5 Dec 1999. DWF.
;       Added NAN keyword. 18 Sept 2000. DWF.
;       Removed check that made minRange less than maxRange to allow ranges to be
;          reversed on axes, etc. 28 Dec 2003. DWF.
;-

FUNCTION Scale_Vector, vector, minRange, maxRange, $
   MAXVALUE=vectorMax, MINVALUE=vectorMin, NAN=nan

On_Error, 1

   ; Check positional parameters.

CASE N_Params() OF
   0: Message, 'Incorrect number of arguments.'
   1: BEGIN
      minRange = 0.0
      maxRange = 1.0
      ENDCASE
   2: BEGIN
      maxRange = 1.0 > (minRange + 0.0001)
      ENDCASE
   ELSE:
ENDCASE

   ; Make sure we are working with floating point numbers.

minRange = Float( minRange )
maxRange = Float( maxRange )

   ; Make sure we have a valid range.
;IF maxRange EQ minRange THEN Message, 'Range max and min are coincidental'
IF maxRange EQ minRange THEN begin
  print, 'Range max and min are coincidental'
  ;MESSAGE, 'Range max and min are coincidental'
  return,vector*0+minRange
ENDIF

   ; Check keyword parameters.

IF N_Elements(vectorMin) EQ 0 THEN vectorMin = Float( Min(vector, NAN=Keyword_Set(nan)) ) $
   ELSE vectorMin = Float( vectorMin )
IF N_Elements(vectorMax) EQ 0 THEN vectorMax = Float( Max(vector, NAN=Keyword_Set(nan)) ) $
   ELSE vectorMax = Float( vectorMax )

   ; Trim vector before scaling.

index = Where(Finite(vector) EQ 1, count)
IF count NE 0 THEN BEGIN
   trimVector = vector
   trimVector[index]  =  vectorMin >  (vector[index]) < vectorMax
ENDIF ELSE trimVector = vectorMin > vector < vectorMax

   ; Calculate the scaling factors.

scaleFactor = [((minRange * vectorMax)-(maxRange * vectorMin)) / $
    (vectorMax - vectorMin), (maxRange - minRange) / (vectorMax - vectorMin)]

   ; Return the scaled vector.

RETURN, trimVector * scaleFactor[1] + scaleFactor[0]

END ;-------------------------------------------------------------------------

