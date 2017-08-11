;+
; NAME:
;	DOUBLES
;
; PURPOSE:
;	Return the subscripts of the non-unique elements in an array,
;       i.e. doubles.
;
;       This procedure automatically sorts the array.  The default is
;       to return only the subscript of one occurance of the non-unique
;       value (i.e. if there are 3 "5"s then only the subscript for one
;       of them is returned; which one is arbitrary).  If all subscripts
;       of all non-unique elements are desired set the /all keyword.
;
;	This function is a simple modification of IDL's uniq function.
;
;       This function can be used to find values that occur in two
;       different arrays.  Just concatenate them.
;
; CATEGORY:
;	Array 
;
; CALLING SEQUENCE:
;	ind = DOUBLES(Array [, all=all, count=count])
;
; INPUTS:
;	Array:	The array to be scanned.  The type and number of dimensions
;		of the array are not important.
;
; OPTIONAL INPUT PARAMETERS:
;      /all     Return subscripts of all non-unique elements of array. 
;               The default is only to return the subscript of one occurance
;               of the non-unique value.
;
; OUTPUTS:
;	An array of indicies into ARRAY is returned containing the subscripts
;       of one occurance of each non-unique element.
;       The expression:
;
;		ARRAY(DOUBLES(ARRAY))
;
;	will be a copy of ARRAY of elements (only one occurance of each) that
;       occur at least twice in ARRRAY (i.e. doubles).
;
;       DBL  This is the actual array of double elements, i.e ARRAY(DOUBLES(ARRAY))
;            This is undefined if count=0.
;
; OPTIONAL OUTPUTS:
;       =COUNT   The count of non-doubles.  This is useful to check for
;                no doubles.
;
; COMMON BLOCKS:
;	None.
;
; MODIFICATION HISTORY:
;       Written by Edward C. Wiebe, 2002-07-23 (Modified IDL's
;       uniq.pro into nonuniq.pro - so this was really written by RSI)
;       Modified by David Nidever   March 2007
;
;-

function DOUBLES, ARRAY, dbl, all=all, count=count

; Check the arguments.
  s = size(ARRAY)
  if (s[0] eq 0) then return, 0		;A scalar

  ; Undefine the DBL array
  undefine,dbl

  ; Sort automatically
  idx = sort(array)

  q = array[idx]

  ; This gets only the UNIQUE doubles
  if not keyword_set(all) then begin
    indices = where(q eq shift(q,-1) and q ne shift(q,1), count)
  endif else begin
    indices = where(q eq shift(q,-1) or q eq shift(q,1), count)
  endelse

  ; No doubles
  if (count eq 0) then return,-1

  ; Doubles
  dblind = idx[indices]

  ; Sort them
  dblind = dblind[sort(dblind)]

  ; Getting the doubles
  dbl = array[dblind]

  return,dblind

end
