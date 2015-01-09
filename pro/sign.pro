Function        SIGN, number
;+
;NAME:
;       SIGN
;CALLING SEQUENCE:
;       Result = SIGN( num )
;PURPOSE:
;       Returns SIGN (-1, 0 or +1) of NUM.  If NUM is an array, SIGN
;       returns
;       an array of -1, 0 and 1. 
;OUTPUT:
;       Result of function is an array (depending on the original
;       parameter) of
;       type INT.
;EXAMPLE:
;       x = 0.1*Sign(Y)
;HISTORY
;       Written by      J. D. Offenberg, Hughes-STX,    Nov 4, 1991
;-
on_error,2

nelems = n_elements(Number)

SGN = intarr(nelems)
QUE = where(number NE 0)
IF QUE(0) GE 0 then $
        SGN(QUE) = NUMBER(QUE)/abs(NUMBER(QUE))
return, SGN
end
