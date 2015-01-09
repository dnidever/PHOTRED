function readpar,array,keyword,stp=stp,error=error

;+
; 
; READPAR
;
; This allows you to get a parameter value from
; a 2xN array of keyword/value pairs.
; This is similar to getting keyword values from
; headers with SXPAR.PRO.
;
; INPUTS:
;  array    A 2xN array of keyword-value pairs
;  keyword  A keyword string for which to return the value.
;            Case insensitive.
;
; OUTPUTS:
;  value    The value corresponding to the input keyword
;            is output.  If none is found then '0' is returned.
;            If the keyword exists in array but has not value
;            then an empty string '' is returned.
;  =error   The error message if there was one, else undefined
;
; USAGE:
;  IDL>value = readpar(setup,'MOSAIC')
;
; By D. Nidever    Oct. 2007
;-

; Error Handling
;------------------
; Establish error handler. When errors occur, the index of the  
; error is returned in the variable Error_status:  
CATCH, Error_status 

;This statement begins the error handler:  
if (Error_status ne 0) then begin 
   print,'READPAR ERROR: ', !ERROR_STATE.MSG  
   error = !ERROR_STATE.MSG            ; There was an error
   CATCH, /CANCEL 
   return,-1
endif


n = n_elements(array)
if n eq 0 then return,'-1'                    ; must exist
sz = size(array)
if sz[0] ne 2 or sz[1] ne 2 then return,'-1'  ; must be 2xN
if size(array,/type) ne 7 then return,'-1'    ; must be string
if n_elements(keyword) eq 0 then return,'-1'  ; keyword must exist

; Looking for keyword
keyword2 = strlowcase(strtrim(keyword[0],2))
keys = reform(array[0,*])
values = reform(array[1,*])
gd = where(strlowcase(keys) eq keyword2,ngd)

if keyword_set(stp) then stop

if ngd eq 0 then return,'0'
value = strtrim(values[gd[0]],2)             ; returning the first value
return,value

end
