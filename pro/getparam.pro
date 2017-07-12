;+
;
; GETPARAM
;
; Function to retrieve a parameter from the PHOTRED setup file.
; A default value can be set.
;
; INPUTS:
;  var        The variable if it was previously set or passed.
;               If this variable is defined and given, then this
;               is returned and the setup file is not queried.
;  name       Name of the variable in the setup file.
;  setup      PHOTRED setup array.
;  default    Default value to use if the parameter is not found
;               in the setup file.
;  logfile    Log file to write message to.
;  /bool      Parameter is a boolean, so for it to 0 or 1.
;
; OUTPUTS:
;  out        The parameter value.  If "var" is input, then out=var.
;               If the parameter is found in the setup file then this
;               value is returned (potentially modified by /bool).
;               If "var" is not input and the parameter is not found
;               in the setup file then the default value is returned.
;
; USAGE:
;  IDL>hyperthread = getparam(hyperthread,'hyperthread',setup,'0',logfile,/bool)
;
; By Antonio Dorta  July 2017
;-

function getparam, var, name, setup, default, logfile, bool=bool

  ; the variable was already given, use it
  if n_elements(var) ne 0 then begin 
    out = var

  ; Get the parameter from the setup file
  endif else begin
    data = READPAR(setup,name)
    ; Parameter not set, use the default value
    if data eq '0' or data eq '-1' or data eq '' then begin 
      out = default
    endif else begin
      out = data
    endelse
  endelse

  if keyword_set(bool) then if out eq '0' then out = 0 else out = 1
    
  if n_elements(logfile) ne 0 then printlog,logfile,"PARAM "+name+": "+strcompress(out)

  return,out

end
