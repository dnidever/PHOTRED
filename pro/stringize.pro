function stringize,num,ndec=ndec,nocomma=nocomma,unc=unc,$
        sigfig=sigfig,length=length,sci=sci

;+
;
; STRINGIZE
;
; This program helps format numbers for output.  For
; numbers larger than 1000 it puts commas after every 
; third digit.  For numbers less than 1 it converts it
; to decimal notation from scientific notation.
;
; INPUT
;  num      The number to be formatted
;
; OPTIONAL INPUT
;  ndec     The number of desired decimal places
;  length   minimum length desired
;  /unc     formats number for uncertainty.  If num > 1 then ndec=0,
;           if num < 0.001 then it only gives one sig.fig., and if
;           0.001 < num < 1 it gives two sig.figs.
;  sigfig   the number of significant figures.  not completely
;           operational yet.
;  /nocomma don't add any commas
;  /sci     Use scientific notation
;
; PROCEDURES USED
;  strtrim0.pro   strips '0' off the end of a string
;  add.pro        adds two "string" numbers together
;
;-

; printing the syntax
if not keyword_set(num) and n_elements(num) eq 0 then begin
  print,' Syntax
  print,' Result = STRINGIZE( NUMBER [, NDEC=value] [, /NOCOMMA] [, /UNC]'
  print,'                     [, SIGFIG=value] [, LENGTH=value])'
  return,''
endif

; Array of numbers input
If n_elements(num) gt 1 then begin

  n = n_elements(num)
  str_num = string(num)
  for i=0LL,n-1 do str_num(i) = stringize(num(i),ndec=ndec,nocomma=nocomma,unc=unc,$
        sigfig=sigfig,length=length,sci=sci)

Endif Else begin

; remembing the original inputs
old_num = num
if keyword_set(ndec) then old_ndec = ndec
if keyword_set(nocomma) then old_nocomma = nocomma
if keyword_set(unc) then old_unc =  unc
if keyword_set(sigfig) then old_sigfig = sigfig
if keyword_set(length) then old_length = length

; initial parameter setups
if n_elements(ndec) eq 0 and n_elements(length) gt 0 and num eq 0.0 then ndec = length-2 
if n_elements(ndec) eq 0 and n_elements(length) eq 0 and num eq 0.0 then ndec = 2
if n_elements(ndec) eq 0 and abs(num) lt 1.d then ndec = 12	;default of 12 decimal places
if n_elements(ndec) eq 0 and abs(num) ge 1.d then ndec = 0	;default of 0 decimal places
if n_elements(sigfig) eq 0 then sigfig = 0
;if ndec eq 0 and sigfig ne 0 then ndec = sigfig + 1
;if not keyword_set(ndec) then ndec=0.     ;no decimal places
if n_elements(ndec) gt 0 then if ndec lt 0. then ndec=0.
if keyword_set(sigfig) then if sigfig lt 0. then sigfig=0.
if keyword_set(unc) then begin
  if num ge 1. then ndec=0
 ; if num lt 0.01d then ndec=4
 ; if num gt 0.01d and num lt 1. then sigfig=2     ;two sig.figs.
  if num lt 1. then sigfig = 1
endif
if (abs(num) gt 1e8 or abs(num) lt 1e-8) and (num ne 0.) then sci=1   ; use sci.not. for large/small numbers
if keyword_set(sci) and (abs(num) lt 1000. and abs(num) ge 1.) then sci=0  ; don't use sci.not. for normal #s

num = double(abs(num))
;strnum = strtrim(string(num),2)
;num = double(strnum+'d')                 ;make it double

str_num = ''

; Using Scientific notation
If keyword_set(sci) Then Begin

  if not keyword_set(ndec) then ndec=1
  blen = ndec+6
  str_num = strupcase(string(num,format='(E'+strtrim(blen,2)+'.'+strtrim(ndec,2)+')'))
  str_num = strtrim(str_num,2)
  if signs(old_num) eq -1 then str_num='-'+str_num

; NOT Using Scientific Notation
Endif Else Begin

; NUM ge 1
IF abs(num) ge 1.d THEN BEGIN

  pow = alog10(abs(num))
  pow = floor(pow)
  ;strint = string(num,format='(I40)')
  ;strint = strtrim(strint,2)
  ;pow = strlen(strint)
  blen = long( 40 > pow+2 )
  strnum = string(num,format='(F'+strtrim(blen,2)+'.'+strtrim(long(ndec),2)+')')
  strnum = strtrim(strnum,2)
  ;strnum = strtrim0(strnum,/back)
  len = strlen(strnum)
  dot = strpos(strnum,'.')
  intg = strmid(strnum,0,dot)             ;integer digits
  ;dec = strmid(strnum,dot+1,len-dot-1)    ;decimal digits
  dec = strmid(strnum,dot+1,ndec)    ;decimal digits

  if keyword_set(sigfig) then begin
    npow = strlen(intg)
    ndec = sigfig-npow

    if ndec lt 0. then num = num*10.^ndec
    if ndec lt 0. then ndec=0.

    ;redoing this part
    blen = long( 40 > pow+2 )
    strnum = string(num,format='(F'+strtrim(blen,2)+'.'+strtrim(long(ndec),2)+')')
    strnum = strtrim(strnum,2)
    len = strlen(strnum)
    dot = strpos(strnum,'.')
    intg = strmid(strnum,0,dot)             ;integer digits
    dec = strmid(strnum,dot+1,ndec)    ;decimal digits

    if sigfig-npow lt 0 then begin
      intg = intg + strmult('0',abs(sigfig-npow))
      dec=''
    endif

  endif

  ;integer part
  ncom = floor(pow/3.)

  ; adding commas every three places
  if not keyword_set(nocomma) then begin
    old_strnum = strnum
    ;putting in the commas
    for i=0.,ncom-1. do begin
      len = strlen(intg)
      lo = strpos(intg,',')
      if lo eq -1 then lo=len                  ; for the first time
     ; first = strmid(strnum,0,len-(i+1)*3)
     ; last = strmid(strnum,len-(i+1)*3,(i+1)*3)
      first = strmid(intg,0,lo-3)
      last = strmid(intg,lo-3,len-lo+3)
      intg = first + ',' + last
      ;stop
    endfor
  endif

  if ndec eq 0. then str_num = intg
  if ndec gt 0. then str_num = intg + '.' + dec

  if old_num lt 0. then str_num = '-' + str_num  
  if old_num lt 0. then strnum = '-' + strnum  

  ;checking for error
  diff = abs(double(strnum)-double(old_num))
  perdiff = diff/double(num)
  if ndec gt 0. then begin
    ;off = 1/(10.^double(ndec))
    off = 1d-10 > 1/(10.^double(ndec))
    if diff gt off then begin
      print,'ERROR IN STRINGIZE.PRO'
      stop
    endif
  endif

ENDIF  ; num ge 1

; Num lt 1
IF abs(num) lt 1. AND num ne 0.0 THEN BEGIN         ;decimals

  pow = alog10(abs(num))
  pow = floor(pow)

  if keyword_set(sigfig) then ndec = sigfig-pow-1

  str_num = string(num,format='(F40.'+strtrim(long(ndec),2)+')')
  str_num = strtrim(str_num,2)

  if ndec eq 0 then str_num = strtrim(long(str_num),2)   ; rounding

  if old_num lt 0. then str_num = '-' + str_num           ;negative sign

  ;checking for error
  diff = abs(double(str_num)-double(old_num))
  perdiff = diff/double(num)
  if ndec gt 0. then begin
    ;off = 1/(10.^double(ndec))
    off = 1d-10 > 1/(10.^double(ndec))
    if diff gt off then begin
      print,'ERROR IN STRINGIZE.PRO'
      stop
    endif
  endif

  ;trim off any extra 0's off the back.
  if keyword_set(unc) then begin
    str_num = strtrim0(str_num,/back)
  endif

ENDIF   ; num lt 1

;IF IT'S EXACTLY ZERO
if num eq 0.0 then begin
  ;if keyword_set(length) then ndec=length-2 > 0
  zeros = strmult('0',ndec)
  str_num = '0.'+zeros
  if ndec eq 0 then str_num = '0'
endif


; Dealing with the length
if keyword_set(length) then begin
  len = strlen(str_num)
  if len lt length then begin
    front=strmult(' ',length-len)   ;adding it to the front
    str_num = front+str_num
    ;str_num = str_num + strmult('0',length-len)
  endif

endif                  ;length

Endelse; not using sci.not.

; putting back the original inputs
num = old_num
if keyword_set(old_ndec) then ndec = old_ndec
if keyword_set(old_nocomma) then nocomma = old_nocomma
if keyword_set(old_unc) then unc =  old_unc
if keyword_set(old_sigfig) then sigfig = old_sigfig
if keyword_set(old_length) then length = old_length

Endelse ; array input

;stop

return,str_num

end
