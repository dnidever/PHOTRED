
; Bit Octal Mask   Meaning
; 12  '4000'o  Setuid: Set user ID on execution.
; 11  '2000'o  Setgid: Set group ID on execution.
; 10  '1000'o  Turn on sticky bit. See the UNIX documentation on chmod(2) for details.
; 9  '0400'o  Allow read by owner.
; 8  '0200'o  Allow write by owner.
; 7  '0100'o  Allow execute by owner.
; 6  '0040'o  Allow read by group.
; 5  '0020'o  Allow write by group.
; 4  '0010'o  Allow execute by group.
; 3  '0004'o  Allow read by others.
; 2  '0002'o  Allow write by others.
; 1  '0001'o  Allow execute by others.

pro file_chmod,file,mode,a_execute=a_execute,a_read=a_read,a_write=a_write,$
               g_execute=g_execute,g_read=g_read,g_write=g_write,$
               o_execute=o_execute,o_read=o_read,o_write=o_write,$
               u_execute=u_execute,u_read=u_read,u_write=u_write,$
               no_expand_path=no_expand_path,setgid=setgid,setuid=setuid,$
               sticky_bit=sticky_bit

nfile = n_elements(file)
if nfile eq 0 then begin
  print,'Syntax - file,mode,a_execute=a_execute,a_read=a_read,a_write=a_write,'
  print,'     	            g_execute=g_execute,g_read=g_read,g_write=g_write,'
  print,'                   o_execute=o_execute,o_read=o_read,o_write=o_write,'
  print,'                   u_execute=u_execute,u_read=u_read,u_write=u_write,'
  print,'       	    no_expand_path=no_expand_path,setguid=setguid,'
  print,'                   setuid=setuid,sticky_bit=sticky_bit'
  return
endif

nmode = n_elements(mode)
For i=0,nfile-1 do begin

  info = file_info(file)
  if info[i].exists eq 0 then begin
    print,file[i],' NOT FOUND'
    continue
  endif

  ; Get octal file permissions for this file
  ; -f on macs
  ;spawn,['stat','-c','%a',file[i]],out,errout,/noshell
  spawn,['ls','-l',file[i]],out,errout,/noshell
  perm = first_el(strsplit(out[0],' ',/extract))
  ;-rw-r-xrwx
  ;---s------ means x and SETUID
  ;---S------ means no x and SETUID
  ;------s--- means x and SETGID
  ;------S--- means no x and SETGID
  ;---------t means x and sticky bit
  ;---------T means no x and sticky bit
  orig_mode = 0L
  for j=9,1 do begin
    bit = 10-j
    if (strmid(perm,j,1) ne '-' and strmid(perm,j,1) ne 'S' and strmid(perm,j,1) ne 'T') then orig_mode += 2^(bit-1)
  endfor
  if strlowcase(strmid(perm,9,1)) eq 't' then orig_mode += '1000'o   ; sticky big
  if strlowcase(strmid(perm,6,1)) eq 's' then orig_mode += '2000'o   ; SETGID
  if strlowcase(strmid(perm,3,1)) eq 's' then orig_mode += '4000'o   ; SETUID

  ; Use input mode
  if nmode gt 0 then begin
    if nmode eq nfile then dec_mode=mode[i] else dec_mode=mode
  endif else dec_mode = orig_mode

  ; places: user + group + all
  ; numbers: 1-execute (x), 2-write (w), 4-read (r)

  ; --- Adding bits ---
  ;  all three groups (user, group, other)
  if keyword_set(a_execute) then  dec_mode+='111'o
  if keyword_set(a_write) then    dec_mode+='222'o
  if keyword_set(a_read) then     dec_mode+='444'o
  ;  others
  if keyword_set(o_execute) then  dec_mode+='001'o
  if keyword_set(o_write) then    dec_mode+='002'o
  if keyword_set(o_read) then     dec_mode+='004'o
  ;  group
  if keyword_set(g_execute) then  dec_mode+='010'o
  if keyword_set(g_write) then    dec_mode+='020'o
  if keyword_set(g_read) then     dec_mode+='040'o
  ;  user
  if keyword_set(u_execute) then  dec_mode+='100'o
  if keyword_set(u_write) then    dec_mode+='200'o
  if keyword_set(u_read) then     dec_mode+='400'o
  ; sticky bit, setuid, setgid
  if keyword_set(sticky_bit) then dec_mode+='1000'o
  if keyword_set(setgid) then     dec_mode+='2000'o
  if keyword_set(setuid) then     dec_mode+='4000'o

  ; --- Clearing bits ---
  ;  -keyword is input
  ;  -keyword must be set to ZERO
  ;  -permissions were already granted and can be taken away
  ;  all three groups
  if n_elements(a_execute) gt 0 and not keyword_set(a_execute) and (('111'o and dec_mode) eq '111'o) then dec_mode-='111'o
  if n_elements(a_write) gt 0 and not keyword_set(a_write) and (('222'o and dec_mode) eq '222'o) then dec_mode-='222'o
  if n_elements(a_read) gt 0 and not keyword_set(a_read) and (('444'o and dec_mode) eq '444'o) then dec_mode-='444'o
  ;  others
  if n_elements(o_execute) gt 0 and not keyword_set(o_execute) and (('001'o and dec_mode) eq '001'o) then dec_mode-='001'o
  if n_elements(o_write) gt 0 and not keyword_set(o_write) and (('002'o and dec_mode) eq '002'o) then dec_mode-='002'o
  if n_elements(o_read) gt 0 and not keyword_set(o_read) and (('004'o and dec_mode) eq '004'o) then dec_mode-='004'o
  ;  group
  if n_elements(g_execute) gt 0 and not keyword_set(g_execute) and (('010'o and dec_mode) eq '010'o) then dec_mode-='010'o
  if n_elements(g_write) gt 0 and not keyword_set(g_write) and (('020'o and dec_mode) eq '020'o) then dec_mode-='020'o
  if n_elements(g_read) gt 0 and not keyword_set(g_read) and (('040'o and dec_mode) eq '040'o) then dec_mode-='040'o
  ;  user
  if n_elements(u_execute) gt 0 and not keyword_set(u_execute) and (('100'o and dec_mode) eq '100'o) then dec_mode-='100'o
  if n_elements(u_write) gt 0 and not keyword_set(u_write) and (('200'o and dec_mode) eq '200'o) then dec_mode-='200'o
  if n_elements(u_read) gt 0 and not keyword_set(u_read) and (('400'o and dec_mode) eq '400'o) then dec_mode-='400'o
  ; sticky bit, setuid, setgid
  if n_elements(sticky_bit) gt 0 and not keyword_set(sticky_bit) and (('1000'o and dec_mode) eq '1000'o) then dec_mode-='1000'o
  if n_elements(setgid) gt 0 and not keyword_set(setgid) and (('2000'o and dec_mode) eq '2000'o) then dec_mode-='2000'o
  if n_elements(setuid) gt 0 and not keyword_set(setuid) and (('4000'o and dec_mode) eq '4000'o) then dec_mode-='4000'o

  ; Convert from decimal to octal
  octal_mode = string(format='(O)',dec_mode)  

  ; Perform the action
  spawn,['chmod',octal_mode,file],out,errout,/noshell

Endfor

end
