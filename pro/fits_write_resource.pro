;+
;
; FITS_WRITE_RESOURCE
;
; Write data and a header to a FITS file that has a resource file.
;
; INPUTS:
;  file    The filename to write to.
;  data    The data to be written to the FITS file.  Set this 
;             to 0 if you only want to update the header.
;  head    The image header string array.
;
; OUTPUTS:
;  =error  The error message if one occurred.
;
; USAGE:
;  IDL>fits_write_resource,file,data,head
;
; By D. Nidever  Feb 2019
;-

pro fits_write_resource,file,data,head,error=error

undefine,error

;; Not enough inputs
if n_elements(file) eq 0 or (n_elements(data) eq 0 and n_elements(head) eq 0) then begin
  error = 'Not enough inputs'
  print,'Syntax - fits_write_resource,file,data,head,error=error'
  return
endif

;; Resource file name
dir = file_dirname(file)
rfile = dir+'/.'+file

;; We have a resource file
if file_test(rfile) eq 1 then begin
  info = file_info(file)

  ;; FITS file exists, update data and header
  if info.exists eq 1 and info.size gt 1 then begin
     if n_elements(data) eq 1 and data[0] eq 0 then MODFITS,file,0,head else $
       MWRFITS,data,file,head,/create     
  endif

  ;; FITS file does NOT exist and data given
  if info.exists eq 0 and n_elements(data) gt 1 then MWRFITS,data,file,head,/create

  ;; Create stand-alone header and update resource file
  ;; Load the resource file
  READLINE,rfile,rlines,count=nlines,/noblank,comment='#'
  arr = strtrim(strsplitter(rlines,'=',/extract),2)
  names = reform(arr[0,*])
  vals = reform(arr[1,*])
  rstr = create_struct(names[0],vals[0])
  if nlines gt 1 then for k=1,nlines-1 do rstr=create_struct(rstr,names[k],vals[k])
  ;; Does a header already exist
  if tag_exist(rstr,'HEADER') then newhfile=rstr.header else newhfile=file+'.head'
  WRITELINE,newhfile,head
  ;; Update resource file
  newrlines = rlines
  if tag_exist(rstr,'HEADER') then begin
    rind = where(strupcase(names) eq 'HEADER',nrind)
    newrlines[rind] = 'HEADER = '+newhfile
  endif else begin
    PUSH,newrlines,'HEADER = '+newhfile
  endelse
  WRITELINE,rfile,newrlines

;; No resource file, regular FITS file
endif else begin
  MWRFITS,data,file,head,/create
endelse

;stop

end
