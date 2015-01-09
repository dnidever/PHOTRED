function photred_getfilter,file,stp=stp,numeric=numeric,noupdate=noupdate,$
             silent=silent,filtname=filtname0,error=error,fold_case=fold_case

;+
;
; PHOTRED_GETFILTER
;
; This gets filter information for an image
; using the "filter" file.
; The "short" filter names that are returned
; are not necessarily "unique".  The filter names
; "I c6028", "T", "T2" all have a short filter
; name of "T".
;
; If a filter is not found in the filter list "filters"
; then a new short name is created for it (based on the
; first word in the string) and added to the list.
; 
; INPUTS:
;  file      FITS filename
;  /numeric  Return a numeric value instead of letter.
;  =filtname Input the filter name explicitly instead of giving the
;              FITS filename.
;  /noupdate Don't update the "filters" file if the filter is not found.
;  /fold_case  Ignore case, the default is for it to be case sensitive.
;                This is needed for checking IDL structure tags.
;  /silent   Don't print anything to the screen.
;  /stp      Stop at the end of the program.
;
; OUTPUTS:
;  The short filter name is output.
;  =error    The error message if an error occured.
;
; USAGE:
;  IDL>filter = photred_getfilter(fitsfile,numeric=numeric,noupdate=noupdate,
;                                 filtname=filtname,silent=silent,stp=stp,fold_case=fold_case,error=error)
;
; By D.Nidever  February 2008
;-

COMMON photred,setup

undefine,error

nfile = n_elements(file)
nfiltname = n_elements(filtname0)
; Not enough inputs
if nfile eq 0 and nfiltname eq 0 then begin
  print,'Syntax - filter = photred_getfilter(fitsfile,numeric=numeric,noupdate=noupdate,'
  print,'                                    filtname=filtname,silent=silent,stp=stp,error=error,fold_case=fold_case)'
  error = 'Not enough inputs'
  return,''
endif

; More than one FITS filename input
if (nfile gt 1) then begin
  filters = strarr(nfile)
  for i=0,nfile-1 do filters[i] = photred_getfilter(file[i],numeric=numeric,noupdate=noupdate,fold_case=fold_case)
  return,filters
endif

; More than one filter name input
; ONLY if no FITS filename input
if (nfile eq 0 and nfiltname gt 1) then begin
  filters = strarr(nfiltname)
  for i=0,nfiltname-1 do filters[i] = photred_getfilter(filtname=filtname0[i],numeric=numeric,$
                                                        noupdate=noupdate,silent=silent,fold_case=fold_case)
  return,filters
endif

; Are we case sensitive
if n_elements(fold_case) eq 0 then fold_case=0 ; by default, case insensitive

; Does the "filters" file exist?
test = file_test('filters')
if test eq 0 then begin
  scriptsdir = READPAR(setup,'SCRIPTSDIR')
  if scriptsdir eq '-1' or scriptsdir eq '0' then begin
    if not keyword_set(silent) then print,'NO SCRIPTSDIR'
    error = 'NO SCRIPTSDIR'
    return,''
  endif
  FILE_COPY,scriptsdir+'/filters','.',/overwrite
endif

; Load the filters
READLINE,'filters',lines
gd = where(strtrim(lines,2) ne '',ngd)
if ngd eq 0 then begin
  if not keyword_set(silent) then print,'NO FILTERS'
  error = 'NO FILTERS'
  return,''
endif
lines = lines[gd]
arr = strsplitter(lines,"'",/extract)
arr = strtrim(arr,2)
; Should be 2xN
; Removing blank lines
longnames = reform(arr[0,*])
shortnames = reform(arr[1,*])

; Get the filter information from the header
if (nfile gt 0) then begin

  ; Does it have the ".fits" ending
  ext = first_el(strsplit(file,'.',/extract),/last)
  if ext ne 'fits' then begin
    if not keyword_set(silent) then print,file,' IS NOT A FITS FILE'
    error = file+' IS NOT A FITS FILE'
    return,''
  endif

  ; Make sure the file exists
  test = file_test(file)
  if test eq 0 then begin
    if not keyword_set(silent) then print,file,' NOT FOUND'
    error = file+' NOT FOUND'
    return,''
  endif

  head = HEADFITS(file)

  ; Problem with header
  if n_elements(head) eq 1 and strtrim(head[0],2) eq '-1' then begin
    if not keyword_set(silent) then print,file,' - Problem loading FITS header'
    error = file+' - Problem loading FITS header'
    return,''
  endif

  filtname = SXPAR(head,'FILTER')
  filtname = strtrim(filtname,2)
  if strtrim(filtname,2) eq '0' then begin
    if not keyword_set(silent) then $
      print,'NO FILTER INFORMATION IN '+file+' HEADER'
    error = 'NO FILTER INFORMATION IN '+file+' HEADER'
    return,''
  endif

; Get the filter name from "filtname"
endif else begin
  filtname = strtrim(filtname0[0],2)
endelse


; Match filter
;ind = where(names eq filtname,nind)
;ind = first_el(where(strcmp(longnames,filtname,/fold_case) eq 1,nind))
ind = first_el(where(strcmp(longnames,filtname,fold_case=fold_case) eq 1,nind))
;if nind eq 0 then begin
;  ind = where(stregex(longnames,filtname,/boolean,/fold_case) eq 1,nind)
;endif


; Match found
if nind gt 0 then begin
  
  filter = shortnames[ind[0]]

  ; Numeric value
  if keyword_set(numeric) then begin
    ui = uniq(shortnames)
    snames = shortnames[ui]           ; unique short names
    nui = n_elements(ui)
    numnames = strtrim(lindgen(nui)+1,2)   ; numbers for the unique short names
 
    gg = where(snames eq filter,ngg)  ; which short name
    numname = numnames[gg[0]]
    return,numname
  endif


  if keyword_set(stp) then stop

  return,filter

; No match found
endif else begin

  if not keyword_set(silent) then print,'NO FILTER MATCH'


  ; Add it to the "filters" file
  if (not keyword_set(noupdate)) then begin

    ; The IRAF task is called "ccdsubset"
    ;# CCDSUBSET -- Return the CCD subset identifier.
    ;#
    ;# 1. Get the subset string and search the subset record file for the ID string.
    ;# 2. If the subset string is not in the record file define a default ID string
    ;#    based on the first word of the subset string.  If the first word is not
    ;#    unique append a integer to the first word until it is unique.
    ;# 3. Add the new subset string and identifier to the record file.
    ;# 4. Since the ID string is used to generate image names replace all
    ;#    nonimage name characters with '_'.
    ;#
    ;# It is an error if the record file cannot be created or written when needed.

    ; Get first word of the ID string
    dum = strsplit(filtname,' ',/extract)
    newshortname = reform(dum[0])

    ; Is this already a "short" filter name
    ;ind = first_el(where(strcmp(shortnames,newshortname,/fold_case) eq 1,nind))
    ind = first_el(where(strcmp(shortnames,newshortname,fold_case=fold_case) eq 1,nind))

    ; Yes, already a "short" name
    ; Append integer until unique
    if nind gt 0 then begin

      ;newshortname = newshortname+'_'

      ; Loop until we have a unique name
      flag=0
      integer=1
      WHILE (flag eq 0) do begin
        sinteger = strtrim(integer,2)
        ;ind = first_el(where(strcmp(shortnames,newshortname+sinteger,/fold_case) eq 1,nind))
        ind = first_el(where(strcmp(shortnames,newshortname+sinteger),fold_case=fold_case eq 1,nind))
 
        ; Unique
        if nind eq 0 then begin
          newshortname = newshortname+sinteger
          flag=1
        ; Not unique, increment integer
        endif else begin
          integer++
        endelse

      ENDWHILE
    endif

    ; Make sure the variable is okay
    newshortname = IDL_VALIDNAME(newshortname,/convert_all)
  
    ; Make sure it doesn't have any weird characters
    ; such as '*', '!', '$'
    newshortname = REPSTR(newshortname,'*','_')
    newshortname = REPSTR(newshortname,'!','_')
    newshortname = REPSTR(newshortname,'$','_')
    newshortname = REPSTR(newshortname,'__','_')
    newshortname = REPSTR(newshortname,'__','_')

    ; Add new filter to the "filters" file
    newline = "'"+filtname+"'     "+newshortname
    WRITELINE,'filters',newline,/append
    print,'Adding new filter name to "filters" list'
    print,'Filter ID string:  ',filtname
    print,'Filter short name: ',newshortname


    ; Numeric value
    if keyword_set(numeric) then begin

      ; Reload filters
      READLINE,'filters',lines
      gd = where(strtrim(lines,2) ne '',ngd)
      lines = lines[gd]
      arr = strsplitter(lines,"'",/extract)
      arr = strtrim(arr,2)
      ; Should be 2xN
      ; Removing blank lines
      longnames = reform(arr[0,*])
      shortnames = reform(arr[1,*])

      ui = uniq(shortnames)
      snames = shortnames[ui]           ; unique short names
      nui = n_elements(ui)
      numnames = strtrim(lindgen(nui)+1,2)   ; numbers for the unique short names
 
      gg = where(snames eq newshortname,ngg)  ; which short name
      numname = numnames[gg[0]]
      return,numname
    endif

  ; Don't update
  endif else begin
    error = 'NO FILTER MATCH'
    return,''
  endelse

  if keyword_set(stp) then stop

  return,newshortname

endelse

end
