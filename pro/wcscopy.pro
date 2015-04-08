pro wcscopy,file1,file2,stp=stp,distortion=distortion

;+
;
; WCSCOPY
;
; Copy WCS from file1 to file2
;
; INPUTS:
;  file1   The FITS file with the reference WCS that is to be copied
;  file2   The FITS file to be updated.
;  /distortion  Include the distortion terms (for TNX and TPV only for now).
;  /stp    Stop at the end of the program.
;
; OUTPUTS:
;  The WCS of file1 will be copied to file2.
;
; USAGE:
;  IDL>wcscopy,file1,file2
;
; By D. Nidever   Jun 2008
;-

nfile1 = n_elements(file1)
nfile2 = n_elements(file2)

; Not enough inputs
if (nfile1 eq 0 or nfile2 eq 0) then begin
  print,'Syntax - wcscopy,file1,file2,distortion=distortion'
  return
endif

test1 = file_test(file1)
test2 = file_test(file2)

; Files don't exist
if (test1 eq 0) then begin
  print,'FILE ',file1,' NOT FOUND'
  return
endif
if (test2 eq 0) then begin
  print,'FILE ',file2,' NOT FOUND'
  return
endif

; Same files
if file_search(file1,/fully) eq file_search(file2,/fully) then begin
  print,file1,' and ',file2,' ARE THE SAME FILE'
  return
endif


; Loading FILE1
FITS_READ,file1,im1,head1,message=error1

if (error1 ne '') then begin
  print,'ERROR reading in ',file1
  print,error1
  return
endif

; Extracting WCS in FILE1
EXTAST,head1,astr1

; Get distortion terms
if keyword_set(distortion) then begin
  ctype1 = sxpar(head1,'CTYPE1',/silent)
  dum = strsplit(ctype1,'-',/extract)
  wcstype = dum[1]
  case wcstype of  ; WCS type
  'TNX': wcs = hdr2wcstnx(head1)
  'TPV': wcs = hdr2wcstpv(head1)
  else: print,'No distortion terms for '+wcstype
  endcase
endif

if (n_elements(astr1) eq 0) then begin
  print,'NO WCS in ',file1
  return
endif

; Loading FILE2
FITS_READ,file2,im2,head2,message=error2

if (error2 ne '') then begin
  print,'ERROR reading in ',file2
  print,error2
  return
endif

; Putting WCS into FILE2
if n_elements(wcs) eq 0 then begin
  PUTAST,head2,astr1
endif else begin
  case wcstype of  ; WCS type
  'TNX': wcstnx2hdr, head2, wcs
  'TPV': wcstpv2hdr, head2, wcs
  else: print,'No distortion terms for '+wcstype
  endcase
endelse


; Updating FILE2
print,'Updating ',file2
MODFITS,file2,0,head2,errmsg=errmsg

if (n_elements(errmsg) ne 0) then begin
  print,'ERROR updating ',file2
  print,errmsg
endif

;stop

if keyword_set(stp) then stop

end
