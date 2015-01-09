pro stdred_copy,inpdir,inpstand,stp=stp

;+
;
; STDRED_COPY
;
; This copes all of the standard star frames (i.e. SA98, SA110, etc.)
; from the varous nights' directories into the "standards/" directory.
; Only directory names that start with "n" are used, or the names in
; the "inpdir" input.  By default only the four Geisler Washington 
; standard star fields SA98, SA110, SA114, and NGC3680 standard are
; copied unless an additional list of standard star field names is
; input with "inpstand".
;
; INPUTS:
;  inpdir       An input list of directories.  This can be an array of
;                 directory names, a string with wildcards (i.e. "n*"),
;                 or a filename (with a @ prepended) that contains
;                 directory names (i.e. "@dirlist.txt").
;  inpstand     An input list of standard star names.  This can be an
;                 array of  names, or a filename (with a @ prepended)
;                 that contains the standard star names (i.e. "@standards.txt").
;                 The four Geisler fields (SA98, SA110, SA114, and
;                 NGC3680) are always included in the standard star
;                 field list.
;  /stp         Stop at the end of the program
;
; OUTPUTS:
;  All standard star field frames are copied to the "standards/" directory.
;
; USAGE:
;  IDL>stdred_copy,'@dirlist.txt','@standards.txt'
;
; By D. Nidever   May 2008
;-

count=0
ndirs = 0
CD,current=curdir
lencurdir = strlen(curdir)

; If there is no "standards/" directory then make it
if FILE_TEST('standards',/directory) eq 0 then begin
  print,'Making "standards" directory'
  FILE_MKDIR,'standards'
endif

; Using input file
ninpdir = n_elements(inpdir)
if ninpdir gt 0 then begin

  ; Load the input
  LOADINPUT,inpdir,input,count=ninput

  ; We've got something
  if (ninput gt 0) then begin
    dirs = FILE_SEARCH(input,/test_directory,count=ndirs)
 
    ; No directories that exist
    if (ndirs eq 0) then begin
      print,'NO directories that exist'
      print,dirs,' NOT FOUND'
      return
    endif

    print,strtrim(ndirs,2),' directories found'
    print,dirs

  ; No input
  endif else begin
    print,'NO INPUT'
    return
  endelse

; Automatically get the directory names (i.e. n1, n2, etc.)
endif else begin

  print,'Automatically getting the directory names (i.e. n1, n2, etc.)'

  dirs = FILE_SEARCH('n*',/test_directory,count=ndirs)

  ; No directories
  if (ndirs eq 0) then begin
    print,'NO "n*" DIRECTORIES FOUND'
    return 
  endif

  print,strtrim(ndirs,2),' directories found'
  print,dirs

endelse


; Loading standards list
ninpstand = n_elements(inpstand)
nstandlist = 0
if (ninpstand gt 0) then begin
  LOADINPUT,inpstand,standlist,count=nstandlist

  ; Remove empty strings
  if nstandlist gt 0 then begin
    standlist = strlowcase(strtrim(standlist,2))
    gd = where(standlist ne '',ngd)
    if ngd eq 0 then undefine,standlist
    if ngd gt 0 then standlist=standlist[gd]
    nstandlist = n_elements(standlist)
  endif

  if (nstandlist gt 0) then begin
    print,''
    print,'Standards list loaded. ',strtrim(nstandlist,2),' standard(s) found'
    printline,standlist
    print,''
  endif

endif


; Loop through the directories and find any FITS files
;-----------------------------------------------------
FOR i=0,ndirs-1 do begin

  print,''
  print,'==============================================='
  print,'Copying Standard Star Frames from ',dirs[i]
  print,'==============================================='

  files = FILE_SEARCH(dirs[i],'*.fits',count=nfiles,/fully)
  files2 = strmid(files,lencurdir+1,200) ; relative
  base = FILE_BASENAME(files)

  ; Remove any raw or calibration files
  print,'Removing any Raw/, zero, bias, flat, focus or test frames'

  ; raw files
  raw = stregex(files2,'Raw/',/boolean)
  ; zero in name?
  zeroname = stregex(files2,'zero',/fold_case,/boolean)
  ; bias in name?
  biasname = stregex(files2,'bias',/fold_case,/boolean)
  ; flat in name?
  flatname = stregex(files2,'flat',/fold_case,/boolean)
  ; focus in name?
  focusname = stregex(files2,'focus',/fold_case,/boolean)
  ; test in name?
  testname = stregex(files2,'test',/fold_case,/boolean)

  bad = raw + zeroname + biasname + flatname + focusname + testname
  gd = where(bad eq 0,ngd)
  if (ngd eq 0) then goto,BOMB

  ; Only want good ones
  files = files[gd]
  files2 = files2[gd]
  nfiles = ngd

  print,strtrim(nfiles,2),' FITS files found'

print,'-------------------------------------------------------------------------------------------------------------------------'
print,'   TYPE            FILENAME                NEW FILENAME           OBJECT                FILTER    EXPTIME       DATE-OBS'
print,'-------------------------------------------------------------------------------------------------------------------------'


  ; Loop through the FITS files
  For j=0,nfiles-1 do begin

    longfile = files[j]
    relfile = files2[j]
    file = FILE_BASENAME(longfile)
    base = FILE_BASENAME(longfile,'.fits')

    head = HEADFITS(longfile)
    object = SXPAR(head,'OBJECT',count=nobject)
    exptime = -1.0
    exptime = SXPAR(head,'EXPTIME')
    filter = SXPAR(head,'FILTER')
    ut = SXPAR(head,'DATE-OBS')

    if (nobject eq 0) then begin
      print,'OBJECT parameter not found'
      goto,BOMB
    endif


    ; Make sure it's not a calibration frame

    ; zero in object string
    zeroobj = stregex(object,'zero',/fold_case,/boolean)
    ; bias in object string
    biasobj = stregex(object,'bias',/fold_case,/boolean)
    ; flat in object string
    flatobj = stregex(object,'flat',/fold_case,/boolean)
    ; twilight
    twiobj = stregex(object,'twil',/fold_case,/boolean)
    ; sky
    skyobj = stregex(object,'sky',/fold_case,/boolean)
    ; pointing
    pointobj = stregex(object,'pointing',/fold_case,/boolean)
    ; focus
    focusobj = stregex(object,'focus',/fold_case,/boolean)
    ; test
    testobj = stregex(object,'test',/fold_case,/boolean)
    bad2 = zeroobj+biasobj+flatobj+twiobj+skyobj+pointobj+focusobj+testobj

    ;if bad2 eq 1 then print,files2[j],' is a calibration frame. NOT standard star frame.'


    ; Check the four Washington Geisler fields
    stand = where(stregex(object,'sa98',/fold_case,/boolean) OR $
                  stregex(object,'sa 98',/fold_case,/boolean) OR $
                  stregex(object,'sa-98',/fold_case,/boolean) OR $
                  stregex(object,'sa110',/fold_case,/boolean) OR $
                  stregex(object,'sa 110',/fold_case,/boolean) OR $
                  stregex(object,'sa-110',/fold_case,/boolean) OR $
                  stregex(object,'sa114',/fold_case,/boolean) OR $
                  stregex(object,'sa 114',/fold_case,/boolean) OR $
                  stregex(object,'sa-114',/fold_case,/boolean) OR $
                  stregex(object,'n3680',/fold_case,/boolean) OR $
                  stregex(object,'ngc3680',/fold_case,/boolean) OR $
                  stregex(object,'ngc 3680',/fold_case,/boolean) OR $
                  stregex(object,'ngc-3680',/fold_case,/boolean),nstand)

    ; Check the input standards list
    ninputstand = 0
    for k=0,nstandlist-1 do begin
      stand2 = where(stregex(object,strlowcase(standlist[k]),/fold_case,/boolean),nstand2)
      ninputstand = ninputstand + nstand2
    end


    ; Calibration frame
    newfile = '---'
    com = 'Non-standard'
    if zeroobj eq 1 then com='Zero'
    if biasobj eq 1 then com='Bias'
    if flatobj eq 1 then com='Flat'
    if twiobj eq 1 then com='Twilight'
    if skyobj eq 1 then com='Sky'
    if pointobj eq 1 then com='Pointing'
    if focusobj eq 1 then com='Focus'
    if testobj eq 1 then com='Test'

    ; Standard star frame
    if ((nstand gt 0 OR ninputstand gt 0) and bad2 eq 0) then begin
      com = 'Standard'
      newfile = 'standards/'+file
    endif

    ; Print the information
    fmt = '(A-15,A-25,A-25,A-25,A-6,F8.2,A16)'
    print,com,relfile,newfile,object,filter,exptime,ut,format=fmt

    ; This is a standard star frame, COPY IT
    if ((nstand gt 0 OR ninputstand gt 0) and bad2 eq 0) then begin
     ; print,'Copying ',files2[j],' to "standards"'
      FILE_COPY,files[j],'standards',/overwrite,/allow_same
      count++
    endif

    ;stop


  End ; files loop



  BOMB:


END

print,''
print,strtrim(count,2),' Files copied to "standards"'
print,''

if keyword_set(stp) then stop

end
