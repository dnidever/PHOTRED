;+
;
; PHOTRED_DEREDDEN_SINGLE
;
; This dereddens photometry for photred
; See deredden.pro
;
; INPUTS:
;  longfile     Catalog filename to deredden.
;  extfile      Extinction file.
;  =todered     Magnitudes and colors to deredden.
;  =toextadd    Extinction values to add.
;  =catformat   Output table format.
;  /redo        Redo files that were already done.
;  /stp          Stop at the end of the program.
;
; OUTPUTS:
;  The final calibrated and dereddened photometry files
;
; By D.Nidever  Mar 2008
;-
pro photred_deredden_single,longfile,extfile,todered=todered,toextadd=toextadd,$
                            catformat=catformat,redo=redo,stp=stp

  CD,current=curdir

  if n_elements(catformat) eq 0 then catformat='FITS'

  ; Check that DUST_DIR is properly defined in the shell
  dustdir = GETENV('DUST_DIR')
  if (dustdir eq '') then begin
    print,'DUST_DIR TO SCHLEGEL MAPS NOT DEFINED IN SHELL'
    return
  endif

  ; Load the "extinction" file
  ;--------------------------
  print,'Loading the >>extinction<< file'
  ; Load the extinction values
  ; EXTRATIO is A(filter)/E(B-V), so A(filter) = EXTRATIO * E(B-V)
  extstr = IMPORTASCII(extfile,fieldname=['filter','extratio'],comment='#',/noprint)
  nextstr = n_elements(extstr)
  if (nextstr eq 0) then begin
    print,'>>extinction<< file is EMPTY'
    return
  endif


  ; What magnitudes and colors are we dereddening
  ;----------------------------------------------
  print,''
  print,'Getting the MAGNITUDES and COLORS to deredden'
  ntodered = n_elements(todered)
  ; Some mags/colors to deredden
  if (ntodered gt 0) then begin
    todered = strcompress(todered,/remove_all)    ; remove all whitespace
    toderedarr = strsplit(todered,',',/extract)   ; split
    ntoderedarr = n_elements(toderedarr)

    ; Structure type
    dum = {type:0,input:'',mag1:'',mag2:'',extratio1:0.0,extratio2:0.0}

    ; Loop through the mags/colors to deredden
    for i=0,ntoderedarr-1 do begin
      ideredstr = dum
      ideredstr.input = toderedarr[i]

      ; Is it a color or magnitude
      minus = strpos(toderedarr[i],'-')
      ; Magnitude
      if (minus eq -1) then begin
        ideredstr.type = 1
        ideredstr.mag1 = strtrim(toderedarr[i],2)

        ; Get extinction value
        gdext = where(extstr.filter eq ideredstr.mag1,ngdext)
        if (ngdext eq 0) then begin
          print,'TODERED - ',toderedarr[i],'  NO EXTINCTION VALUE FOR FILTER ',ideredstr.mag1
          goto,EXT_BOMB
        endif
        ideredstr.extratio1 = extstr[gdext[0]].extratio
      
        ; Add to structure
        PUSH,toderedstr,ideredstr

      ; Color
      endif else begin
        ideredstr.type = 0

        magarr = strsplit(toderedarr[i],'-',/extract)
        magarr = strtrim(magarr,2)
        nmagarr = n_elements(magarr)

        ; Can only have two magnitudes in the color
        if (nmagarr gt 2) then begin
          print,'TODERED - ',ideredstr.input,'  HAS MORE THAN 2 MAGNITUDES IN THE COLOR'
          goto,EXT_BOMB
        endif

        ideredstr.mag1 = magarr[0]
        ideredstr.mag2 = magarr[1]

        ; Get extinction values, Magnitude1
        gdext1 = where(extstr.filter eq ideredstr.mag1,ngdext1)
        if (ngdext1 eq 0) then begin
          print,'TODERED - ',ideredstr.input,'  NO EXTINCTION VALUE FOR FILTER ',ideredstr.mag1
          goto,EXT_BOMB
        endif
        ideredstr.extratio1 = extstr[gdext1[0]].extratio

        ; Get extinction values, Magnitude2
        gdext2 = where(extstr.filter eq ideredstr.mag2,ngdext2)
        if (ngdext2 eq 0) then begin
          print,'TODERED - ',ideredstr.input,'  NO EXTINCTION VALUE FOR FILTER ',ideredstr.mag2
          goto,EXT_BOMB
        endif
        ideredstr.extratio2 = extstr[gdext2[0]].extratio

        ; Add to structure
        PUSH,toderedstr,ideredstr

      endelse

      EXT_BOMB:

    endfor  ; todered loop

    ; Print out magnitudes/colors to deredden
    ntoderedstr = n_elements(toderedstr)
    for i=0,ntoderedstr-1 do begin
      ; Print header
      if (i eq 0) then begin
        print,''
        print,'Magnitudes/Colors to Deredden'
        print,'------------------------'
        print,'MAG/COLOR  A(M)/E(M1-M2)'
        print,'------------------------'
      endif
      tstr = toderedstr[i]

      ; Magnitude
      if (tstr.type eq 1) then begin
        print,'',tstr.mag1,tstr.extratio1,format='(A3,A-8,F8.3)'
      
      ; Color
      endif else begin
        red = tstr.extratio1-tstr.extratio2
        print,'',tstr.input,red,format='(A3,A-8,F8.3)'
      endelse

      ; Print footer
      if (i eq (ntoderedstr-1)) then print,'------------------------'
    endfor
  
  ; No mags/colors to deredden
  endif else begin
    print,'NO MAGNITUDES/COLORS to DEREDDEN'
  endelse


  ; Extinction values to add
  ;----------------------------
  print,''
  print,'Getting EXTINCTION values to add'
  ntoextadd = n_elements(toextadd)
  ; Some extinction values to add
  if (ntoextadd gt 0) then begin
    toextadd = strcompress(toextadd,/remove_all)    ; remove all whitespace
    toextaddarr = strsplit(toextadd,',',/extract)   ; split
    ntoextaddarr = n_elements(toextaddarr)

    ; Structure type
    dum = {type:0,input:'',mag1:'',mag2:'',extratio1:0.0,extratio2:0.0}

    ; Loop through the mags/colors to deredden
    for i=0,ntoextaddarr-1 do begin
      iextaddstr = dum
      iextaddstr.input = toextaddarr[i]

      ; Is it a color or magnitude
      minus = strpos(toextaddarr[i],'-')
      ; Magnitude
      if (minus eq -1) then begin
        iextaddstr.type = 1
        iextaddstr.mag1 = strtrim(toextaddarr[i],2)

        ; Get extinction value
        gdext = where(extstr.filter eq iextaddstr.mag1,ngdext)
        if (ngdext eq 0) then begin
          print,'TODERED - ',toextaddarr[i],'  NO EXTINCTION VALUE FOR FILTER ',iextaddstr.mag1
          goto,EXT_BOMB2
        endif
        iextaddstr.extratio1 = extstr[gdext[0]].extratio
      
        ; Add to structure
        PUSH,toextaddstr,iextaddstr

      ; Color
      endif else begin
        iextaddstr.type = 0

        magarr = strsplit(toextaddarr[i],'-',/extract)
        magarr = strtrim(magarr,2)
        nmagarr = n_elements(magarr)

        ; Can only have two magnitudes in the color
        if (nmagarr gt 2) then begin
          print,'TODERED - ',iextaddstr.input,'  HAS MORE THAN 2 MAGNITUDES IN THE COLOR'
          goto,EXT_BOMB2
        endif

        iextaddstr.mag1 = magarr[0]
        iextaddstr.mag2 = magarr[1]

        ; Get extinction values, Magnitude1
        gdext1 = where(extstr.filter eq iextaddstr.mag1,ngdext1)
        if (ngdext1 eq 0) then begin
          print,'TODERED - ',iextaddstr.input,'  NO EXTINCTION VALUE FOR FILTER ',iextaddstr.mag1
          goto,EXT_BOMB2
        endif
        iextaddstr.extratio1 = extstr[gdext1[0]].extratio

        ; Get extinction values, Magnitude2
        gdext2 = where(extstr.filter eq iextaddstr.mag2,ngdext2)
        if (ngdext2 eq 0) then begin
          print,'TODERED - ',iextaddstr.input,'  NO EXTINCTION VALUE FOR FILTER ',iextaddstr.mag2
          goto,EXT_BOMB2
        endif
        iextaddstr.extratio2 = extstr[gdext2[0]].extratio

        ; Add to structure
        PUSH,toextaddstr,iextaddstr

      endelse

      EXT_BOMB2:

    endfor  ; todered loop

    ; Print out magnitudes/colors to deredden
    ntoextaddstr = n_elements(toextaddstr)
    for i=0,ntoextaddstr-1 do begin
      ; Print header
      if (i eq 0) then begin
        print,''
        print,'Extinctions to add'
        print,'------------------------'
        print,'MAG/COLOR  A(M)/E(M1-M2)'
        print,'------------------------'
      endif
      tstr = toextaddstr[i]
      ; Magnitude
      if (tstr.type eq 1) then begin
        print,'',tstr.mag1,tstr.extratio1,format='(A3,A-8,F8.3)'
      ; Color
      endif else begin
        red = tstr.extratio1-tstr.extratio2
        print,'',tstr.input,red,format='(A3,A-8,F8.3)'
      endelse
      ; Print footer
      if (i eq (ntoextaddstr-1)) then print,'------------------------'
    endfor

  ; No extinctions to add
  endif else begin
    print,'NO EXTINCTIONS TO ADD'
  endelse


  ; Only adding E(B-V)
  ntoderedstr = n_elements(toderedstr)
  ntoextaddstr = n_elements(toextaddstr)
  if (ntoderedstr eq 0) and (ntoextaddstr eq 0) then print,'ONLY ADDING E(B-V)'



  file = FILE_BASENAME(longfile)
  filedir = FILE_DIRNAME(longfile)
  ending = first_el(strsplit(file,'.',/extract),/last)
  base = FILE_BASENAME(file,'.'+ending)  

  print,''
  print,'============================='
  print,'DEREDDENING ',file
  print,'============================='
  print,''
  print,systime(0)
  
  CD,filedir

  ; Make sure it exists
  test = FILE_TEST(file)
  if test eq 1 then nlines = FILE_LINES(file)
  if (test eq 0) or (nlines eq 0) then begin
    if test eq 0 then print,file,' NOT FOUND'
    if test eq 1 and nlines eq 0 then print,file,' HAS 0 LINES'
    CD,curdir
    return
  endif

  ; Load the data file
  ;--------------------
  print,'Loading file'
  str = PHOTRED_READFILE(file,meta,count=nstr)

  ; Checking that we've got coordinates
  tags = TAG_NAMES(str)
  raind = where(tags eq 'RA',nraind)
  decind = where(tags eq 'DEC',ndecind)
  if (nraind eq 0) or (ndecind) eq 0 then begin
    print,file,' DOES NOT HAVE RA/DEC'
    CD,curdir
    return
  endif

  ;; Make new structure schema
  ;;---------------------------
  schema = str[0]
  struct_assign,{dum:''},schema
  schema = create_struct(schema,'EBV',0.0)
  ; ADDING Dereddened Magnitudes/Colors
  for j=0,ntoderedstr-1 do begin
    tstr = toderedstr[j]
    ; Magnitude
    if (tstr.type eq 1) then begin
      ; Check that the magnitude exists
      magind = where(tags eq strupcase(tstr.mag1)+'MAG',nmagind)
      ; Adding the tag to the structure
      if nmagind gt 0 then begin
        newtagname = strupcase(tstr.mag1)+'0'
        newtagname = IDL_VALIDNAME(newtagname,/convert_all)    ; make sure it's a valid IDL name
        gdnewtagname = where(tags eq newtagname,ngdnewtagname)
        if (ngdnewtagname eq 0) then schema = create_struct(schema,newtagname,0.0)
      endif
    ; Color
    endif else begin
      ; Check that magnitude1 exists
      mag1ind = where(tags eq strupcase(tstr.mag1)+'MAG',nmag1ind)
      ; Check that magnitude2 exists
      mag2ind = where(tags eq strupcase(tstr.mag2)+'MAG',nmag2ind)
      ; Adding the tag to the structure
      if nmag1ind gt 0 and nmag2ind gt 0 then begin
        newtagname = strupcase(tstr.mag1)+strupcase(tstr.mag2)+'0'
        newtagname = IDL_VALIDNAME(newtagname,/convert_all)    ; make sure it's a valid IDL name
        gdnewtagname = where(tags eq newtagname,ngdnewtagname)
        if (ngdnewtagname eq 0) then schema = create_struct(schema,newtagname,0.0)
      endif
    endelse
  endfor  ; dereddened mag loop
  ; ADDING EXTINCTIONS
  for j=0,ntoextaddstr-1 do begin
    tstr = toextaddstr[j]
    ; Magnitude
    if (tstr.type eq 1) then begin
      ; Check that the magnitude exists
      magind = where(tags eq strupcase(tstr.mag1)+'MAG',nmagind)
      if (nmagind gt 0) then begin
        ; Adding the tag to the structure
        newtagname = 'A'+strupcase(tstr.mag1)
        newtagname = IDL_VALIDNAME(newtagname,/convert_all)    ; make sure it's a valid IDL name
        gdnewtagname = where(tags eq newtagname,ngdnewtagname)
        if (ngdnewtagname eq 0) then schema = create_struct(schema,newtagname,0.0)
      endif
    ; Color
    endif else begin
      ; Check that magnitude1 exists
      mag1ind = where(tags eq strupcase(tstr.mag1)+'MAG',nmag1ind)
      ; Check that magnitude2 exists
      mag2ind = where(tags eq strupcase(tstr.mag2)+'MAG',nmag2ind)
      if nmag1ind gt 0 and nmag2ind gt 0 then begin
        ; Adding the tag to the structure
        newtagname = 'E'+strupcase(tstr.mag1)+strupcase(tstr.mag2)
        newtagname = IDL_VALIDNAME(newtagname,/convert_all)    ; make sure it's a valid IDL name
        gdnewtagname = where(tags eq newtagname,ngdnewtagname)
        if (ngdnewtagname eq 0) then schema = create_struct(schema,newtagname,0.0)
      endif
    endelse
  endfor  ; extinction loop

  ;; Creating new structure and copying over existing data
  orig = str
  str = REPLICATE(schema,nstr)
  STRUCT_ASSIGN,orig,str,/nozero
  undefine,orig
  tags = TAG_NAMES(str)


  ; Getting E(B-V) from the Schlegel maps
  ;----------------------------------------
  ; Convert to galactic coordinates
  GLACTC,str.ra,str.dec,2000.,glon,glat,1,/degree
  ; Get the E(B-V) reddening from the Schlegel maps
  ebv = DUST_GETVAL(glon,glat,/interp,/noloop)
  str.EBV = EBV



  ; There can be problems with the dereddening if there
  ; are multiple observations per band.  These will have
  ; names of M1, M2, etc. and won't match the filter names.

  ;---------------------------------------
  ; ADDING Dereddened Magnitudes/Colors
  ;---------------------------------------
  ; Loop through the toderedstr structure
  for j=0,ntoderedstr-1 do begin

    tstr = toderedstr[j]

    ; Magnitude
    ;----------
    if (tstr.type eq 1) then begin

      ; Check that the magnitude exists
      ;magind = where(tags eq strupcase(tstr.mag1),nmagind)
      magind = where(tags eq strupcase(tstr.mag1)+'MAG',nmagind)
      if (nmagind eq 0) then begin
        print,'NO Magnitude ',tstr.mag1
        goto,TODERED_BOMB
      endif

      ; Getting the magnitude
      mag = str.(magind[0])

      ; Getting tag name
      newtagname = strupcase(tstr.mag1)+'0'
      newtagname = IDL_VALIDNAME(newtagname,/convert_all)    ; make sure it's a valid IDL name
      deredmagind = where(tags eq newtagname,nderedmagind)

      ; Dereddening
      red = tstr.extratio1 * str.ebv
      deredmag = mag - red
      ; Only deredden stars with decent magnitudes
      bdmag = where(mag gt 50.,nbdmag)
      if (nbdmag gt 0) then deredmag[bdmag] = 99.99
      ; Adding to the structure
      str.(deredmagind) = deredmag

    ; Color
    ;------
    endif else begin

      ; Check that magnitude1 exists
      ;mag1ind = where(tags eq strupcase(tstr.mag1),nmag1ind)
      mag1ind = where(tags eq strupcase(tstr.mag1)+'MAG',nmag1ind)
      if (nmag1ind eq 0) then begin
        print,'NO Magnitude ',tstr.mag1
        goto,TODERED_BOMB
      endif

      ; Check that magnitude2 exists
      ;mag2ind = where(tags eq strupcase(tstr.mag2),nmag2ind)
      mag2ind = where(tags eq strupcase(tstr.mag2)+'MAG',nmag2ind)
      if (nmag2ind eq 0) then begin
        print,'NO Magnitude ',tstr.mag2
        goto,TODERED_BOMB
      endif

      ; Getting the magnitude
      mag1 = str.(mag1ind[0])
      mag2 = str.(mag2ind[0])

  ; If it is E(B-V) then DO NOT redo it.!!

      ; Getting tag name
      newtagname = strupcase(tstr.mag1)+strupcase(tstr.mag2)+'0'
      newtagname = IDL_VALIDNAME(newtagname,/convert_all)    ; make sure it's a valid IDL name
      deredcolind = where(tags eq newtagname,nderedmagind)

      ; Dereddening
      red1 = tstr.extratio1 * str.ebv
      red2 = tstr.extratio2 * str.ebv
      deredcol = (mag1 - red1) - (mag2 - red2)
      ; Only deredden stars with decent magnitudes
      bdmag = where(mag1 gt 50. or mag2 gt 50.,nbdmag)
      if (nbdmag gt 0) then deredcol[bdmag] = 99.99
      ; Adding to the structure
      str.(deredcolind) = deredcol
    endelse

    TODERED_BOMB:
  endfor


  ;---------------------------------------
  ; ADDING EXTINCTIONS
  ;---------------------------------------
  ; Loop through the toderedstr structure
  for j=0,ntoextaddstr-1 do begin
    tstr = toextaddstr[j]

    ; Magnitude
    ;----------
    if (tstr.type eq 1) then begin

      ; Check that the magnitude exists
      ;magind = where(tags eq strupcase(tstr.mag1),nmagind)
      magind = where(tags eq strupcase(tstr.mag1)+'MAG',nmagind)
      if (nmagind eq 0) then begin
        print,'NO Magnitude ',tstr.mag1
        goto,TODERED_BOMB2
      endif

      ; Getting tag name
      newtagname = 'A'+strupcase(tstr.mag1)
      newtagname = IDL_VALIDNAME(newtagname,/convert_all)    ; make sure it's a valid IDL name
      extind = where(tags eq newtagname,nextind)

      ; Extinction
      str.(extind) = tstr.extratio1 * str.ebv

    ; Color
    ;------
    endif else begin

      ; Check that magnitude1 exists
      ;mag1ind = where(tags eq strupcase(tstr.mag1),nmag1ind)
      mag1ind = where(tags eq strupcase(tstr.mag1)+'MAG',nmag1ind)
      if (nmag1ind eq 0) then begin
        print,'NO Magnitude ',tstr.mag1
        goto,TODERED_BOMB2
      endif

      ; Check that magnitude2 exists
      ;mag2ind = where(tags eq strupcase(tstr.mag2),nmag2ind)
      mag2ind = where(tags eq strupcase(tstr.mag2)+'MAG',nmag2ind)
      if (nmag2ind eq 0) then begin
        print,'NO Magnitude ',tstr.mag2
        goto,TODERED_BOMB2
      endif

      ; Getting tag name
      newtagname = 'E'+strupcase(tstr.mag1)+strupcase(tstr.mag2)
      newtagname = IDL_VALIDNAME(newtagname,/convert_all)    ; make sure it's a valid IDL name
      extind = where(tags eq newtagname,nextind)

      ; Dereddening
      str.(extind) = tstr.extratio1*str.ebv - tstr.extratio2*str.ebv
    endelse

    TODERED_BOMB2:

  endfor

  ; WRITE new file
  ;------------------
  ; Output the structure to the DERED file
  deredfile = filedir+'/'+base+'.dered'
  print,'New file with extinctions is: ',deredfile
  if catformat eq 'ASCII' then begin
    PRINTSTR,str,deredfile
    if n_elements(meta) gt 0 then PRINTSTR,meta,deredfile+'.meta'
  endif else begin
    MWRFITS,str,deredfile,/create
    if n_elements(meta) gt 0 then MWRFITS,meta,deredfile,/silent
  endelse

  ; Check that the file DERED file is there
  deredtest = FILE_TEST(deredfile)
  if (deredtest eq 0) then begin
    print,deredfile,' NOT FOUND'
  endif

  CD,curdir

  if keyword_set(stp) then stop

end
