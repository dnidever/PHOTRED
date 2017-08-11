PRO fakered_allframe_wrapper

  args = command_line_args()
  n_args = n_elements(args)
  img = ""
  scriptsdir = ""
  irafdir    = ""
  finditer   = ""
  detectprog = ""
  tile       = ""

  print,n_args, " arguments received!" 
  ; Check all arguments
  for i=0,n_args-1 do begin 
    print,"arg",i, ": ", args[i]
    ; Remove all quote symbols "
    args[i] = STRJOIN(STRSPLIT(args[i], '"', /EXTRACT), '')
    ; Get subarguments (separeted by commas)
    subargs = strsplit(args[i],",", /EXTRACT)
    n_subargs = n_elements(subargs)
    if n_subargs gt 2 and subargs[0] eq "allframe" then $
      ; If first argument is procedure "allframe", second arg must be the image
      img = subargs[1]
    ; Process remaining subargs (they should be KeyWords)
    for j=2, n_elements(subargs)-1 do begin
      ;print,"Param",j, ": ", subargs[j]
      ; Separate keywords to get name=value 
      kw = strsplit(subargs[j],"=", /EXTRACT)
      ; Check if there are two elements: name=value)
      n_kw = n_elements(kw)
      if n_kw eq 2 then begin
        ; Process each keyword
        case strtrim(kw[0],2) of
          "scriptsdir": scriptsdir = strtrim(kw[1],2)
          "irafdir"   : irafdir    = strtrim(kw[1],2) 
          "finditer"  : finditer   = strtrim(kw[1],2) 
          "detectprog": detectprog = strtrim(kw[1],2) 
          "tile"      : begin 
                          ; Get the tile type (ignore everything before ':')
                          tile_aux = strsplit(kw[1],":", /EXTRACT)
                          if n_elements(tile_aux) eq 2 then $
                            ; Remove extra '}'
                            tile = strtrim(strjoin(strsplit(tile_aux[1], '}', /EXTRACT), ''), 2)
                        end
          else        : print,"Ignoring: ", kw[1]
         endcase
      endif
    endfor
  endfor
      
  ; RUN ALLFRAME!!    
  print,"RUNNING: allframe,",img,",scriptsdir=",scriptsdir,", irafdir=",irafdir,     $
        ", finditer=", finditer,", detectprog=",detectprog,", tile={type:",tile,"}, /fake"
  allframe, img, scriptsdir=scriptsdir, irafdir=irafdir, finditer=finditer, detectprog=detectprog, tile={type:tile}, /fake


END
