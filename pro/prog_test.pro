;=================================================================

function which_find_routine, proname, _REF_EXTRA=_extra
; LOOKS FOR A MATCH BETWEEN ROUTINE INFORMATION AND
; AN IDL MODULE NAME...
; CLEVERLY, COMPILATION OF WHICH GUARANTEES THAT THERE WILL ALWAYS
; BE AT LEAST ONE PROCEDURE (WHICH) AND FUNCTION (WHICH_FIND_ROUTINE)...
compile_opt idl2, hidden
return, strmatch(routine_info(_EXTRA=_extra), proname, /FOLD_CASE)
end; which_find_routine

;=================================================================

function prog_test, proname, silent=silent, embedded=embedded
;+
; NAME:
;       PROG_TEST
;
; PURPOSE: 
;       This tests whether a specified program(s) exists.  It looks
;       at the already compiled programs and programs in the IDL
;       !path.
;       This will not find the program if it is imbedded in a file
;       by another name.
;
; CALLING SEQUENCE:
;       test = PRO_TEST(file)
;
; INPUTS:
;       FILE - file name to search for.  The suffix .pro will be
;              appended if not included.  This can be a list.
;
; KEYWORD PARAMETERS:
;       /SILENT  Don't print anything to the screen
;       /EMBEDDED  Search for embedded programs.  This takes a *long* time.
;
; OUTPUTS:
;       TEST -  1 if the program exists and 0 if it doesn't.
;
; COMMON BLOCKS:
;       None.
;
; RESTRICTIONS: 
;       The IDL !path is searched for file names that are simply the
;       module (in IDL documentation, "module" and "routine" are used
;       interchangeably) name with a ".pro" suffix appended to them.
;       A module stored inside a file whose name is different than the
;       module name (followed by a ".pro") will not be found UNLESS
;       that module happens to be the currently-resolved module!
;       E.g., if the module "pro test_proc" lives in a file named
;       "dumb_name.pro", then it will not be found:
; 
; PROCEDURES CALLED:
;       STRSPLIT(), WHICH_FIND_ROUTINE()
;
; EXAMPLES:
;       IDL>test = prog_test('program.pro')
;       IDL>print,test
;                  1
;
; MODIFICATION HISTORY:
;   30 May 2003  WHICH.PRO written by Tim Robishaw, Berkeley
;   17 Feb 2004  Fixed oddity where user tries to call a function as
;                if it were a procedure, thus listing the module in both
;                the Compiled Functions and Compiled Procedures list.
;   30 Oct 2007  Modified by D. Nidever to make prog_test.pro
;-

on_error, 2
resolve_routine, 'strsplit', /IS_FUN, /NO_RECOMPILE

if (N_params() lt 1) then begin
    message, 'syntax: test = prog_test( proname (suffix .pro assumed))', /INFO
    return,-1
endif

; Multiple files input
nproname = n_elements(proname)
if nproname gt 1 then begin
  test = lonarr(nproname)
  for i=0L,nproname-1 do test[i] = prog_test(proname[i],silent=silent)
  return,test
endif


; IF .PRO SUFFIX INCLUDED, DROP IT...
proname = strtrim(proname,2)
if strmatch(proname,'*.pro', /FOLD_CASE) $
  then proname = strmid(proname,0,strlen(proname)-4)

; SEARCH THE CURRENTLY-COMPILED PROCEDURES AND FUNCTIONS FIRST...
pindx = where(which_find_routine(proname),presolved)
findx = where(which_find_routine(proname,/FUNCTIONS),fresolved)

; IF PROCEDURE OR FUNCTION WAS FOUND, IS IT UNRESOLVED...
punresolved = total(which_find_routine(proname,/UNRESOLVED))
funresolved = total(which_find_routine(proname,/UNRESOLVED,/FUNCTIONS))

;if (presolved and not punresolved) OR $
;   (fresolved and not funresolved) then begin
;
;    ; THE PROCEDURE OR FUNCTION WAS FOUND...
;    resolved_routine = (presolved AND not fresolved) ? $
;      (routine_info(/SOURCE))[pindx].PATH : $
;      (routine_info(/SOURCE,/FUNCTIONS))[findx].PATH
;
;    print, 'Currently-Compiled Module '+strupcase(proname)+' in File:'
;    print, resolved_routine, format='(A,%"\N")'
;
;endif $
;else print, strupcase(proname), format='("Module ",A," Not Compiled.",%"\N")'

test = 0L
if (presolved or fresolved or punresolved or funresolved) then test=1L

; Not resolved
; maybe could also check for errors in resolve_routine.pro
if test eq 0 then begin

  ; EXTRACT THE !PATH INTO A STRING ARRAY...
  path = strsplit(!path, ':', /EXTRACT)

  ; GET RID OF "." IF USER INCLUDES THIS IN PATH...
  path = path[where(path ne '.')]

  ; SEARCH CURRENT DIRECTORY, EVEN IF NOT IN IDL PATH...
  cd, CURRENT=current
  if (total(strmatch(path,current)) eq 0) then path = [current,path]

  ; ADD THE FILENAME TO EACH PATH DIRECTORY...
  filenames = path + '/' + proname + '.pro'

  ; DOES ANY SUCH FILE EXIST IN THE CURRENT PATH...
  file_exists = where(file_test(filenames), N_exists)

  if n_exists gt 0 then test=1L

end


; Check for embedded files
; This is very time intensive
if test eq 0 and keyword_set(embedded) then begin

  ; EXTRACT THE !PATH INTO A STRING ARRAY...
  path = strsplit(!path, ':', /EXTRACT)

  ; GET RID OF "." IF USER INCLUDES THIS IN PATH...
  path = path[where(path ne '.')]

  ; SEARCH CURRENT DIRECTORY, EVEN IF NOT IN IDL PATH...
  cd, CURRENT=current
  if (total(strmatch(path,current)) eq 0) then path = [current,path]


  npath = n_elements(path)
  N_exists = 0
  i=0
  while (N_exists eq 0 and i lt npath-1) do begin

    cd,path[i]

    ; Check for procedures    
    spawn,'grep -i "pro '+proname+'" *pro',pro_out,pro_errout
    dum = where(stregex(pro_out,proname,/boolean,/fold_case) eq 1,N_exists)

    if N_exists eq 0 then begin
      spawn,'grep -i "function '+proname+'" *pro',func_out,func_errout
      dum = where(stregex(func_out,proname,/boolean,/fold_case) eq 1,N_exists)
    endif

    i++
  end

  if N_exists gt 0 then test=1

endif


; Print status
if test eq 0 and not keyword_set(silent) then $
  print,strupcase(proname)+'.PRO DOES NOT EXIST'




;; IF THERE IS NO SUCH FILE THEN SPLIT...
;if (N_exists eq 0) then begin
;    if (N_elements(resolved_routine) eq 0) then $
;        message, proname + '.pro not found on IDL !path.', /INFO
;    ;return
;    stop
;endif
;
;; PULL OUT ALL THE FILES THAT EXIST...
;filenames = filenames[file_exists]
;
;; TAKE RESOLVED ROUTINE OUT OF THE LIST...
;if (N_elements(resolved_routine) gt 0) then begin
;
;    ; GET THE INDICES OF THE UNRESOLVED ROUTINES...
;    file_exists = where(strmatch(filenames,resolved_routine) eq 0, N_exists)
;
;    ; WAS THE RESOLVED ROUTINE THE ONLY ONE...
;    ;if (N_exists eq 0) then return
;
;    filenames = filenames[file_exists]
;endif
;
; PRINT THE REMAINING ROUTINES...
;print, 'Other Files Containing Module '+strupcase(proname)+' in IDL !path:'
;print, transpose(filenames)
;print

return,test

end
