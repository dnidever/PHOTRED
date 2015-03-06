pro pbs_makescript,input,dir=dir0,name=name,scriptname=scriptname,$
                    idle=idle,prefix=prefix,hyperthread=hyperthread

;+
;
; PBS_MAKESCRIPT
;
; This makes PBS scripts
;
; INPUTS:
;  input    The command to execute.  Can be an array.  idlbatch will
;             be used if it's an IDL command.  If the command is a
;             series of unix commands separated by commas then these
;             will be put on separate lines.
;  =dir     The directory to put the PBS script in.
;  =name    The name to call the PBS script (without the '.sh' ending)
;  /idle    This is an IDL command, otherwise a SHELL command.
;  =prefix  The prefix to use for the PBS script name. "pr" by default.
;  /hyperthread  Not on a PBS server but one with multiple hyperthreaded
;                  processors.
;
; OUTPUTS:
;  PBS scripts output to the direcotires and with the names
;  specified.
;  =scriptname   The absolute names of the scripts
;
; USAGE:
;  IDL>pbs_makescript,input,dir=dir,name=name,idle=idle,scriptname=scriptname,$
;                     prefix=prefix,hyperthread=hyperthread
;
; By D.Nidever   February 2008
;-

ninput = n_elements(input)
ndir = n_elements(dir0)
nname = n_elements(name)

; Not enough inputs
if ninput eq 0 then begin
  print,'Syntax - pbs_makescript,input,dir=dir,name=name,scriptname=scriptname,'
  print,'                        idle=idle,prefix=prefix,hyperthread=hyperthread'
  return
endif

; Not enough directories input
if ndir gt 0 and ndir ne ninput then begin
  print,'INPUT and DIRECTORIES are of different size'
  return
endif


; Current directory
CD,current=curdir

; No directories input
if ndir gt 0 then dir = dir0
if ndir eq 0 then dir = REPLICATE(curdir,ninput)
if ndir eq 1 then dir = REPLICATE(dir0,ninput)    ; multiple commands in same dir

; Make names
if nname eq 0 then begin
  name = strarr(ninput)
  pre = 'pr'
  if n_elements(prefix) gt 0 then pre=prefix[0]
  for i=0,ninput-1 do name[i]=maketemp(pre)   ; I sped this up
  ;for i=0,ninput-1 do name[i]=file_basename(mktemp(pre))
endif

; Make scriptnames
scriptname = dir+'/'+name+'.sh'


; Script loop
FOR i=0,ninput-1 do begin

  base = name[i]
  sname = dir[i]+'/'+base+'.sh'

  ;------
  ; PBS
  ;------
  If not keyword_set(hyperthread) then begin

    ; IDL command
    if keyword_set(idle) then begin

      ; Making an IDL batch file
      bname = dir[i]+'/'+base+'.batch'
      WRITELINE,bname,input[i]

      ; The execution command
      cmd = 'idl < '+base+'.batch'

    ; SHELL command
    endif else begin

      ; The execution command
      cmd = input[i]

      ; If there are commas in the line then break it up into multiple
      ; lines
      commas = strpos(cmd,';')
      if commas[0] ne -1 then begin
        temp = strsplit(cmd,';',/extract)
        cmd = temp
      endif

    endelse


    ; Make the command
    ;----------------------
    undefine,lines
    push,lines,'#!/bin/sh'
    push,lines,'#PBS -l nodes=1:ppn=1'
    push,lines,'#PBS -l walltime=96:00:00'
    push,lines,'#PBS -o '+base+'.report.out'
    push,lines,'#PBS -e '+base+'.error.out'
    ;push,lines,'##PBS -m abe'
    push,lines,'#PBS -M dln5q@virginia.edu'
    push,lines,'#PBS -V'
    push,lines,''
    push,lines,'echo Running on host `hostname`'
    push,lines,'echo Time is `date`'
    push,lines,'echo "Nodes used for this job:"'
    push,lines,'echo "------------------------"'
    push,lines,'cat $PBS_NODEFILE'
    push,lines,'echo "------------------------"'
    push,lines,''
    push,lines,'cd '+dir[i]
    for j=0,n_elements(cmd)-1 do push,lines,cmd[j]
    ;push,lines,cmd
    push,lines,''
    push,lines,'# print end time'
    push,lines,'echo'
    push,lines,'echo "Job Ended at `date`"'
    push,lines,'echo'

    ; Writing the file
    WRITELINE,scriptname[i],lines
    ;print,'Made PBS script. To run: >>qsub '+base+'.sh<<'

    ; Print info
    ;print,'Command: ',cmd
    print,'PBS script written to: ',scriptname[i]

  ;----------------
  ; Hyperthreaded
  ;----------------
  Endif else begin


    ; Just make batch file
    ; treat shell and idl the same

    ; The execution command
    cmd = input[i]
   
    ; If there are commas in the line then break it up into multiple
    ; lines
    commas = strpos(cmd,';')
    if commas[0] ne -1 then begin
      temp = strsplit(cmd,';',/extract)
      cmd = temp
    endif

    ; IDL files should end in .batch
    if keyword_set(idle) then scriptname[i]=dir[i]+'/'+base+'.batch'

    ; Writing the file
    WRITELINE,scriptname[i],cmd

    ; Make SHELL scripts executable
    if not keyword_set(idle) then FILE_CHMOD,scriptname[i],'755'o

    ; Print info
    print,'HYPERTHREAD script written to: ',scriptname[i]

  Endelse


ENDFOR

; Erase the temporary files that MAKETEMP makes
file_delete,name,'\/allow

; Go back to original directory
CD,curdir

;stop

end
