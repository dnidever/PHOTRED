;+
;
; PHOTRED_TILEGROUPS
;
; This program takes files and creates groups based
; on the tiling scheme.
;
; INPUTS:
;  files      The FITS files to group.
;  tilestr    The structure with information on the tiling scheme.
;  =logfile   The name of the log file.
;
; OUTPUTS:
;  groupstr   The structure with information for each group.
;  =count     The number of groups.
;  =error     The error message if one occurred.
;
; USAGE:
;  IDL>photred_tilegroups,files,tilestr,groupstr
;
; By D.Nidever  Jan 2017
;-

pro photred_tilegroups,files,tilestr,groupstr,count=count,error=error,logfile=logfile

undefine,groupstr,error,count

; Not enough inputs
if n_elements(files) eq 0 or n_elements(tilestr) eq 0 then begin
  print,'Syntax - photred_tilegroups,files,tilestr,groupstr,count=count,error=error,logfile=logfile'
  error = 'Not enough inputs'
  return
endif

if n_elements(logfile) eq 0 then logfile=-1

; Gather information about the files
nfiles = n_elements(files)
PHOTRED_GATHERFILEINFO,files,filestr

;; Figure out which images overlap the tiles
overlapdata = replicate({filenum:-1L,tilenum:-1L},nfiles*tilestr.ntiles)
ocount = 0LL
for i=0,nfiles-1 do begin
  HEAD_ADXY,tilestr.head,filestr[i].vertices_ra,filestr[i].vertices_dec,fx,fy,/deg
  ;; Loop over the tiles
  for j=0,tilestr.ntiles-1 do begin
    tile1 = tilestr.tiles[j]
    tx = [tile1.x0,tile1.x1,tile1.x1,tile1.x0]
    ty = [tile1.y0,tile1.y0,tile1.y1,tile1.y1] 
    overlap = DOPOLYGONSOVERLAP(fx,fy,tx,ty)
    if overlap eq 1 then begin
      overlapdata[ocount].filenum = i
      overlapdata[ocount].tilenum = j
      ocount++
    endif
  endfor  
endfor
; Trim excess overlapdata
overlapdata = overlapdata[0:ocount-1]
noverlapdata = n_elements(overlapdata)

; One "group" for each tile that has overlaps
uitilenum = uniq(overlapdata.tilenum,sort(overlapdata.tilenum))
utilenum = overlapdata[uitilenum].tilenum
ngroups = n_elements(utilenum)

; Group loop
groupstr = replicate({tilenum:0L,tilename:'',nfiles:0L,files:ptr_new()},ngroups)
for i=0,ngroups-1 do begin
  MATCH,overlapdata.tilenum,utilenum[i],ind1,ind2,count=nmatch,/sort
  overlapdata1 = overlapdata[ind1]
  gfiles1 = files[overlapdata1.filenum]
  groupstr[i].tilenum = tilestr.tiles[utilenum[i]].num
  groupstr[i].tilename = tilestr.tiles[utilenum[i]].name
  groupstr[i].nfiles = nmatch
  groupstr[i].files = ptr_new(gfiles1)
endfor
count = ngroups

end
