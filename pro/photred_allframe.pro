;+
;
; PHOTRED_ALLFRAME
;
; This runs allframe for photred
; Can be on pleione or just a single machine depending
; on where this is launched from.
; See daophot_setup.pro for how to make the pleione scripts
;
; INPUTS:
;  =tiletype  The type of tiling scheme to use in the combination.
;               The default is "ORIG'.
;  /redo      Redo files that were already done.
;  /fake      Run for artificial star tests.
;  /stp       Stop at the end of the program.
;
; OUTPUTS:
;  The ALLFRAME/SEXTRACTOR MAG file
;
; By D.Nidever  Feb 2008
;-

pro photred_allframe,redo=redo,fake=fake,nmulti=nmulti,tiletype=tiletype

COMMON photred,setup

print,''
print,'########################'
print,'RUNNING PHOTRED_ALLFRAME'
print,'########################'
print,''

CD,current=curdir

; Does the logs/ directory exist?
testlogs = file_test('logs',/directory)
if testlogs eq 0 then FILE_MKDIR,'logs'


; Log files
;----------
thisprog = 'ALLFRAME'
logfile = 'logs/'+thisprog+'.log'
logfile = FILE_EXPAND_PATH(logfile)  ; want absolute filename
if file_test(logfile) eq 0 then SPAWN,'touch '+logfile,out


; Print date and time to the logfile
printlog,logfile,''
printlog,logfile,'Starting PHOTRED_'+thisprog+'  ',systime(0)


; Check that all of the required programs are available
progs = ['allframe','allfprep','readline','readlist','readpar','photred_getinput','photred_updatelists',$
         'photred_loadsetup','undefine','loadmch','push','pbs_daemon','pbs_checkstat','pbs_makescript',$
         'loadraw','arr2str','mktemp','printlog','first_el','importascii','strsplitter','add_tag','iraf_imshift',$
         'iraf_imalign','iraf_imcombine','mkopt','randomize','touchzero','writecol','writeline','maketemp',$
         'stress','combine_structs','getpixscale','imfwhm','loadals','loadinput','printline','writeals',$
         'mad','rndint','strep','head_xyad','hdr2wcstnx','wcstnx_xy2rd','parsetnx','wcstnxcor','xieta2rd',$
         'fiximages','iraf_run','check_iraf','ia_trim', 'getparam']
test = PROG_TEST(progs)
if min(test) eq 0 then begin
  bd = where(test eq 0,nbd)
  printlog,logfile,'SOME NECESSARY PROGRAMS ARE MISSING:'
  printlog,logfile,progs[bd]
  return
endif


; Check that the DAOPHOT/ALLSTAR/SEXTRACTOR/ALLFRAME programs exist
SPAWN,'which daophot',out,errout
daophotfile = FILE_SEARCH(out,count=ndaophotfile)
if (ndaophotfile eq 0) then begin
  print,'DAOPHOT PROGRAM NOT AVAILABLE'
  return
endif
SPAWN,'which allstar',out,errout
allstarfile = FILE_SEARCH(out,count=nallstarfile)
if (nallstarfile eq 0) then begin
  print,'ALLSTAR PROGRAM NOT AVAILABLE'
  return
endif
SPAWN,'which allframe',out,errout
allframefile = FILE_SEARCH(out,count=nallframefile)
;allframefile = FILE_SEARCH('/net/halo/bin/allframe.2004.fixed',count=nallframefile)
;allframefile = FILE_SEARCH('/net/halo/bin/allframe.2008',count=nallframefile)
if (nallframefile eq 0) then begin
  print,'ALLFRAME PROGRAM NOT AVAILABLE'
  return
endif
SPAWN,'which sex',out,errout
sexfile = FILE_SEARCH(out,count=nsexfile)
if (nsexfile eq 0) then begin
  print,'SEXTRACTOR PROGRAM NOT AVAILABLE'
  return
endif


; LOAD THE SETUP FILE if not passed
;-----------------------------------
; This is a 2xN array.  First colume are the keywords
; and the second column are the values.
; Use READPAR.PRO to read it
if n_elements(setup) eq 0 then begin
  PHOTRED_LOADSETUP,setup,count=count
  if count lt 1 then return
endif


; Are we redoing?
doredo = READPAR(setup,'REDO')
if keyword_set(redo) or (doredo ne '-1' and doredo ne '0') then redo=1

; Telescope, Instrument
telescope = READPAR(setup,'TELESCOPE')
telescope = strupcase(strtrim(telescope,2))
if (telescope eq '0' or telescope eq '' or telescope eq '-1') then begin
  printlog,logfile,'NO TELESCOPE FOUND.  Please add to >>photred.setup<< file'
  return
endif
; INSTRUMENT
instrument = READPAR(setup,'INSTRUMENT')
instrument = strupcase(strtrim(instrument,2))
if (instrument eq '0' or instrument eq '' or instrument eq '-1') then begin
  printlog,logfile,'NO INSTRUMENT FOUND.  Please add to >>photred.setup<< file'
  return
endif

; ALLFRAME Source Detection Program
alfdetprog = READPAR(setup,'ALFDETPROG')
if alfdetprog eq '0' or alfdetprog eq '' or alfdetprog eq '-1' then $
   alfdetprog='sextractor'
alfdetprog = strlowcase(alfdetprog)
if alfdetprog eq 'dao' then alfdetprog='daophot'
if alfdetprog eq 'sex' then alfdetprog='sextractor'
if alfdetprog ne 'sextractor' and alfdetprog ne 'daophot' then alfdetprog='sextractor'

; Hyperthread?
hyperthread = READPAR(setup,'hyperthread')
if hyperthread ne '0' and hyperthread ne '' and hyperthread ne '-1' then hyperthread=1
if strtrim(hyperthread,2) eq '0' then hyperthread=0
if n_elements(hyperthread) eq 0 then hyperthread=0

; Scaling images to be combined?
alfnocmbimscale = READPAR(setup,'alfnocmbimscale')
if strtrim(alfnocmbimscale,2) eq '1' then alfnocmbimscale=1 else undefine,alfnocmbimscale

; Trim combined images to overlapping region
alftrimcomb = READPAR(setup,'alftrimcomb')
if strtrim(alftrimcomb,2) eq '1' then alftrimcomb=1 else undefine,alftrimcomb

; Fields to exclude from ALLFRAME processing
alfexclude = READPAR(setup,'alfexclude')
if alfexclude eq '0' or alfexclude eq '' or alfexclude eq '-1' then undefine,alfexclude

; Use common sources file from the reference image
alfusecmn = READPAR(setup,'alfusecmn')
if alfusecmn eq '0' or alfusecmn eq '' or alfusecmn eq '-1' then undefine,alfusecmn

; Type of combination tiling scheme to use
if n_elements(tiletype) eq 0 then begin
  alftiletype = READPAR(setup,'alftiletype')
  if alftiletype eq '0' or alftiletype eq '' or alftiletype eq '-1' then tiletype='ORIG' else tiletype=alftiletype
endif

; MCHUSETILES
mchusetiles = READPAR(setup,'MCHUSETILES')
if mchusetiles eq '0' or mchusetiles eq '' or mchusetiles eq '-1' then undefine,mchusetiles
if n_elements(mchusetiles) gt 0 then tiletype='TILES'
tilesep = '+'
;tilesep = '.'

; Catalog format to use
catformat = READPAR(setup,'catformat')
if catformat eq '0' or catformat eq '' or catformat eq '-1' then catformat='ASCII'
if catformat ne 'ASCII' and catformat ne 'FITS' then catformat='ASCII'

; Temporary working directory
workdir = READPAR(setup,'WORKDIR',count=nworkdir)
if nworkdir eq 0 then undefine,workdir

; Get the IRAF directory from the setup file
;-------------------------------------------
irafdir = READPAR(setup,'IRAFDIR')
if irafdir eq '0' or irafdir eq '-1' then irafdir = '~/iraf/'
irafdir = FILE_SEARCH(irafdir,/fully_qualify,count=nirafdir)
irafdir = irafdir[0]
if nirafdir eq 0 then begin
  print,'NO IRAF DIRECTORY'
  return
endif

; Get the scripts directory from setup
scriptsdir = READPAR(setup,'SCRIPTSDIR')
if scriptsdir eq '' then begin
  printlog,logfile,'NO SCRIPTS DIRECTORY'
  return
endif
pythonbin = READPAR(setup,'pythonbin')
htcondor = READPAR(setup,'htcondor')

; LOAD THE "imagers" FILE
;----------------------------
printlog,logfile,'Loading imager information'
imagerstest = FILE_TEST(scriptsdir+'/imagers')
if (imagerstest eq 0) then begin
  printlog,logfile,'NO >>imagers<< file in '+scriptsdir+'  PLEASE CREATE ONE!'
  return
endif
; The columns need to be: Telescope, Instrument, Naps, separator
imagers_fieldnames = ['telescope','instrument','observatory','namps','separator']
imagers_fieldtpes = [7,7,7,3,7]
imagers = IMPORTASCII(scriptsdir+'/imagers',fieldnames=imagers_fieldnames,$
                      fieldtypes=imagers_fieldtypes,comment='#')
imagers.telescope = strupcase(strtrim(imagers.telescope,2))
imagers.instrument = strupcase(strtrim(imagers.instrument,2))
imagers.observatory = strupcase(strtrim(imagers.observatory,2))
singleind = where(imagers.namps eq 1,nsingle)
if nsingle gt 0 then imagers[singleind].separator = ''
if (n_tags(imagers) eq 0) then begin
  printlog,logfile,'NO imagers in '+scriptsdir+'/imagers'
  return
endif

; What IMAGER are we using??
;---------------------------
ind_imager = where(imagers.telescope eq telescope and imagers.instrument eq instrument,nind_imager)
if nind_imager eq 0 then begin
  printlog,logfile,'TELESCOPE='+telescope+' INSTRUMENT='+instrument+' NOT FOUND in >>imagers<< file'
  return
endif
thisimager = imagers[ind_imager[0]]
; print out imager info
printlog,logfile,''
printlog,logfile,'USING IMAGER:'
printlog,logfile,'Telescope = '+thisimager.telescope
printlog,logfile,'Instrument = '+thisimager.instrument
printlog,logfile,'Namps = '+strtrim(thisimager.namps,2)
printlog,logfile,"Separator = '"+thisimager.separator+"'"
printlog,logfile,''

; Run CHECK_IRAF.PRO to make sure that you can run IRAF from IDL
;---------------------------------------------------------------
CHECK_IRAF,iraftest,irafdir=irafdir
if iraftest eq 0 then begin
  print,'IRAF TEST FAILED.  EXITING'
  return
endif



;###################
; GETTING INPUTLIST
;###################
; INLIST         MCH files
; OUTLIST        MAG files
; SUCCESSLIST    MCH files

; Get input
;-----------
precursor = 'MATCH'
lists = PHOTRED_GETINPUT(thisprog,precursor,redo=redo,ext='mch')
ninputlines = lists.ninputlines


; Exclude fields from ALLFRAME processing
If ninputlines gt 0 and n_elements(alfexclude) gt 0 then begin
  inputlines = lists.inputlines
  inputfiles = file_basename(inputlines)   

  ; Splitting the comma-delimited list
  fields_exclude = strupcase(strtrim(strsplit(alfexclude,',',/extract),2))
  nfields_exclude = n_elements(fields_exclude)
  printlog,logfile,' '
  printlog,logfile,'Excluding fields from ALLFRAME processing: ',fields_exclude

  ; Check if there are any matches with the input list
  toexcludearr = lonarr(ninputlines)
  for i=0,nfields_exclude-1 do begin
    indmatch = where(stregex(inputfiles,'^'+fields_exclude[i]+'-',/boolean) eq 1,nindmatch)
    if nindmatch gt 0 then toexcludearr[indmatch]=1
  endfor
  excludefilesind = where(toexcludearr eq 1,nexcludefiles,comp=keepfilesind,ncomp=nkeepfiles)
  if nexcludefiles eq 0 then printlog,logfile,'No files to exclude'

  ; Some files to exclude
  if nexcludefiles gt 0 then begin
    printlog,logfile,'Excluding ',strtrim(nexcludefiles,2),' input files from ALLFRAME processing'

    ; ---Adding Excluded files to logs/ASTROM.inlist---
    ;  reading in ASTROM.inlist
    READLIST,'logs/ASTROM.inlist',astinputlines,/exist,/unique,/fully,count=nastinputlines,logfile=logfile,/silent
    ; Add new files to the inputlist
    PUSH,astinputlines,inputlines[excludefilesind]
    nastinputlines = n_elements(astinputlines)
    ; Remove redundant names
    ui = UNIQ(astinputlines,sort(astinputlines))
    ui = ui[sort(ui)]
    astinputlines = astinputlines[ui]
    nastinputlines = n_elements(astinputlines)
    ; Write file
    WRITELINE,'logs/ASTROM.inlist',astinputlines
    printlog,logfile,strtrim(nastinputlines,2),' files moved to logs/ASTROM.inlist'

    ; ---Removing excluded files from input list and logs/ALLFRAME.inlist---
    if nkeepfiles gt 0 then begin
      REMOVE,excludefilesind,inputlines    
      ; Write file
      WRITELINE,'logs/ALLFRAME.inlist',inputlines
      ninputlines = n_elements(inputlines)
      printlog,logfile,strtrim(nexcludefiles,2),' excluded files removed from ALLFRAME.inlist'
      printlog,logfile,strtrim(ninputlines,2),' non-excluded input files left in ALLFRAME.inlist'
    endif else begin
      ; All files removed
      FILE_DELETE,'logs/ALLFRAME.inlist'
      TOUCHZERO,'logs/ALLFRAME.inlist' 
    endelse

    ; ---Reload the input file---
    lists = PHOTRED_GETINPUT(thisprog,precursor,redo=redo,ext='mch')
    ninputlines = lists.ninputlines
    
  endif   ; files to exclude
Endif  ; allframe exclude files


; No files to process
;---------------------
if ninputlines eq 0 then begin
  printlog,logfile,'NO FILES TO PROCESS'
  return
endif

inputlines = lists.inputlines



; Unique directories
inputdirs = FILE_DIRNAME(inputlines)
uidirs = uniq(inputdirs,sort(inputdirs))
uidirs = uidirs[sort(uidirs)]
dirs = inputdirs[uidirs]
ndirs = n_elements(dirs)



; Copy the scripts to the directories
;-------------------------------------
; Checking that the scripts exist
scripts = ['getpsf.sh','photo.opt','apcor.opt','lstfilter.py','goodpsf.pro','allframe.opt',$
           'default.sex','default.param','default.nnw','default.conv']
nscripts = n_elements(scripts)
for i=0,nscripts-1 do begin
  scriptfile = FILE_SEARCH(scriptsdir+'/'+scripts[i],count=nscriptfile)
  if (nscriptfile eq 0) then begin
    print,scriptsdir+'/'+scripts[i],' NOT FOUND'
    return
  endif
endfor
; Copy the scripts to the directories
for i=0,ndirs-1 do begin
  FILE_COPY,scriptsdir+'/'+scripts,dirs[i],/overwrite
endfor




;####################################################
;#  PROCESSING THE FILES
;####################################################
printlog,logfile,''
printlog,logfile,'-----------------------'
printlog,logfile,'PROCESSING THE FILES'
printlog,logfile,'-----------------------'
printlog,logfile,''
printlog,logfile,systime(0)

successarr = intarr(ninputlines)-1         ; 0-bad, 1-good
undefine,outputarr

;--------------------------------
; Checking necessary input files
;--------------------------------
For i=0,ninputlines-1 do begin

  longfile = inputlines[i]
  file = FILE_BASENAME(longfile)
  filedir = FILE_DIRNAME(longfile)

  CD,filedir

  ;CD,mchdirlist[i]
  ;
  ;fil = mchbaselist[i]

  ; Check that the MCH exists
  test = FILE_TEST(file)
  if test eq 1 then nlines=FILE_LINES(file)

  ; NO or BAD MCH file
  if (test eq 0) or (nlines eq 0) then begin
    successarr[i]=0
    if test eq 0 then printlog,logfile,file,' NOT FOUND'
    if test eq 1 and nlines eq 0 then printlog,logfile,file,' HAS 0 LINES'
    goto,BOMB
  endif

  ; Check that the RAW file exists
  mchbase = file_basename(file,'.mch')
  test = FILE_TEST(mchbase+'.raw')
  if test eq 1 then nlines=FILE_LINES(mchbase+'.raw')

  ; NO or BAD RAW file
  if (test eq 0) or (nlines eq 0) then begin
    successarr[i]=0
    if test eq 0 then printlog,logfile,mchbase,'.raw NOT FOUND'
    if test eq 1 and nlines eq 0 then printlog,logfile,mchbase,'.raw HAS 0 LINES'
    goto,BOMB
  endif

  ; Load the MCH file
  LOADMCH,file,files,trans

  ; Check that the individual FITS, ALS, and OPT files exist
  ;----------------------------------------------------------
  nfiles = n_elements(files)
  For j=0,nfiles-1 do begin

    base = FILE_BASENAME(files[j],'.als')
  
    ; Checking FITS file
    fitstest = FILE_TEST(base+'.fits') OR FILE_TEST(base+'.fits.fz')
    if fitstest eq 0 then begin
      printlog,logfile,base+'.fits/.fits.fz NOT FOUND'
      successarr[i] = 0
    endif


    ; Make sure that the FITS files are FLOAT
    ;----------------------------------------
    ; Make sure that BITPIX = -32 otherwise this can cause problems
    if file_test(base+'.fits') eq 0 and file_test(base+'.fits.fz') eq 1 then begin
      fpack = 1
      head = PHOTRED_READFILE(base+'.fits.fz',exten=1,/header)
    endif else begin
      fpack = 0
      head = PHOTRED_READFILE(base+'.fits',/header)
    endelse
    bitpix = long(SXPAR(head,'BITPIX',/silent))
    if (bitpix eq 8 or bitpix eq 16 or bitpix eq 32) and (fpack eq 0) then begin
      printlog,logfile,'BIXPIX = ',strtrim(bitpix,2),'.  Making image FLOAT'

      ; Read in the image
      im = PHOTRED_READFILE(base+'.fits',head,error=error)

      ; Make sure BZERO=0
      bzero = sxpar(head,'BZERO',count=nbzero,/silent)
      if nbzero gt 0 then sxaddpar,head,'BZERO',0.0

      ; Write the FLOAT image
      ; Write the FLOAT image
      if n_elements(error) eq 0 and size(im,/type) lt 4 then $
        FITS_WRITE_RESOURCE,base+'.fits',float(im),head

      ; There was a problem reading the image
      if n_elements(error) gt 0 then begin
        printlog,logfile,'PROBLEM READING IN ',base+'.fits'
        successarr[i] = 0
      endif
    endif

    ; Checking OPT file
    opttest = FILE_TEST(base+'.opt')
    if opttest eq 0 then begin
      printlog,logfile,base+'.opt NOT FOUND'
      successarr[i] = 0
    endif

    ; Checking ALS.OPT file
    alsopttest = FILE_TEST(base+'.als.opt')
    if alsopttest eq 0 then begin
      printlog,logfile,base+'.als.opt NOT FOUND'
      successarr[i] = 0
    endif

    ; Checking AP file
    aptest = FILE_TEST(base+'.ap')
    if aptest eq 0 then begin
      printlog,logfile,base+'.ap NOT FOUND'
      successarr[i] = 0
    endif

    ; Checking ALS file
    alstest = FILE_TEST(base+'.als')
    if alstest eq 0 then begin
      printlog,logfile,base+'.als NOT FOUND'
      successarr[i] = 0
    endif

    ; Checking LOG file
    logtest = FILE_TEST(base+'.log')
    if logtest eq 0 then begin
      printlog,logfile,base+'.log NOT FOUND'
      successarr[i] = 0
    endif

  endfor  ; files in MCH file loop

  ; Is everything okay for this file
  if (successarr[i] ne 0) then successarr[i]=1

  BOMB:

  CD,curdir

endfor  ; loop through MCH files


; Files with okay option files
gd = where(successarr eq 1,ngd)
if (ngd eq 0) then begin
  print,'NO GOOD MCH FILES TO PROCESS'

  ; UPDATE the Lists
  PUSH,failurelist,inputlines
  PHOTRED_UPDATELISTS,lists,failurelist=failurelist,setupdir=curdir,/silent

  return
endif


; These are the files to process
procbaselist = FILE_BASENAME(inputlines[gd])
procdirlist = FILE_DIRNAME(inputlines[gd])
nprocbaselist = n_elements(procbaselist)

; How many FIND iterations do we want
finditer = READPAR(setup,'FINDITER')
finditer = strtrim(finditer,2)
if finditer eq '0' or finditer eq '-1' or finditer eq '' then finditer='2'   ; 2->1 on 7/22/20


; Have some been done already??
;------------------------------
If not keyword_set(redo) then begin

  ; Checking possible output files 
  donearr = lonarr(nprocbaselist)
  for i=0,nprocbaselist-1 do begin
    base = FILE_BASENAME(procbaselist[i],'.mch')
    magfile = procdirlist[i]+'/'+base+'.mag'
    ; Check that this file has a MAG file
    magtest = FILE_TEST(magfile)
    if magtest eq 1 then maglines=FILE_LINES(magfile) else maglines=0
    ; Done properly?
    if (maglines gt 3) then donearr[i]=1
  endfor

  ; Some done already. DO NOT REDO
  bd = where(donearr eq 1,nbd)
  if nbd gt 0 then begin
    ; Print out the names
    printlog,logfile,''
    for i=0,nbd-1 do $
      printlog,logfile,procdirlist[bd[i]]+'/'+procbaselist[bd[i]]+' ALLFRAME ALREADY DONE'

    ; Add these to the "success" and "outlist" list
    PUSH,successlist,procdirlist[bd]+'/'+procbaselist[bd]
    PUSH,outlist,procdirlist[bd]+'/'+procbaselist[bd]+'.als'
    PHOTRED_UPDATELISTS,lists,outlist=outlist,successlist=successlist,$
                        failurelist=failurelist,setupdir=curdir,/silent

    ; Remove them from the arrays
    if nbd lt nprocbaselist then REMOVE,bd,procbaselist,procdirlist
    if nbd eq nprocbaselist then UNDEFINE,procbaselist,procdirlist
    nprocbaselist = n_elements(procbaselist)

    printlog,logfile,''
    printlog,logfile,'REMOVING '+strtrim(nbd,2)+' files from INLIST.  '+$
                     strtrim(nprocbaselist,2)+' files left to PROCESS'
    printlog,logfile,''
  endif

  ; No files to run
  if nprocbaselist eq 0 then begin
    printlog,logfile,'NO FILES TO PROCESS'
    return
  endif

Endif ; some done already?


;-----------------------------
; Running ALLFRAME.PRO
;-----------------------------

printlog,logfile,''
printlog,logfile,'Running ALLFRAME on '+strtrim(n_elements(procbaselist),2)+' files'
printlog,logfile,''
printlog,logfile,systime(0)

;; Use tiles
if tiletype eq 'TILES' then begin
   
  ;; Loop over the fields in this directory                                                                             
  dum = strsplitter(procbaselist,'-',/extract)
  allfields = reform(dum[0,*])
  uifields = uniq(allfields,sort(allfields))
  ufields = allfields[uifields]
  nfields = n_elements(ufields)
  undefine,cmd
  For i=0,nfields-1 do begin
    thisfield = ufields[i]
    gdfieldmch = where(allfields eq thisfield,ngdfieldmch)
    fieldmchfiles = procbaselist[gdfieldmch]

    ;; Load the Field tile file
    ;;  the files for each tile are in F#/F#-TXX/ but the tile file is
    ;;  in directory F#/
    tilefile = procdirlist[gdfieldmch[0]]+'/'+thisfield+'.tiling'
    if file_test(tilefile) eq 0 then tilefile = file_dirname(procdirlist[gdfieldmch[0]])+'/'+thisfield+'.tiling'
    PHOTRED_LOADTILEFILE,tilefile,tilestr

    For j=0,ngdfieldmch-1 do begin
      mchfile = fieldmchfiles[j]  ; F1-00507801_01+T1.mch
      thistile = (strsplit(file_basename(mchfile,'.mch'),tilesep,/extract))[1]
      tilename = thisfield+'-'+thistile
      tind = where(tilestr.tiles.name eq tilename,ntind)
      tstr = tilestr.tiles[tind[0]]
      ; Create the TILE structure
      ;; Offset XREF and YREF for this tile
      tile = {type:'WCS',naxis:long([tilestr.nx,tilestr.ny]),cdelt:double([tilestr.xstep,tilestr.ystep]),$
              crpix:double([tilestr.xref-(tstr.x0-1),tilestr.yref-(tstr.y0-1)]),$
              crval:double([tilestr.cenra,tilestr.cendec]),ctype:['RA---TAN','DEC--TAN'],$
              xrange:[1,tstr.nx],yrange:[1,tstr.ny],$
              nx:tstr.nx,ny:tstr.ny}
      ; Create the string representation of the TILE structure
      stile = "{type:'WCS',naxis:["+strtrim(tilestr.nx,2)+"L,"+strtrim(tilestr.ny,2)+"L],"+$
              "cdelt:["+strdouble(tilestr.xstep)+","+strdouble(tilestr.ystep)+"],"+$
              "crpix:["+strdouble(tilestr.xref-(tstr.x0-1))+","+strdouble(tilestr.yref-(tstr.y0-1))+"],"+$
              "crval:["+strdouble(tilestr.cenra)+","+strdouble(tilestr.cendec)+"],"+$
              "ctype:['RA---TAN','DEC--TAN'],"+$
              "xrange:[1L,"+strtrim(tstr.nx,2)+"L],"+$
              "yrange:[1L,"+strtrim(tstr.ny,2)+"L],"+$
              "nx:"+strtrim(tstr.nx,2)+"L,ny:"+strtrim(tstr.ny,2)+"L}"
      ; Create the string representation of the THISIMAGER structure
      simager = "{telescope:'"+thisimager.telescope+"',instrument:'"+thisimager.instrument+"',"+$
                "observatory:'"+thisimager.observatory+"',namps:"+strtrim(thisimager.namps,2)+","+$
                "separator:'"+thisimager.separator+"'}"
      ; Make commands for allframe
      cmd1 = "allframe,'"+mchfile+"'"+',setupdir="'+curdir+'",scriptsdir="'+scriptsdir+'",irafdir="'+irafdir+'",finditer='+finditer+$
            ",detectprog='"+alfdetprog+"',catformat='"+catformat+"',tile="+stile+",imager="+simager
      if keyword_set(alfnocmbimscale) then cmd1+=",/nocmbimscale"
      if keyword_set(alftrimcomb) then cmd1+=",/trimcomb"
      if keyword_set(alfusecmn) then cmd1+=",/usecmn"
      if keyword_set(fake) then cmd1+=",/fake"
      if n_elements(workdir) gt 0 then cmd1+=",workdir='"+workdir+"'"
      PUSH,cmd,cmd1
    Endfor  ; mch file loop
  Endfor  ; field loop

;; WCS or ORIG
endif else begin
  ; Make commands for allframe
  cmd = "allframe,'"+procbaselist+"'"+',setupdir="'+curdir+'",scriptsdir="'+scriptsdir+'",irafdir="'+irafdir+'",finditer='+finditer+$
        ",detectprog='"+alfdetprog+"',catformat='"+catformat+"',tile={type:'"+strupcase(tiletype)+"'}"
  if keyword_set(alfnocmbimscale) then cmd+=",/nocmbimscale"
  if keyword_set(alftrimcomb) then cmd+=",/trimcomb"
  if keyword_set(alfusecmn) then cmd+=",/usecmn"
  if keyword_set(fake) then cmd+=",/fake"
  if n_elements(workdir) gt 0 then cmd+=",workdir='"+workdir+"'"
endelse

; Getting NMULTI from setup file if not input
if n_elements(nmulti) eq 0 then begin
  nmulti = READPAR(setup,'NMULTI')
  nmulti = long(nmulti)

  ; Use NMULTI_ALLFRAME if set
  nmultiallframe = READPAR(setup,'NMULTI_ALLFRAME')
  if nmultiallframe ne '0' and nmultiallframe ne '' and nmultiallframe ne '-1' then nmulti=long(nmultiallframe)
  nmulti = nmulti > 1  ; must be >=1
endif


; What host
host = getenv('HOST')
pleione = stregex(host,'pleione',/boolean,/fold_case)
hyades = stregex(host,'hyades',/boolean,/fold_case)


; Running on multiple machines
if (((pleione eq 1) or (hyades eq 1)) or keyword_set(hyperthread) or (htcondor ne 0))  then begin

  if htcondor ne 0 then begin
    ; Get SETUP to see if we are using IDL VM and/or we have some submit commands
    htcondor_cmd   = getparam(htcondor_cmd, 'htcondor_cmd',   setup, '', logfile)
    htcondor_idlvm = getparam(htcondor_cmd, 'htcondor_idlvm', setup, '', logfile)

    ; Build the HTCondor submit commands (Use |=| for assignations and |;| to separate commands
    ; Add user's HTCondor commands and needed environment variables (PATH, LD_LIBRARY_PATH and HOME, and those needed by IRAF)
    htcondor_cmd  = "#shared|;|" + htcondor_cmd
    htcondor_cmd += "|;|environment|=|PATH=" + getenv('PATH')            $
                 +       ";LD_LIBRARY_PATH=" + getenv('LD_LIBRARY_PATH') $
                 +                  ";iraf=" + getenv('iraf')            $
                 +              ";IRAFARCH=" + getenv('IRAFARCH')        $
                 +                 ";SHELL=" + getenv('SHELL')           $
                 +                  ";HOME=" + getenv('HOME')            $
                 +           ";IDL_STARTUP=" + getenv('IDL_STARTUP')
  endif else begin
    htcondor_cmd   = '0' 
    htcondor_idlvm = '0' 
  endelse

  ; Submit the jobs to the daemon
  PBS_DAEMON,cmd,procdirlist,/idle,nmulti=nmulti,prefix='alf',hyperthread=hyperthread, $
             pythonbin=pythonbin,scriptsdir=scriptsdir,htcondor=htcondor_cmd,htc_idlvm=htcondor_idlvm,waittime=15,/cdtodir


; Normal, single jobs
endif else begin

  ; Looping through the files
  FOR i=0,nprocbaselist-1 do begin
    cd,procdirlist[i]
    ALLFRAME,procbaselist[i],scriptsdir=scriptsdir,finditer=finditer,detectprog=alfdetprog
    cd,curdir
  ENDFOR

endelse


; IT WOULD BE BETTER TO UPDATE THE LISTS
; AFTER EACH FILE IS PROCESSED!!!!!
; CAN PBS_DAOMON DO THAT?????

;-------------------
; Checking OUTPUTS
;-------------------

; Loop through all files in mchbaselist
for i=0,ninputlines-1 do begin

  longfile = inputlines[i]
  file = FILE_BASENAME(longfile)
  filedir = FILE_DIRNAME(longfile)

  ; Only check "good" files
  if (successarr[i] eq 1) then begin

    CD,filedir

    base = FILE_BASENAME(file,'.mch')
    magfile = base+'.mag'
    ; Check that this file has a MAG file
    magtest = FILE_TEST(magfile)

    ; We have a MAG file
    if (magtest eq 1) then begin
      PUSH,outputarr,filedir+'/'+magfile
    endif else begin
      successarr[i]=0
      printlog,logfile,magfile,' NOT FOUND'
    endelse

  endif

  CD,curdir

endfor



;##########################################
;#  UPDATING LIST FILES
;##########################################
undefine,outlist,successlist,failurelist

; Success List
ind = where(successarr eq 1,nind)
if nind gt 0 then successlist = inputlines[ind] else UNDEFINE,successlist

; Output List
; Creating the new output array, ALS files
noutputarr = n_elements(outputarr)
if (noutputarr gt 0) then begin
  outlist = outputarr
endif else UNDEFINE,outlist

; Failure List
bd = where(successarr eq 0,nbd)
if (nbd gt 0) then begin
  failurelist = inputlines[bd]
endif else UNDEFINE,failurelist

PHOTRED_UPDATELISTS,lists,outlist=outlist,successlist=successlist,$
                    failurelist=failurelist,setupdir=curdir


;;######################
;;  CLEANING UP
;;######################
if keyword_set(clean) then begin
  printlog,logfile,'CLEANING UP.  CLEAN='+strtrim(clean,2)

  ;; Only clean up for successful files
  READLIST,curdir+'/logs/ALLFRAME.success',mchfiles,/unique,/fully,setupdir=curdir,count=nmchfiles,logfile=logfile,/silent
  for i=0,mchfiles-1 do begin
    dir1 = file_dirname(mchlist[i])
    base1 = file_basename(mchlist[i])
    ;; _comb
    ;; lst, lst1, lst2, lst1.chi, grp, nst, lst2.chi, plst.chi, psfini.ap
    ;; nei, als.inp, a.fits, cmn.log, cmn.coo, cmn.ap, cmn.lst,
    ;; _sub.fits, _sub.cat, _sub.als, _all.coo, makemag
    FILE_DELETE,dir1+'/'+base1+'_comb'+['.lst','.lst1','.lst2','.lst1.chi','.lst2.chi','.grp','.nst','.plst.chi',$
                                         '.nei','.als.inp','.cmn.log','.cmn.coo','.cmn.ap','.cmn.lst','a.fits',$
                                         'a.fits.fz','_sub.fits','_sub.cat','_sub.als','_all.coo','.makemag'],/allow
    ;; If CLEAN=2, then also remove coo, ap, cat, s.fits, mask.fits as well
    if clean ge 2 then FILE_DELETE,dir1+'/'+base1+['.coo','.ap','.cat','s.fits','s.fits.fz','.mask.fits'],/allow
  endfor
endif


printlog,logfile,'PHOTRED_ALLFRAME Finished  ',systime(0)

if keyword_set(stp) then stop

end
