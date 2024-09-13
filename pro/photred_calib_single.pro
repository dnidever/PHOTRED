;+
;
; PHOTRED_CALIB_SINGLE
;
; This does the calibration (i.e runs photcalib.pro) for photred
; See photcalib_prep.pro and photcalib.pro
;
; INPUTS:
;  longfile            Input filename
;  imagerfile          Filename of imager information
;  transfile           Filename of transformation equations.
;  apcorfile           Filename of apcor corrections.
;  =catformat          Type of output catalog format to use.
;  =telescope          Telescope name.
;  =instrument         Instrument name.
;  =observatory        Observatory nam.
;  /ddo51radialoffset  Apply DDO51 radial offset.
;  =fieldra            Field center RA used for ddo51 radial offset
;  =fielddec           Field center DEC used for ddo51 radial offset
;  /mchusetiles        Use tiles.
;  /keepinstr          Keep instrumental magnitudes.
;  /avgmag             Add average magnitudes.
;  /avgonlymag         Only keep average magnitudes.
;  /redo               Redo files that were already done.
;  /stp                Stop at the end of the program.
;
; OUTPUTS:
;  The calibrated photometry PHOT files from PHOTCALIB
;
; By D.Nidever  Mar 2008
;-

pro photred_calib_single,longfile,imagerfile,transfile,apcorfile,catformat=catformat,$
                         telescope=telescope,instrument=instrument,observatory=observatory,$
                         ddo51radialoffset=ddo51radialoffset,fieldra=fieldra,fielddec=fielddec,$
                         mchusetiles=mchusetiles,keepinstr=keepinstr,avgmag=avgmag,$
                         avgonlymag=avgonlymag,redo=redo,stp=stp

  if n_elements(catformat) eq 0 then catformat='FITS'
  tilesep = '+'

  file = file_basename(longfile)
  filedir = file_dirname(longfile)

  CD,current=curdir

  if n_elements(telescope) eq 0 or n_elements(instrument) eq 0 then begin
    print,'Need telescope and instrument'
  endif

  if strlowcase(telescope) eq 'blanco' then observatory='ctio'
  if strlowcase(telescope) eq 'swope' then observatory='lco'
  if strlowcase(telescope) eq 'magellan' then observatory='lco'
  if strlowcase(telescope) eq 'lbt' then observatory='mgio'

  if keyword_set(ddo51radialoffset) and (n_elements(fieldra) eq 0 or n_elements(fielddec) eq 0) then begin
    print,'Need FIELDRA and FIELDDEC for /ddo51radialoffset'
    return
  endif

  ; LOAD THE "imagers" FILE
  ;----------------------------
  print,'Loading imager information'
  ; The columns need to be: Telescope, Instrument, Naps, separator
  imagers_fieldnames = ['telescope','instrument','observatory','namps','separator']
  imagers_fieldtpes = [7,7,7,3,7]
  imagers = IMPORTASCII(imagerfile,fieldnames=imagers_fieldnames,$
                        fieldtypes=imagers_fieldtypes,comment='#')
  imagers.telescope = strupcase(strtrim(imagers.telescope,2))
  imagers.instrument = strupcase(strtrim(imagers.instrument,2))
  imagers.observatory = strupcase(strtrim(imagers.observatory,2))
  singleind = where(imagers.namps eq 1,nsingle)
  if nsingle gt 0 then imagers[singleind].separator = ''
  if (n_tags(imagers) eq 0) then begin
    print,'NO imagers in '+scriptsdir+'/imagers'
    return
  endif

  ; What IMAGER are we using??
  ;---------------------------
  ind_imager = where(imagers.telescope eq strupcase(telescope) and $
                     imagers.instrument eq strupcase(instrument),nind_imager)
  if nind_imager eq 0 then begin
    print,'TELESCOPE='+telescope+' INSTRUMENT='+instrument+' NOT FOUND in >>imagers<< file'
    return
  endif
  thisimager = imagers[ind_imager[0]]
  ; print out imager info
  print,''
  print,'USING IMAGER:'
  print,'Telescope = '+thisimager.telescope
  print,'Instrument = '+thisimager.instrument
  print,'Namps = '+strtrim(thisimager.namps,2)
  print,"Separator = '"+thisimager.separator+"'"
  print,''

  ; Load the transformation equations
  READ_TRANS,transfile,trans,/silent

  ; Getting the aperture correction structure
  apcor = IMPORTASCII(apcorfile,fieldnames=['name','value'],/noprint)
  ; Remove the 'a.del' endings for the names
  apcor_orig = apcor
  apcor.name = repstr(apcor.name,'a.del','')  ; base names


  ; Are we keeping the INSTRUMENTAL magnitudes 
  if n_elements(keepinstr) eq 0 then keepinstr=0
  ; Are we AVERAGING multiple observations
  if n_elements(avgmag) eq 0 then avgmag=0
  if n_elements(avgonlymag) eq 0 then avgonlymag=0
  ; If BOTH are set then do AVGONLY
  if keyword_set(avgmag) and keyword_set(avgonlymag) then avgmag=0
  ; Printing PHOTCALIB Settings
  if keyword_set(keepinstr) then print,$
    'KEEPING INSTRUMENTAL MAGNITUDES IN OUTPUT FILES'
  if keyword_set(avgmag) then print,$
    'PUTTING INDIVIDUAL MAGNITUDES AND AVERAGE MAGNITUDES IN OUTPUT FILES'
  if keyword_set(avgonlymag) then print,$
    'PUTTING *ONLY* AVERAGE MAGNITUDES IN OUTPUT FILES'


  file = FILE_BASENAME(longfile)
  filedir = FILE_DIRNAME(longfile)
  base = FILE_BASENAME(file,'.ast')
  mchfile = base+'.mch'

  print,'CALIBRATING ',longfile
  print,systime(0)

  CD,filedir

  ; Check that the AST and MCH files exist
  ;---------------------------------------
  asttest = FILE_TEST(longfile)
  if asttest eq 1 then nastlines = FILE_LINES(longfile) else nastlines=0
  mchtest = FILE_TEST(mchfile)
  if mchtest eq 1 then nmchlines = FILE_LINES(mchfile) else nmchlines=0
  if (nmchlines eq 0) or (nastlines eq 0) then begin
    if asttest eq 0 then print,file+' NOT FOUND'
    if asttest eq 1 and nastlines eq 0 then print,file+' HAS 0 LINES'
    if mchtest eq 0 then print,mchfile,' NOT FOUND'
    if mchtest eq 1 and nmchlines eq 0 then print,mchfile,' HAS 0 LINES'
    CD,curdir
    return
  endif

  ; Check that the individual FITS files exist
  ; and appear in the apcor.lst
  ; the 'a.del' endings were already removed
  LOADMCH,mchfile,alsfiles
  nalsfiles = n_elements(alsfiles)
  for j=0,nalsfiles-1 do begin
    ialsfile = alsfiles[j]
    ibase = FILE_BASENAME(ialsfile,'.als')
    ifitsfile = ibase+'.fits'
    if file_test(ifitsfile) eq 0 then ifitsfile=ibase+'.fits.fz'
    ifitstest = FILE_TEST(ifitsfile)
    igd = where(apcor.name eq ibase,nigd)
    ;; TILES, strip off the tile portion at the end
    if nigd eq 0 and keyword_set(mchusetiles) and stregex(ibase,'\'+tilesep+'T',/boolean) eq 1 then begin
      itilebase = (strsplit(ibase,tilesep,/extract))[0]
      igd = where(apcor.name eq itilebase,nigd)
    endif
    if (ifitstest eq 0) or (nigd eq 0) then begin
      if ifitstest eq 0 then print,ifitsfile,' NOT FOUND'
      if nigd eq 0 then print,ialsfile,' NOT in apcor.lst'
      CD,curdir
      return
    endif
  endfor


  ; DDO51 Radial Offset Correction
  ;-------------------------------
  undefine,photfile
  If keyword_set(ddo51radoffset) and (telescope eq 'BLANCO') and (instrument eq 'MOSAIC') then begin

    ; Check that we have DDO51
    alsbases = FILE_BASENAME(alsfiles,'.als')
    alsfitsfiles = alsbases+'.fits'
    bdfits = where(file_test(alsfitsfiles) eq 0,nbdfits)
    if nbdfits gt 0 then alsfitsfiles[bdfits]=alsbases[bdfits]+'.fits.fz'
    filters = PHOTRED_GETFILTER(alsfitsfiles)

    ; We have a DDO51 filter
    gd_ddo51 = where(filters eq 'D',ngd_ddo51)
    if (ngd_ddo51 gt 0) then begin

      ; Load the AST file
      ast = IMPORTASCII(longfile,/header,/noprint)
      tags = tag_names(ast)

      ; Make sure we have RA/DEC in the AST structure
      if TAG_EXIST(ast,'RA') eq 0 then begin
        print,file+' DOES NOT HAVE RA.  CANNOT DO DDO51 Radial Offset Correction'
        CD,curdir
        return
      endif
      if TAG_EXIST(ast,'DEC') eq 0 then begin
        print,file+' DOES NOT HAVE DEC.  CANNOT DO DDO51 Radial Offset Correction'
        CD,curdir
        return
      endif

      print,'Applying DDO51 Radial Offset Correction'

      ; Convert from RA/DEC to X/Y
      ROTSPHCEN,ast.ra,ast.dec,fieldra,fielddec,xi,eta,/gnomic
      ; The MOSAIC camera is oriented so that N is to the right.
      ; So xi -> YB, and eta -> XB.
      yb = xi*3600./0.26   ; convert from deg to pixels
      xb = eta*3600./0.26
      rad = sqrt(xb^2.0 + yb^2.0)
      xb = xb + 4096
      yb = yb + 4096

      ; Calculate the offset
      ; expr = 'P[0]*exp(-0.5*(X-P[1])^2.0/P[2]^2.0)+P[3]+P[4]*X'
      ; 0.0632479 593.665 600.155 0.0186573 -9.44991e-06
      ddo51_radoffset = 0.0632479*exp(-0.5*(rad-593.665)^2.0/(600.155^2.0)) + 0.0186573 - 9.44991e-06*rad
      ddo51_radoffset = -ddo51_radoffset     ; we want to remove the offset by addition

      ; Add to the structure
      ADD_TAG,ast,'RPIX',0.0,ast
      ADD_TAG,ast,'XB',0.0,ast
      ADD_TAG,ast,'YB',0.0,ast
      ADD_TAG,ast,'DDO51_RADOFFSET',0.0,ast
      ast.rpix = rad
      ast.xb = xb
      ast.yb = yb
      ast.ddo51_radoffset = ddo51_radoffset


      ; Now correct the magnitudes
      for j=0,ngd_ddo51-1 do begin
        ; The columns are: ID, X, Y, MAG1, MAG1ERR, MAG2, MAG2ERR, ...
        thismag = gd_ddo51[j]
        thiscol = 3+2*thismag
        ;thiscol = where(tags eq 'MAG'+strtrim(thismag+1,2),nthiscol)

        ; Applying the correction
        ast.(thiscol) = ast.(thiscol) + ddo51_radoffset

        ; bad values are still bad
        bdval = where(ast.(thiscol) gt 90,nbdval)
        if nbdval gt 0 then ast[bdval].(thiscol)=99.9999

      end ; DDO51 filter loop

      ; Write the new file  
      photfile = base+'.temp'
      PRINTSTR,ast,photfile,/silent

    endif else begin
      print,'NO DDO51 Filter.  CANNOT APPLY DDO51 Radial Offset'
      undefine,photfile
    endelse

  Endif  ; DDO51 radial offset


  ; Make the input file with PHOTRED_PHOTCALIB_PREP
  ;------------------------------------------------
  inpfile = base+'.input'
  print,'Making input file: ',inpfile
  PHOTRED_PHOTCALIB_PREP,mchfile,apcor,inpfile,error=error,imager=thisimager,$
                         observatory=observatory,photfile=photfile

  ; Make sure the input file exists
  inptest = FILE_TEST(inpfile)
  if (inptest eq 0) or (n_elements(error) ne 0) then begin
    if (inptest eq 0) then print,inpfile,' NOT FOUND' else $
        print,inpfile,' ERROR'
    CD,curdir
    return
  endif


  ;*****************************************************
  ; DOES THE TRANS FILE HAVE THE MAGNITUDES WE NEED???
  ;*****************************************************
  ; I think PHOTCALIB deals with this

  ; Check info there's CHIP information in the trans eqns.
  ;  CHIP=-1 means there's no info
  chipinfo = 0
  if tag_exist(trans,'CHIP') then begin
    gchip = where(trans.chip ge 1,ngchip)
    if ngchip gt 0 then chipinfo=1
  endif

  ; CHIP-SPECIFIC transformation equations
  if chipinfo eq 1 then begin
    inptransfile = ''
    ext = first_el(strsplit(base,thisimager.separator,/extract),/last)
    chip = long(ext)
    ;ind = where(trans.chip eq chip,nind)
    ind = where(trans.chip eq chip or trans.chip eq -1,nind)      ; keep trans without chip info as well

    ; Nothing for this chip
    if nind eq 0 then begin
      print,'NO transformation information for CHIP=',strtrim(chip,2)
      CD,curdir
      return
    endif
    inptrans = trans[ind]

  ; Global transformation equations
  endif else begin
    undefine,inptrans
    inptransfile = transfile
  endelse

  ; Run PHOTCALIB
  ;---------------
  print,'Calibrating photometry for ',file
  PHOTCALIB,inpfile,inptransfile,inptrans=inptrans,average=avgmag,keepinstrumental=keepinstr,$
            onlyaverage=avgonlymag,catformat=catformat,logfile=logfile,/header


  ; Check that the PHOT file exists
  ;-----------------------------------
  photfile = base+'.phot'
  phottest = FILE_TEST(photfile)
  if phottest eq 1 then nlines = FILE_LINES(photfile) else nlines=0
  if (phottest eq 0) or (nlines eq 0) then begin
    if phottest eq 0 then print,photfile,' NOT FOUND'
    if phottest eq 1 and nlines eq 0 then print,photfile,' HAS 0 LINES'
  endif

  CD,curdir

  if keyword_set(stp) then stop

end
