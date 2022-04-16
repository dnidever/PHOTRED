#!/usr/bin/env python

import os
import time
import numpy as np


def combine(filename,tile=None,scriptsdir=None,logfile=None,irafdir=None,
            satlevel=6e4,nocmbimscale=False,fake=False,usecmn=False,imager=None):
    """
    This combines/stacks images for ALLFRAME. 
 
    You need to have run daophot, allstar, daomatch and daomaster 
    already.  There needs to be a fits, opt, als.opt, ap and als 
    file for each file in the MCH file.  There also needs to be an 
    associated RAW file for the MCH file. 
 
    Parameters
    ----------
    filename : str
       The MCH filename
    tile : dict
       Tile information on the sky.
    nocmbimscale : boolean, optional
       Don't scale the images when combining them.  Not 
        recommended, but the old way of doing it.  Bright 
        stars can be missed this way.  Default is False.
    scriptsdir : str
       The directory that contains all of the necessary scripts. 
    irafdir : str
       The IRAF home directory. 
    satlevel : float, optional
       Saturation level.  Default is 6e4.
    logfile : str
       A logfile to print to output to. 
    fake : boolean, optiona
       Run for artificial star tests.  Default is False.
    usecmn : boolean, optional
       Use the individual cmn.lst files to construct a 
         cmn.lst file for the combined image. 
    imager : dict
       Imager structure with basic information. 
 
    Returns
    -------
    The combined image and mask: 
    FILEBASE_comb.fits       combined image 
    FILEBASE_comb.bpm.fits   bad pixel mask 
    FILEBASE_comb.mask.fits  SExtractor weight map (very similar to bpm) 
    FILEBASE_comb.mch        The transformations of the individual frames 
                               to the combined frame. 
    maskdatalevel : float
       The "bad" data level above the highest "good" value 
         in the combined image. 
    filestr : table
       Information on all of the files used for the 
         stacked iamge. 
 
    Example
    -------

    combine('ccd1001.mch',scriptsdir='/net/home/dln5q/daophot/')
 
 
    By D.Nidever   February 2008 
    Automation of steps and scripts by J.Ostheimer and Rachael Beaton 
    Major upgrade of the resampling and combination code, Oct 2016  
    Translated to Python by D. Nidever,  April 2022
    """ 

    global photred,setup 

    # Logfile 
    if keyword_set(logfile): 
        logf=logfile 
    else: 
        logf=-1 
     
    # Getting scripts directory and iraf directory 
    nsetup = len(setup) 
    if nsetup > 0: 
        scriptsdir = setup['SCRIPTSDIR']
        irafdir = setup['IRAFDIR']
     
    # No irafdir 
    if irafdir is None:
        raise ValueError('IRAFDIR NOT INPUT')
    # No irafdir 
    if scriptsdir is None:
        raise ValueError('SCRIPTSDIR NOT INPUT')    
     
    # Run CHECK_IRAF.PRO to make sure that you can run IRAF from IDL 
    #--------------------------------------------------------------- 
    if check_iraf(iraftest,irafdir=irafdir)==False:
        raise ValueError('IRAF TEST FAILED.  EXITING')

    # Check if the scripts exist in the current directory 
    scripts = ['getpsf.sh','photo.opt','apcor.opt','lstfilter.py','goodpsf.pro','allframe.opt',
               'default.sex','default.param','default.nnw','default.conv'] 
    nscripts = len(scripts) 
    # Loop through the scripts 
    for i in range(nscripts): 
        exists = os.path.exists(scriptsdir+'/'+scripts[i])
        if exits:
            size = os.path.size(scriptsdir+'/'+scripts[i])
        else:
            size = 0
        curexists = os.path.exists(scripts[i])
        if curexists:
            cursize = os.path.size(scripts[i])
        else:
            cursize = 0
         
        # No file 
        if exists == False or size == 0: 
            raise ValueError(scriptsdir+'/'+scripts[i],' NOT FOUND or EMPTY')
         
        # Check if the two files are the same size, if not copy it 
        if isize != cursize:
            if os.path.exists(scripts[i]): os.remove(scripts[i])
            shutil.copyfile(scriptsdir+'/'+scripts[i],scripts[i])
     
    # FILENAME 
    mchfile = os.path.basename(filename) 
    mchdir = os.path.dirname(filename) 
    mchbase = os.path.splitext(filename)[0]
     
    # CD to the directory
    curdir = os.path.abspath(os.getcwd())
    os.chdir(mchdir)



    
    # Check that the mch, als, and opt files exist 
    mchtest = os.path.exists(mchfile) 
    if mchtest == 0: 
        printlog,logf,mchfile,' NOT FOUND' 
        return 
     
    # Checking RAW file 
    rawtest = os.path.exists(mchbase+'.raw') 
    if rawtest == 0: 
        printlog,logf,mchbase+'.raw NOT FOUND' 
        return 
     
     
    ############################################ 
    # CHECK NECESSARY FILES 
     
    # Load the MCH file 
    LOADMCH,mchfile,files,trans,magoff 
    nfiles = len(files) 
     
    # FAKE, check that we have all the files that we need 
    if keyword_set(fake): 
        # weights, scale, zero, comb_psf, _shift.mch 
        chkfiles = mchbase+['.weights','.scale','.zero','_comb.psf','_comb.mch'] 
        bdfiles , = np.where(os.path.exists(chkfiles) == 0,nbdfiles) 
        if nbdfiles > 0: 
            error = 'FAKE.  Some necessary files not found. '+strjoin(chkfiles[bdfiles],' ') 
            printlog,logf,error 
            return 
     
     
    # Final files 
    base = os.path.basename(files,'.als') 
    fitsfiles = base+'.fits' 
    outfiles = base+'.shft.fits' 
    outmaskfiles = base+'.mask.shft.fits' 
     
     
    # Gather information on all of the files 
    # photred_gatherfileinfo.pro can do most of this 
    printlog,logf,'Gathering file information' 
    ntrans = len(trans[0,:]) 
    filestr = replicate({fitsfile:'',catfile:'',nx:0L,ny:0L,trans:dblarr(ntrans),magoff:fltarr(2),head:ptr_new(),                      vertices_ra:dblarr(4),vertices_dec:dblarr(4),pixscale:0.0,saturate:0.0,                      background:0.0,comb_zero:0.0,comb_scale:0.0,comb_weights:0.0,                      resampfile:'',resampmask:'',resamptrans:dblarr(ntrans),resamptransrms:0.0},nfiles) 
    filestr.fitsfile = fitsfiles 
    filestr.catfile = files 
    filestr.trans = transpose(trans) 
    filestr.magoff = magoff 
    filestr.resampfile = outfiles 
    filestr.resampmask = outmaskfiles 
    for i in range(nfiles): 
        im1 = PHOTRED_READFILE(filestr[i].fitsfile,head1) 
        filestr[i].head = ptr_new(head1) 
        filestr[i].nx = sxpar(head1,'NAXIS1') 
        filestr[i].ny = sxpar(head1,'NAXIS2') 
        HEAD_XYAD,head1,[0,filestr[i].nx-1,filestr[i].nx-1,0],[0,0,filestr[i].ny-1,filestr[i].ny-1],vra,vdec,/degree 
        filestr[i].vertices_ra = vra 
        filestr[i].vertices_dec = vdec 
        GETPIXSCALE,'',pixscale,head=head1 
        filestr[i].pixscale = pixscale 
        saturate = sxpar(head1,'SATURATE',count=nsaturate,/silent) 
        if nsaturate == 0 : 
            saturate=50000L 
        filestr[i].saturate = saturate 
        gdpix , = np.where(im1 < saturate,ngdpix,ncomp=nbdpix) 
        background = np.median(im1[gdpix]) 
        filestr[i].background = background 
     
     
    ################################################## 
    # Create default reference frame if TILE not input 
    if len(tileinp) > 0: 
        tile=tileinp 
    else: 
        tile={type:'WCS'} 
    if tile.type == 'WCS' and n_tags(tile) == 1: 
        printlog,logf,'Creating TILE projection' 
        # The default projection is a tangent plane centered 
        # at halfway between the ra/dec min/max of all of 
        # the images.  The mean pixel scale is used. 
        #  near RA=0 line 
        if range(filestr.vertices_ra) > 180: 
            vertices_ra = filestr.vertices_ra 
            over , = np.where(vertices_ra > 180,nover,comp=under,ncomp=nunder) 
            if nover > 0 : 
                vertices_ra[over]-=360 
            rar = minmax(vertices_ra) 
            cenra = np.mean(rar) 
        else: 
            rar = minmax(filestr.vertices_ra) 
            cenra = np.mean(rar) 
        decr = minmax(filestr.vertices_dec) 
        cendec = np.mean(decr) 
        pixscale = np.mean(filestr.pixscale) 
        # Set up the tangent plane projection 
        step = pixscale/3600.0d0 
        delta_dec = range(decr) 
        delta_ra = range(rar)*cos(cendec/!radeg) 
        nx = ceil(delta_ra*1.01/step) 
        ny = ceil(delta_dec*1.01/step) 
        xref = nx/2 
        yref = ny/2 
        MKHDR,tilehead,fltarr(5,5) 
        SXADDPAR,tilehead,'NAXIS1',nx 
        SXADDPAR,tilehead,'CDELT1',step 
        SXADDPAR,tilehead,'CRPIX1',xref+1L 
        SXADDPAR,tilehead,'CRVAL1',cenra 
        SXADDPAR,tilehead,'CTYPE1','RA---TAN' 
        SXADDPAR,tilehead,'NAXIS2',ny 
        SXADDPAR,tilehead,'CDELT2',step 
        SXADDPAR,tilehead,'CRPIX2',yref+1L 
        SXADDPAR,tilehead,'CRVAL2',cendec 
        SXADDPAR,tilehead,'CTYPE2','DEC--TAN' 
        EXTAST,tilehead,tileast 
        tileast.equinox = 2000 
         
        printlog,logf,'RA range = [',str(rar[0],2),',',str(rar[1],2),'] deg' 
        printlog,logf,'DEC range = [',str(decr[0],2),',',str(decr[1],2),'] deg' 
        printlog,logf,'Central RA = ',str(cenra,2) 
        printlog,logf,'Central DEC = ',str(cendec,2) 
        printlog,logf,'NX = ',str(nx,2) 
        printlog,logf,'NY = ',str(ny,2) 
         
        # Create the TILE structure 
        tile = {type:'WCS',naxis:int([nx,ny]),cdelt:double([step,step]),crpix:double([xref+1L,yref+1L]),          crval:double([cenra,cendec]),ctype:['RA--TAN','DEC--TAN'],          head:tilehead,ast:tileast,xrange:[0,nx-1],yrange:[0,ny-1],nx:nx,ny:ny} 
     
    # Check that the TILE is valid 
    if allframe_validtile(tile,error=tilerror) == 0: 
        error = tilerror 
        if not keyword_set(silent) : 
            printlog,logf,error 
        return 
     
    # Add HEAD/AST to TILE if needed 
    if tag_exist(tile,'HEAD') == 0: 
        MKHDR,tilehead,fltarr(5,5) 
        SXADDPAR,tilehead,'NAXIS1',tile.naxis[0] 
        SXADDPAR,tilehead,'CDELT1',tile.cdelt[0] 
        SXADDPAR,tilehead,'CRPIX1',tile.crpix[0] 
        SXADDPAR,tilehead,'CRVAL1',tile.crval[0] 
        SXADDPAR,tilehead,'CTYPE1',tile.ctype[0] 
        SXADDPAR,tilehead,'NAXIS2',tile.naxis[1] 
        SXADDPAR,tilehead,'CDELT2',tile.cdelt[1] 
        SXADDPAR,tilehead,'CRPIX2',tile.crpix[1] 
        SXADDPAR,tilehead,'CRVAL2',tile.crval[1] 
        SXADDPAR,tilehead,'CTYPE2',tile.ctype[1] 
        if tag_exist(tile,'CD'): 
            SXADDPAR,tilehead,'CD1_1',tile.cd[0,0] 
            SXADDPAR,tilehead,'CD1_2',tile.cd[0,1] 
            SXADDPAR,tilehead,'CD2_1',tile.cd[1,0] 
            SXADDPAR,tilehead,'CD2_2',tile.cd[1,1] 
        EXTAST,tilehead,tileast 
        tileast.equinox = 2000 
        tile = CREATE_STRUCT(tile,'HEAD',tilehead,'AST',tileast) 
    if tag_exist(tile,'AST') == 0: 
        EXTAST,tile.head,tileast 
        tileast.equinox = 2000 
        tile = CREATE_STRUCT(tile,'AST',ast) 
    # Add XRANGE/YRANGE 
    if tag_exist(tile,'XRANGE') == 0: 
        tile = CREATE_STRUCT(tile,'XRANGE',int([0,tile.naxis[0]-1]),'NX',tile.naxis[0]) 
    if tag_exist(tile,'YRANGE') == 0: 
        tile = CREATE_STRUCT(tile,'YRANGE',int([0,tile.naxis[1]-1]),'NY',tile.naxis[1]) 
     
     
    ############################################ 
    # STEP 1: IMALIGN PREP 
     
    printlog,logf,'------------------------' 
    printlog,logf,'STEP 1: GETTING WEIGHTS' 
    printlog,logf,'------------------------' 
     
    #----------------------------------- 
    # Computs Weights 
    if not keyword_set(fake): 
        ALLFRAME_GETWEIGHTS,mchfile,weights,scales,sky,imager=imager,logfile=logf#,raw2=raw2 
        invscales = 1.0/scales 
        bdscale , = np.where(scales < 1e-5 or invscales > 900,nbdscale) 
        if nbdscale > 0: 
            scales[bdscale] = 1.0 
            invscales[bdscale] = 1.0 
            weights[bdscale] = 0.0 
        weightfile = mchbase+'.weights' 
        WRITECOL,weightfile,weights,fmt='(F10.6)' 
        scalefile = mchbase+'.scale' 
        WRITECOL,scalefile,invscales,fmt='(F10.6)'# want to scale it UP 
        zerofile = mchbase+'.zero' 
        WRITECOL,zerofile,-sky,fmt='(F12.4)'# want to remove the background, set to 1st frame 
         
        # FAKE, use existing ones 
    else: 
        weightfile = mchbase+'.weights' 
        scalefile = mchbase+'.scale' 
        zerofile = mchbase+'.zero' 
        READCOL,weightfile,weights,format='F',/silent 
        READCOL,scalefile,invscales,format='F',/silent 
        scales = 1.0/invscales 
        READCOL,zerofile,sky,format='F',/silent 
        sky = -sky 
    # Stuff the information into the FILEINFO structure 
    filestr.comb_weights = weights 
    filestr.comb_scale = invscales 
    filestr.comb_zero = -sky 
     
     
    ############################################ 
    # STEP 2: Resample the images 
     
    printlog,logf,'---------------------------' 
    printlog,logf,'STEP 2: Resample the images' 
    printlog,logf,'---------------------------' 
     
    CASE tile.type of 
         
        # --- WCS projection on the sky --- 
        'WCS': begin 
         
        # Convert x/y and ra/dec grid for entire tile image 
        #  that we'll be using, ~15sec 
        xb = (lindgen(tile.nx)+tile.xrange[0])#replicate(1,tile.ny) 
        yb = replicate(1,tile.nx)#(lindgen(tile.ny)+tile.yrange[0]) 
        HEAD_XYAD,tile.head,xb,yb,rab,decb,/deg 
         
        # Loop through the files 
        for i in range(nfiles): 
            im1 = PHOTRED_READFILE(filestr[i].fitsfile,head1) 
             
            # Make the mask 
            mask = bytarr(filestr[i].nx,filestr[i].ny)+1 
            gdpix , = np.where(im1 < filestr[i].saturate,ngdpix,comp=bdpix,ncomp=nbdpix) 
            # 0-bad, 1-good 
            if nbdpix > 0: 
                mask[bdpix] = 0 
                im1[bdpix] = filestr[i].background 
             
            # Get X/Y range for this image in the final coordinate system 
            HEAD_ADXY,tile.head,filestr[i].vertices_ra,filestr[i].vertices_dec,vx,vy,/deg 
            xout = [floor(min(vx))-2 > tile.xrange[0], ceil(max(vx))+2 < tile.xrange[1]] 
            xoutrel = xout-tile.xrange[0]# relative to xrange[0] 
            nxout = xout[1]-xout[0]+1 
            yout = [floor(min(vy))-2 > tile.yrange[0], ceil(max(vy))+2 < tile.yrange[1]] 
            youtrel = yout-tile.yrange[0]# relative to yrange[0] 
            nyout = yout[1]-yout[0]+1 
            rr = rab[xoutrel[0]:xoutrel[1],youtrel[0]:youtrel[1]] 
            dd = decb[xoutrel[0]:xoutrel[1],youtrel[0]:youtrel[1]] 
            ALLFRAME_ADXYINTERP,head1,rr,dd,xx,yy,nstep=10 
             
            # The x/y position to bilinear need to be in the original system, ~1sec 
            rim = BILINEAR(im1,xx,yy,missing=filestr[i].background) 
            rmask = BILINEAR(mask,xx,yy,missing=0) 
             
            # Contruct final image 
            fim = fltarr(tile.nx,tile.ny)+filestr[i].saturate 
            fim[xoutrel[0]:xoutrel[1],youtrel[0]:youtrel[1]] = rim 
            fmask = bytarr(tile.nx,tile.ny) 
            fmask[xoutrel[0]:xoutrel[1],youtrel[0]:youtrel[1]] = rmask 
             
            # Contruct the final header 
            fhead = head1 
            # Delete any previous WCS keywords 
            sxdelpar,fhead,['CRVAL1','CRVAL2','CRPIX1','CRPIX2','CDELT1','CDELT2','CTYPE1','CTYPE2'] 
            sxdelpar,fhead,['CD1_1','CD1_2','CD2_1','CD2_2'] 
            pvind , = np.where(stregex(strmid(fhead,0,5),'PV[0-9]_+[0-9]',/boolean) == 1,npvind) 
            if npvind > 0 : 
                remove,pvind,fhead 
            # Add the new WCS 
            PUTAST,fhead,tile.ast 
            sxaddpar,fhead,'NAXIS1',tile.nx 
            sxaddpar,fhead,'NAXIS2',tile.ny 
            sxaddpar,fhead,'BPM',filestr[i].resampmask 
            printlog,logf,filestr[i].resampfile,' ['+str(xout[0],2)+':'+str(xout[1],2)+','+                  str(yout[0],2)+':'+str(yout[1],2)+']' 
            MWRFITS,fim,filestr[i].resampfile,fhead,/create 
            mhead = fhead 
            sxaddpar,mhead,'BITPIX',8 
            sxdelpar,mhead,'BPM' 
            MWRFITS,fmask,filestr[i].resampmask,mhead,/create 
            # this takes about ~37-50 sec for a 2kx4k image. 
            #  now it takes ~5 sec for a 2kx4k image 
 
     
    # --- Pixel based --- 
    'PIXEL': begin 
     
    # Expand the images to sizes that will allow all of the shifts 
    hd1 = PHOTRED_READFILE(fitsfiles[0],/header) 
    nx = sxpar(hd1,'NAXIS1') 
    ny = sxpar(hd1,'NAXIS2') 
    #  left, down, right, up 
    pix_expand = [ abs(floor(min(xshift)) < 0), abs(floor(min(yshift)) < 0),                  ceil(max(xshift)) > 0 , ceil(max(yshift)) > 0 ] 
    outmaskfiles = os.path.dirname(tempfits)+'/'+os.path.basename(tempfits,'.fits')+'.mask.shft.fits' 
    printlog,logf,'Expanding images by [',strjoin(str(pix_expand,2),','),'] pixels' 
    LOADMCH,shiftmch+'.mch',files2,trans2 
    nxf = nx+pix_expand[0]+pix_expand[2] 
    nyf = ny+pix_expand[1]+pix_expand[3] 
    xx1 = lindgen(nx)#replicate(1,ny) 
    yy1 = replicate(1,nx)#lindgen(ny) 
    for i in range(nfiles): 
        # Image 
        tim = PHOTRED_READFILE(tempfits[i],thead) 
        background = np.median(tim) 
        out = trans_coo(xx1[:],yy1[:],reform(trans[i,:])) 
        xx2 = xx1*0. 
        yy2 = yy1*0. 
        xx2[:] = reform(out[0,:]) + pix_expand[0]# shift to expanded grid 
        yy2[:] = reform(out[1,:]) + pix_expand[1]# shift to expanded grid 
        triangulate,xx2,yy2,tr,b# triangulate 
        xout = lindgen(nxf) 
        yout = lindgen(nyf) 
        tim2 = TRIGRID(xx2,yy2,tim, tr, XOUT = xout, YOUT = yout, missing=background) 
        #tim2 = fltarr(nxf,nyf)+background  ; set out of bounds pixels to background 
        #tim2[pix_expand[0]:pix_expand[0]+nx-1,pix_expand[1]:pix_expand[1]+ny-1]=tim 
        thead2 = thead 
        sxaddpar,thead2,'NAXIS1',nxf 
        sxaddpar,thead2,'NAXIS2',nyf 
        sxaddpar,thead2,'CRPIX1',sxpar(thead2,'CRPIX1')+pix_expand[0] 
        sxaddpar,thead2,'CRPIX2',sxpar(thead2,'CRPIX2')+pix_expand[1] 
        MWRFITS,tim2,outfiles[i],thead2,/create 
        # Mask 
        mfile = os.path.basename(tempfits[i],'.fits')+'.mask.fits' 
        mim = PHOTRED_READFILE(mfile,mhead) 
        mim2 = TRIGRID(xx2,yy2,mim, tr, XOUT = xout, YOUT = yout, missing=background) 
        #mim2 = fltarr(nxf,nyf)   ; out of bounds pixels set to 0=bad 
        #mim2[pix_expand[0]:pix_expand[0]+nx-1,pix_expand[1]:pix_expand[1]+ny-1]=mim 
        mhead2 = mhead 
        sxaddpar,mhead2,'NAXIS1',nxf 
        sxaddpar,mhead2,'NAXIS2',nyf 
        sxaddpar,mhead2,'CRPIX1',sxpar(mhead2,'CRPIX1')+pix_expand[0] 
        sxaddpar,mhead2,'CRPIX2',sxpar(mhead2,'CRPIX2')+pix_expand[1] 
        #outmaskfile = FILE_DIRNAME(outfiles[i])+'/'+FILE_BASENAME(outfiles[i],'.fits')+'.mask.shft.fits' 
        MWRFITS,mim2,outmaskfiles[i],mhead2,/create 
 
else: import pdb; pdb.set_trace(),tile.type+' not implemented yet' 
 
 
 
# Creating new MCH file for the combined file 
if not keyword_set(fake): 
print('Deriving new transformation equations for the resampled coordinate system' 
for i in range(nfiles): 
     
    # Convert X/Y of this system into the combined reference frame 
    #  The pixel values are 1-indexed like DAOPHOT uses. 
    ngridbin = 50 
    nxgrid = filestr[i].nx / ngridbin 
    nygrid = filestr[i].ny / ngridbin 
    xgrid = (lindgen(nxgrid)*ngridbin+1)#replicate(1,nygrid) 
    ygrid = replicate(1,nxgrid)#(lindgen(nygrid)*ngridbin+1) 
    HEAD_XYAD,(*filestr[i].head),xgrid-1,ygrid-1,ragrid,decgrid,/deg 
    HEAD_ADXY,tile.head,ragrid,decgrid,refxgrid,refygrid,/deg 
    refxgrid += 1# convert 0-indexed to 1-indexed 
    refygrid += 1 
     
    # Now fit the transformation 
    xdiff = refxgrid-xgrid 
    ydiff = refygrid-ygrid 
    xmed = np.median([xdiff],/even) 
    ymed = np.median([ydiff],/even) 
    # Fit rotation with linear fits if enough points 
    coef1 = robust_poly_fitq(ygrid,xdiff,1)# fit rotation term 
    coef1b = dln_poly_fit(ygrid,xdiff,1,measure_errors=xdiff*0+0.1,sigma=coef1err,/bootstrap) 
    coef2 = robust_poly_fitq(xgrid,ydiff,1)# fit rotation term 
    coef2b = dln_poly_fit(xgrid,ydiff,1,measure_errors=ydiff*0+0.1,sigma=coef2err,/bootstrap) 
    #theta = mean([-coef1[1],coef2[1]]) 
    WMEANERR,[-coef1[1],coef2[1]],[coef1err[1],coef2err[1]],theta,thetaerr 
     
    # [xoff, yoff, cos(th), sin(th), -sin(th), cos(th)] 
    #trans = [xmed, ymed, 1.0, 0.0, 0.0, 1.0] 
    trans = [xmed, ymed, 1.0-theta**2, theta, -theta, 1.0-theta**2] 
    # Adjust Xoff, Yoff with this transformation 
    xyout = trans_coo(xgrid,ygrid,trans) 
    trans[0] += np.median([refxgrid - xyout[0,:]],/even) 
    trans[1] += np.median([refygrid - xyout[1,:]],/even) 
     
    # Fit full six parameters if there are enough stars 
    fa = {x1:(refxgrid)(*),y1:(refygrid)(*),x2:(xgrid)(*),y2:(ygrid)(*)} 
    initpar = trans 
    fpar = MPFIT('trans_coo_dev',initpar,functargs=fa, perror=perror, niter=iter, status=status,                 bestnorm=chisq,:f=dof, autoderivative=1, /quiet) 
    trans = fpar 
     
    diff = trans_coo_dev(fpar,x1=refxgrid,y1=refygrid,x2=xgrid,y2=ygrid) 
    rms = sqrt(np.mean(diff**2.)) 
    filestr[i].resamptrans = trans 
    filestr[i].resamptransrms = rms 
     
    # The output is: 
    # filename, xshift, yshift, 4 trans, mag offset, magoff sigma 
    format = '(A2,A-30,A1,2A10,4A12,F9.3,F8.4)' 
    # In daomaster.f the translations are 10 digits with at most 4 
    # decimal places (with a leading space), the transformation 
    # coefficients are 12 digits with at most 9 decimal places. 
    # Need a leading space to separate the numbers. 
    strans = ' '+[str(string(trans[0:1],format='(F30.4)'),2),                 str(string(trans[2:5],format='(F30.9)'),2)] 
    newline = STRING("'",filestr[i].catfile,"'", strans, filestr[i].magoff[0], rms, format=format) 
    PUSH,mchfinal,newline 
     
    # Printing the transformation 
    printlog,logf,format='(A-20,2A10,4A12,F9.3,F8.4)',filestr[i].catfile,strans,filestr[i].magoff[0],rms 
# Write to the new MCH file 
combmch = mchbase+'_comb.mch' 
WRITELINE,combmch,mchfinal 
 
# FAKE, use existing one 
else: 
combmch = mchbase+'_comb.mch' 
# don't need to load the information 
 
 
############################################ 
# STEP 5: COMBINE IMAGES 
printlog,logf,'-------------------' 
printlog,logf,'STEP 5: IMCOMBINE' 
printlog,logf,'-------------------' 
 
# The imcombine input file 
resampfile = mchbase+'.resamp' 
WRITELINE,resampfile,filestr.resampfile 
 
# SCALE the images for combining 
#------------------------------- 
if not keyword_set(nocmbimscale): 
 
# Put BPM mask names in file headers 
#  these will be used by IMCOMBINE 
#for i=0,nfiles-1 do begin 
#  head = headfits(outfiles[i]) 
#  sxaddpar,head,'BPM',outmaskfiles[i] 
#  modfits,outfiles[i],0,head 
#endfor 
 
# Combine the frames WITH scaling/offset/masking, for the bright stars 
#printlog,logf,'Creating SCALED image' 
combfile = mchbase+'_comb.fits' 
os.remove(combfile,/allow 
os.remove(mchbase+'_comb.bpm.pl',/allow 
IRAF_IMCOMBINE,'@'+resampfile,combfile,combine='average',reject='avsigclip',                 weight='@'+weightfile,rdnoise='!rdnoise',gain='!gain',                 irafdir=irafdir,error=imcombineerror2,scale='@'+scalefile,zero='@'+zerofile,                 masktype='badvalue',maskvalue=0,bpmasks=mchbase+'_comb.bpm' 
 
if len(imcombineerror2) != 0: 
    printlog,logf,'ERROR in IRAF_IMCOMBINE' 
    printlog,logf,imcombineerror2 
    error = imcombineerror2 
    return 
 
# Convert BPM mask from PL to FITS 
os.remove(mchbase+'_comb.bpm.fits',/allow 
undefine,lines 
cd,current=curdir 
push,lines,'print("")'# first line will be ignored 
push,lines,'cd '+curdir 
push,lines,'imcopy '+mchbase+'_comb.bpm.pl '+mchbase+'_comb.bpm.fits' 
push,lines,'logout' 
tempfile = mktemp('tiraf') 
WRITELINE,tempfile,lines 
IRAF_RUN,tempfile,irafdir,silent=silent,out=out,error=error 
 
# Delete temporary scripts and PL file 
os.remove([tempfile,mchbase+'_comb.bpm.pl'],/allow 
 
 
# Fix the rdnoise and background/sky level and saturate 
#  the bad pixels for DAOPHOT 
#------------------------------------------------------ 
 
# 10/02/12 
# THE IMAGES ARE (1) ZERO-SUBTRACTED, (2) SCALED, AND (3) WEIGHT AVERAGED 
# The algorithm is: 
# 1.) add zero-level correction.  im = im+zero 
# 2.) scale the images.  im = im*scale 
# 3.) take weighted average.  combim=total(weight*im) 
#      there is also clipping that takes place during the averaging 
# The final RDNOISE is essentially: comb_rdnoise = sqrt(total((weights*rdnoise*scale)^2)) 
# A gain that changes from frame to frame could be problematic, 
# but this shouldn't happen since it's the same chip from the same night. 
 
# IMCOMBINE wants rdnoise in electrons and gain in electrons/DN. 
# DAOPHOT expects rdnoise in DN.  That's why mkopt converts 
#  it with the gain.  So we are fine.  The header should have 
#  rdnoise in ELECTRONS. 
 
# page 63-65 of daophot2.pdf shows how you need to modify rdnoise/gain 
# when averaging/summing frames. in observing/mosaic/. 
 
# Load the IMCOMBINE output combined file and BPM 
combim = PHOTRED_READFILE(combfile,combhead) 
badmask = PHOTRED_READFILE(mchbase+'_comb.bpm.fits',maskhead)# 0-good, 1-bad 
 
# Fix the gain 
# For N averaged frames gain(N)=N*gain(1) 
# Leave the gain as is!  We are scaling everything to the reference 
# and using its gain.  It's nearly impossible to figure out the real 
# gain since we are scaling the images and then taking a weighted 
# average with outlier rejection.  Find a gain that properly 
# describes/follows Poisson errors for the final combined image is 
# difficult/impossible.  But that's okay.  This is just for source 
# detection and DAOPHOT FIND just cares about the noise in the 
# background.  We just need to ensure that the sky and rdnoise 
# are correct. 
 
# Fix the rdnoise 
# The final RDNOISE is essentially: comb_rdnoise = sqrt(total((weights*rdnoise*scale)^2)) 
rdnoisearr = fltarr(nfiles) 
for i in range(nfiles): 
    rdnoisearr[i] = PHOTRED_GETRDNOISE(base[i]+'.fits') 
#  the "scales" array here is actually 1/scales used by IMCOMBINE. 
rdnoise = sqrt(np.sum((weights*rdnoisearr/scales)**2)) 
rdnoise = rdnoise > 0.01# must be >=0.01 or it will be 0.00 in the opt file and daophot will crash 
dummy = PHOTRED_GETRDNOISE(combfile,keyword=rdnoisekey)# get keyword 
sxaddpar,combhead,rdnoisekey,rdnoise 
 
# Fix the sky 
# DAOPHOT FIND computes the random error per pixel in ADU as 
# noise = sqrt( sky level/gain + rdnoise^2) 
# So it assumes that the noise in the background is sqrt(sky/gain) 
# in ADU.  We need to set the sky level so this is correct. 
# The final noise should be 
# final sky noise = sqrt(total((weights*scale*sqrt(sky/gain))^2)) 
# So the final sky level should be 
# final sky = total((weights*scale*sqrt(sky/gain))^2)*gain 
gain = PHOTRED_GETGAIN(combfile,keyword=gainkey) 
comb_sky = np.sum((weights*sqrt((sky>0)/gain)/scales)**2)*gain 
# the "scales" array here is actually 1/scale 
combim += float(comb_sky)# keep it float 
 
 
# set the maximum to a "reasonable" level 
# Rescale the image and increase the gain 
if max(combim) > 50000: 
    rescale = 50000./max(combim) 
    combim = combim*rescale 
    sxaddpar,combhead,gainkey,gain/rescale 
    # rdnoise does NOT get modified since it's in electrons 
    # we just need to modify the gain which takes you from ADUs to electrons 
 
maskdatalevel = max(combim) + 10000# set "bad" data level above the highest "good" value 
combim2 = combim*(1-badmask) + maskdatalevel*badmask# set bad pixels to maskdatalevel 
sxaddpar,combhead,'SATURATE',maskdatalevel 
MWRFITS,combim2,combfile,combhead,/create# fits_write can create an empty PDU 
 
# Create the weight map for Sextractor using the BPM output by IMCOMBINE 
#  bad only if bad in ALL images 
weightmap = -2.0*float(badmask == 1) + 1.0 
combweightfile = mchbase+'_comb.mask.fits' 
MWRFITS,weightmap,combweightfile,whead,/create 
 
# NO SCALING of the images for combining 
#--------------------------------------- 
else: 
 
combfile = mchbase+'_comb.fits' 
os.remove(combfile,/allow 
IRAF_IMCOMBINE,'@'+resampfile,combfile,combine='average',reject='avsigclip',                 weight='@'+weightfile,rdnoise='!rdnoise',gain='!gain',                 irafdir=irafdir,error=imcombineerror2 
 
if len(imcombineerror2) != 0: 
    printlog,logf,'ERROR in IRAF_IMCOMBINE' 
    printlog,logf,imcombineerror2 
    error = imcombineerror2 
    return 
 
# Fix the rdnoise and background/sky level and saturate 
#  the bad pixels for DAOPHOT 
#------------------------------------------------------ 
 
# See the explanations for all these steps above!! 
 
# Load the IMCOMBINE output combined file and BPM 
combim = PHOTRED_READFILE(combfile,combhead) 
badmask = PHOTRED_READFILE(mchbase+'_comb.bpm.fits',maskhead)# 0-good, 1-bad 
 
# Fix the rdnoise 
# The final RDNOISE is essentially: comb_rdnoise = sqrt(total((weights*rdnoise)^2)) 
rdnoisearr = fltarr(nfiles) 
for i in range(nfiles): 
    rdnoisearr[i] = PHOTRED_GETRDNOISE(base[i]+'.fits') 
rdnoise = sqrt(np.sum((weights*rdnoisearr)**2)) 
rdnoise = rdnoise > 0.01# must be >=0.01 or it will be 0.00 in the opt file and daophot will crash 
dummy = PHOTRED_GETRDNOISE(combfile,keyword=rdnoisekey)# get keyword 
sxaddpar,combhead,rdnoisekey,rdnoise 
 
# Fix the sky 
# So the final sky level should be 
# final sky = total((weights*scale*sqrt(sky/gain))^2)*gain 
gain = PHOTRED_GETGAIN(combfile,keyword=gainkey) 
comb_sky = np.sum((weights*sqrt(sky/gain))**2)*gain 
# the "scales" array here is actually 1/scale 
combim += comb_sky 
 
 
# set the maximum to a "reasonable" level 
# Rescale the image and increase the gain 
if max(combim) > 50000: 
    rescale = 50000./max(combim) 
    combim = combim*rescale 
    sxaddpar,combhead,gainkey,gain/rescale 
    # rdnoise does NOT get modified since it's in electrons 
    # we just need to modify the gain which takes you from ADUs to electrons 
 
 
# Making Sextractor "weight" map file 
#------------------------------------ 
# masks have 0-bad, 1-good. 
# anything with less than 1.0 is considered bad 
# weight map, -1 is bad, +1 is good 
# "bpm" is the SUM of the bad pixel masks 
# consider a pixel bad that is bad in ANY image 
weightmap = -2.0*float(bpm < nfiles) + 1. 
combweightfile = mchbase+'_comb.mask.fits' 
FITS_WRITE,combweightfile,weightmap,whead 
 
#--------------------------------------------- 
# SATURATE BAD pixels in the COMBINED IMAGE 
# DAOPHOT needs to have the bad pixels "saturated", 
# SExtractor will know which pixels are bad from the "weight" map. 
# 
# We could skip the fiximage.pro step but we still need the 
# individual bpm masks and setting the bad pixels to the background 
# probably helps in the IMALIGN/IMCOMBINE steps. 
printlog,logf,'' 
printlog,logf,'"Saturating" bad pixels in the COMBINED image' 
printlog,logf,'' 
 
 
badmask = float(weightmap < 0.5) 
maskdatalevel = max(combim) + 10000# set "bad" data level above the highest "good" value 
combim2 = combim*(1.0-badmask) + maskdatalevel*badmask# set bad pixels to 100,000 
sxaddpar,combhead,'SATURATE',maskdatalevel 
FITS_WRITE,combfile,combim2,combhead 
 
# no scaling of images for combining 
 
# Add TILETYPE to the combined image 
combhead = PHOTRED_READFILE(combfile,/header) 
sxaddpar,combhead,'AFTILTYP',tile.type 
MODFITS,combfile,0,combhead 
 
# Delete the resampled images 
os.remove(filestr.resampfile,/allow,/quiet 
os.remove(filestr.resampmask,/allow,/quiet 
 
# Make the common source file for the combined image 
#--------------------------------------------------- 
if keyword_set(usecmn): 
print('Combining COMMON SOURCE files for the combined image.' 
# Loop through the files and convert to coordinates to the comined file 
undefine,allcmn 
for i in range(nfiles): 
    cmnfile1 = os.path.basename(filestr[i].fitsfile,'.fits')+'.cmn.lst' 
    if os.path.exists(cmnfile1) == 1: 
        cmn1 = IMPORTASCII(cmnfile1,fieldnames=['id','x','y','mag','err','sky','skysig','sharp','round','round2'],                         skipline=3,/silent) 
        READLINE,cmnfile1,cmnlines1 
        coohead1 = cmnlines1[0:1] 
        ncmn1 = len(cmn1) 
        # Get coordinates on the resampled/combined image grid 
        out = trans_coo(cmn1.x,cmn1.y,filestr[i].resamptrans) 
        newx = reform(out[0,:]) 
        newy = reform(out[1,:]) 
        cmn1.x = newx 
        cmn1.y = newy 
        if len(allcmn) == 0: 
            allcmn = cmn1 
        else: 
            # Remove any duplicates 
            SRCMATCH,allcmn.x,allcmn.y,cmn1.x,cmn1.y,2.0,ind1,ind2,count=nmatch 
            if nmatch > 0: 
                if nmatch < ncmn1: 
                    remove,ind2,cmn1 
                else: 
                    undefine,cmn1 
            if len(cmn1) > 0 : 
                push,allcmn,cmn1 
if len(allcmn) > 0: 
    WRITECOL,mchbase+'_comb.cmn.lst',allcmn.id,allcmn.x,allcmn.y,allcmn.mag,allcmn.err,allcmn.sky,allcmn.skysig,             allcmn.sharp,allcmn.round,allcmn.round2,fmt='(I7,2F9.2,3F9.3,F9.2,3F9.3)' 
    WRITELINE,mchbase+'_comb.cmn.lst',[coohead1,''],/prepend# prepend the COO header 
 else printlog,logf,'No common file to combine' 
     
    # Don't use common file 
else: 
    # Make sure it doesn't exist otherwise it will be used 
    os.remove(mchbase+'_comb.cmn.lst',/allow 
 
BOMB: 
 
if keyword_set(stp) : 
    import pdb; pdb.set_trace() 
 
 
#!/usr/bin/env python

import os
import time
import numpy as np

#+ 
# 
# ALLFRAME_COMBINE_ORIG 
# 
# This is the old version of the code that combines images ALLFRAME. 
# 
# You need to have run daophot, allstar, daomatch and daomaster 
# already.  There needs to be a fits, opt, als.opt, ap and als 
# file for each file in the MCH file.  There also needs to be an 
# associated RAW file for the MCH file. 
# 
# INPUTS: 
#  file            The MCH filename 
#  =nocmbimscale   Don't scale the images when combining them.  Not 
#                    recommended, but the old way of doing it.  Bright 
#                    stars can be missed this way. 
#  /combtrim       Trim the combined images to the overlapping region. 
#                    This used to be the default, but now the default 
#                    is to keep the entire original region. 
#  =scriptsdir     The directory that contains all of the necessary scripts. 
#  =irafdir        The IRAF home directory. 
#  =logfile        A logfile to print to output to. 
#  /fake           Run for artificial star tests. 
#  /usecmn         Use the cmn.lst file of the reference image for the 
#                    combined image. 
#  /stp            Stop at the end of the program 
# 
# OUTPUTS: 
#  The combined image and mask: 
#    FILEBASE.comb.fits       combined image 
#    FILEBASE.comb.bpm.fits   bad pixel mask 
#    FILEBASE.comb.mask.fits  SExtractor weight map (very similar to bpm) 
#    FILEBASE.comb.mch        The transformations of the individual frames 
#                               to the combined frame. 
#  =maskdatalevel  The "bad" data level above the highest "good" value 
#                    in the combined image. 
#  =error          The error message, if there was one, else undefined 
# 
# USAGE: 
#  IDL>allframe_combine_orig,'ccd1001.mch',scriptsdir='/net/home/dln5q/daophot/',finditer=2 
# 
# 
# By D.Nidever   February 2008 
# Automation of steps and scripts by J.Ostheimer and Rachael Beaton 
#- 
 
 
def allframe_combine_orig,file,stp=stp,scriptsdir=scriptsdir,error=error,logfile=logfile,                          irafdir=irafdir,satlevel=satlevel,nocmbimscale=nocmbimscale,                          trimcomb=trimcomb,fake=fake,maskdatalevel=maskdatalevel,                          xoff=xoff,yoff=yoff,usecmn=usecmn 
     
    COMMON photred,setup 
     
    undefine,error 
     
    # Not enough inputs 
    nfile = len(file) 
    if (nfile == 0): 
        print('Syntax - allframe_combine_orig,file,stp=stp,scriptsdir=scriptsdir,satlevel=satlevel,' 
        print('                  ,nocmbimscale=nocmbimscale,error=error,logfile=logfile,' 
        print('                  irafdir=irafdir,trimcomb=trimcomb,usecmn=usecmn,fake=fake' 
        return 
     
    # Logfile 
    if keyword_set(logfile): 
        logf=logfile 
    else: 
        logf=-1 
     
    # Error Handling 
    #------------------ 
    # Establish error handler. When errors occur, the index of the 
    # error is returned in the variable Error_status: 
    CATCH, Error_status 
     
    #This statement begins the error handler: 
    if (Error_status != 0): 
        print('ALLFRAME_COMBINE_ORIG ERROR: ', !ERROR_STATE.MSG 
        error = !ERROR_STATE.MSG 
        CATCH, /CANCEL 
        return 
     
    # Saturation level 
    if len(satlevel) == 0 : 
        satlevel=6e4 
     
    # Scaling of the images to be combined 
    if len(nocmbimscale) == 0 : 
        nocmbimscale=0 
     
    # Getting scripts directory and iraf directory 
    nsetup = len(setup) 
    if nsetup > 0: 
        scriptsdir = READPAR(setup,'SCRIPTSDIR') 
        irafdir = READPAR(setup,'IRAFDIR') 
     
     
    # No irafdir 
    if len(scriptsdir) == 0: 
        printlog,logf,'SCRIPTSDIR NOT INPUT' 
        error = 'SCRIPTSDIR NOT INPUT' 
        return 
     
    # Run CHECK_IRAF.PRO to make sure that you can run IRAF from IDL 
    #--------------------------------------------------------------- 
    CHECK_IRAF,iraftest,irafdir=irafdir 
    if iraftest == 0: 
        print('IRAF TEST FAILED.  EXITING' 
        return 
     
     
    # No scriptsdir 
    if len(scriptsdir) == 0: 
        printlog,logf,'SCRIPTSDIR NOT INPUT' 
        error = 'SCRIPTSDIR NOT INPUT' 
        return 
    # Check if the scripts exist in the current directory 
    scripts = ['getpsf.sh','photo.opt','apcor.opt','lstfilter.py','goodpsf.pro','allframe.opt',           'default.sex','default.param','default.nnw','default.conv'] 
    nscripts = len(scripts) 
    # Loop through the scripts 
    for i in range(nscripts): 
        info = FILE_INFO(scriptsdir+'/'+scripts[i]) 
        curinfo = FILE_INFO(scripts[i]) 
         
        # No file 
        if info.exists == 0 or info.size == 0: 
            printlog,logf,scriptsdir+'/'+scripts[i],' NOT FOUND or EMPTY' 
            error = scriptsdir+'/'+scripts[i]+' NOT FOUND or EMPTY' 
            return 
         
        # Check if the two files are the same size, if not copy it 
        if info.size != curinfo.size: 
            FILE_COPY,info.name,curinfo.name,/overwrite 
# scripts loop 
     
     
    printlog,logf,'Combining images in ',file 
     
    # FILENAME 
    mchfile = os.path.basename(file) 
    mchdir = os.path.dirname(file) 
    mchbase = os.path.basename(file,'.mch') 
     
     
    # CD to the directory 
    cd,current=curdir 
    cd,mchdir 
     
     
    # Check that the mch, als, and opt files exist 
    mchtest = os.path.exists(mchfile) 
    if mchtest == 0: 
        printlog,logf,mchfile,' NOT FOUND' 
        return 
     
    # Checking RAW file 
    rawtest = os.path.exists(mchbase+'.raw') 
    if rawtest == 0: 
        printlog,logf,mchbase+'.raw NOT FOUND' 
        return 
     
     
    ############################################ 
    # CHECK NECESSARY FILES 
     
    # Load the MCH file 
    LOADMCH,mchfile,files,trans 
     
    # Check that the fits, als, opt, and psf files exist 
    nfiles = len(files) 
    for i in range(nfiles): 
        #dir = file_dirname(mchfile) 
        base = os.path.basename(files[i],'.als') 
         
        # Checking FITS file 
        fitstest = os.path.exists(base+'.fits') 
        if fitstest == 0: 
            printlog,logf,base+'.fits NOT FOUND' 
            return 
         
        # Checking OPT file 
        opttest = os.path.exists(base+'.opt') 
        if opttest == 0: 
            printlog,logf,base+'.opt NOT FOUND' 
            return 
         
        # Checking ALS.OPT file 
        alsopttest = os.path.exists(base+'.als.opt') 
        if alsopttest == 0: 
            printlog,logf,base+'.als.opt NOT FOUND' 
            return 
         
        # Checking AP file 
        aptest = os.path.exists(base+'.ap') 
        if aptest == 0: 
            printlog,logf,base+'.ap NOT FOUND' 
            return 
         
        # Checking ALS file 
        alstest = os.path.exists(base+'.als') 
        if alstest == 0: 
            printlog,logf,base+'.als NOT FOUND' 
            return 
         
     
    # FAKE, check that we have all the files that we need 
    if keyword_set(fake): 
        # weights, scale, zero, comb_psf, _shift.mch 
        chkfiles = mchbase+['.weights','.scale','.zero','_comb.psf','_shift.mch'] 
        bdfiles , = np.where(os.path.exists(chkfiles) == 0,nbdfiles) 
        if nbdfiles > 0: 
            error = 'FAKE.  Some necessary files not found. '+strjoin(chkfiles[bdfiles],' ') 
            printlog,logf,error 
            return 
     
     
     
    ############################################ 
    # STEP 1: IMALIGN PREP 
     
    printlog,logf,'' 
    printlog,logf,'Step A: Getting Weights' 
    printlog,logf,'-----------------------' 
     
    #----------------------------------- 
    # Computs Weights 
    if not keyword_set(fake): 
        ALLFRAME_GETWEIGHTS,mchfile,weights,scales,sky#,raw2=raw2 
        invscales = 1.0/scales 
        bdscale , = np.where(scales < 1e-5 or invscales > 900,nbdscale) 
        if nbdscale > 0: 
            scales[bdscale] = 1.0 
            invscales[bdscale] = 1.0 
            weights[bdscale] = 0.0 
        weightfile = mchbase+'.weights' 
        WRITECOL,weightfile,weights,fmt='(F10.6)' 
        scalefile = mchbase+'.scale' 
        WRITECOL,scalefile,invscales,fmt='(F10.6)'# want to scale it UP 
        zerofile = mchbase+'.zero' 
        WRITECOL,zerofile,-sky,fmt='(F12.4)'# want to remove the background, set to 1st frame 
         
        # FAKE, use existing ones 
    else: 
        weightfile = mchbase+'.weights' 
        scalefile = mchbase+'.scale' 
        zerofile = mchbase+'.zero' 
        READCOL,weightfile,weights,format='F',/silent 
        READCOL,scalefile,invscales,format='F',/silent 
        scales = 1.0/invscales 
        READCOL,zerofile,sky,format='F',/silent 
        sky = -sky 
     
     
    #--------------------------------------- 
    # Get X/Y translations using DAOMASTER 
    #  NO ROTATION ONLY SHIFTS 
    #  Need these shifts for IMSHIFT 
    shiftmch = mchbase+'_shift' 
    if not keyword_set(fake): 
        print('Measuring X/Y shifts' 
        FILE_COPY,mchbase+'.mch',shiftmch+'.mch',/overwrite,/allow 
        # Make the DAOMASTER script 
        undefine,cmdlines 
        PUSH,cmdlines,'#!/bin/csh' 
        PUSH,cmdlines,'set input=${1}' 
        PUSH,cmdlines,'daomaster <<:NE' 
        PUSH,cmdlines,'${input}.mch' 
        PUSH,cmdlines,'1,1,1' 
        PUSH,cmdlines,'99.' 
        PUSH,cmdlines,'2' 
        PUSH,cmdlines,'10' 
        PUSH,cmdlines,'5' 
        PUSH,cmdlines,'4' 
        PUSH,cmdlines,'3' 
        PUSH,cmdlines,'2' 
        PUSH,cmdlines,'1' 
        PUSH,cmdlines,'0' 
        PUSH,cmdlines,'n' 
        PUSH,cmdlines,'n' 
        PUSH,cmdlines,'n' 
        PUSH,cmdlines,'n' 
        PUSH,cmdlines,'y' 
        PUSH,cmdlines,'' 
        PUSH,cmdlines,'' 
        PUSH,cmdlines,'n' 
        PUSH,cmdlines,'n' 
        PUSH,cmdlines,'n' 
        PUSH,cmdlines,'DONE' 
        tempscript = MKTEMP('daomaster')# absolute filename 
        WRITELINE,tempscript,cmdlines 
        FILE_CHMOD,tempscript,'755'o 
        # Run DAOMASTER 
        cmd2 = tempscript+' '+shiftmch 
        SPAWN,cmd2,out2,errout2 
        # Remove temporary DAOMASTER script 
        os.remove(tempscript,/allow_non 
    LOADMCH,shiftmch+'.mch',files2,trans2 
     
    xshift = reform(trans2[:,0]) 
    yshift = reform(trans2[:,1]) 
    xyshifts = [[xshift],[yshift]] 
    printlog,logf,'Image shifts' 
    for i in range(nfiles): 
        printlog,logf,files[i],xshift[i],yshift[i] 
     
     
    #----------------------------------- 
    # Create imalign prep files 
    # This is done by preimalign_k.sh 
    # Need an input list of fits files 
    # Need an output list of fits files 
    # Shift file 
    base = os.path.basename(files,'.als') 
    fitsfiles = base+'.fits' 
    outfiles = base+'.shft.fits' 
    infile = mchbase+'.inlist' 
    outfile = mchbase+'.outlist' 
    #WRITELINE,infile,fitsfiles   ; this is done below now with the temp files 
    WRITELINE,outfile,outfiles 
    # Remove outfiles 
    os.remove(outfiles,/allow 
    # shift list 
    shiftfile = mchbase+'.shift' 
    WRITECOL,shiftfile,xshift,yshift,fmt='(2F15.4)' 
     
     
    # Make temporary files for bad pixel fixing and combining 
    #  FIXIMAGES doesn't work properly on the shifted images 
    #  because the interpolation can bring the bad pixel values down 
    #  below the saturation threshold, and we don't want to touch 
    #  the original images. 
    tempfits = base+'.temp.fits' 
    FILE_COPY,fitsfiles,tempfits,/overwrite,/allow 
    WRITELINE,infile,tempfits 
     
     
     
    ############################################ 
    # STEP B: FIX BAD PIXELS 
    printlog,logf,'' 
    printlog,logf,'Step B: Fixing bad pixels' 
    printlog,logf,'-------------------------' 
    FIXIMAGES,'@'+infile,satlevel=satlevel#6e4 
    # This also makes the FILE.mask.fits files for each image 
     
    # Find the maximum saturation level 
    satlevelarr = fltarr(nfiles) 
    for i in range(nfiles): 
        #head = headfits(base[i]+'.fits') 
        im = PHOTRED_READFILE(base[i]+'.fits',head) 
        saturate = sxpar(head,'SATURATE',count=nsaturate,/silent) 
        if nsaturate == 0 : 
            saturate=max(im)-1000. 
        satlevelarr[i] = saturate 
    maxsatlevel = max(satlevelarr) 
     
     
    ############################################ 
    # STEP C: IMALIGN 
    #  This figures out the X/Y-shifts between the images 
    #  and creates the shifted images (".shft.fits") 
    printlog,logf,'' 
    printlog,logf,'Step C: IMALIGN' 
    printlog,logf,'---------------' 
    reffile = mchbase+'.fits' 
     
    # IMALIGN basically is a script that runs: 
    #  IMCENTROID - to compute the shifts 
    #  IMSHIFT - shifts the images 
    #  IMCOPY - trims the images 
     
     
    # First, shift the images 
    printlog,logf,'Shifting the images' 
    IRAF_IMSHIFT,'@'+infile,'@'+outfile,shifts_file=shiftfile,interp_type='linear',             boundary_type='constant',constant=0,irafdir=irafdir,error=imshifterror 
    if len(imshifterror) != 0: 
        printlog,logf,'ERROR in IRAF_IMSHIFT' 
        printlog,logf,imshifterror 
        error = imshifterror 
        return 
     
    #stop 
     
    # Trim the images 
    if keyword_set(trimcomb): 
         
        # Calculate the trim section 
        hd = PHOTRED_READFILE(reffile,/header) 
        xsize = lonarr(nfiles)+sxpar(hd,'NAXIS1',/silent) 
        ysize = lonarr(nfiles)+sxpar(hd,'NAXIS2',/silent) 
        IA_TRIM,xshift,yshift,xsize,ysize,trimsection 
        xoff = trimsection[0]-1 
        yoff = trimsection[2]-1 
         
        # Trim the shifted images 
        printlog,logf,'Trimming the shifted images' 
        xstart = trimsection[0]-1 
        ximport pdb; pdb.set_trace() = trimsection[1]-1 
        ystart = trimsection[2]-1 
        yimport pdb; pdb.set_trace() = trimsection[3]-1 
         
        for i in range(nfiles): 
            im = PHOTRED_READFILE(outfiles[i],head) 
            newim = im[xstart:ximport pdb; pdb.set_trace(),ystart:yimport pdb; pdb.set_trace()] 
            MWRFITS,newim,outfiles[i],head,/create,/silent 
        # could also use IRAF_IMCOPY here instead 
         
         
        # Don't trim the images 
    else: 
        hd = PHOTRED_READFILE(reffile,/header) 
        xsize = sxpar(hd,'NAXIS1',/silent) 
        ysize = sxpar(hd,'NAXIS2',/silent) 
        trimsection = [1,xsize,1,ysize] 
        xoff = 0 
        yoff = 0 
     
    # Delete the temporary FITS files 
    os.remove(tempfits,/allow 
     
     
     
    ############################################ 
    # STEP D: MAKE BAD PIXEL/WEIGHT MAP 
    # in how many images does the pixel need to be bad?? 
    # 1. shift the masks (created by FIXIMAGES.PRO) 
    #     using the shifts from IRAF_IMALIGN 
    # 2. trim the masks 
    # 3. combine the masks 
    printlog,logf,'' 
    printlog,logf,'Step D: Making Bad Pixel Mask' 
    printlog,logf,'-----------------------------' 
     
    # Make lists 
    maskfiles = os.path.dirname(tempfits)+'/'+os.path.basename(tempfits,'.fits')+'.mask.fits' 
    outmaskfiles = os.path.dirname(tempfits)+'/'+os.path.basename(tempfits,'.fits')+'.mask.shft.fits' 
    maskinfile = mchbase+'.maskinlist' 
    maskoutfile = mchbase+'.maskoutlist' 
    maskshiftsfile = mchbase+'.maskshifts' 
    WRITELINE,maskinfile,maskfiles 
    WRITELINE,maskoutfile,outmaskfiles 
    os.remove(outmaskfiles,/allow 
    strxyshifts = strarr(nfiles) 
    for i in range(nfiles): 
        strxyshifts[i] = strjoin(reform(xyshifts[i,:]),'  ') 
    WRITELINE,maskshiftsfile,strxyshifts 
     
    # Run IMSHIFT 
    #  set boundary to 0=BAD 
    printlog,logf,'Shifting masks' 
    undefine,iraflines 
    push,iraflines,'print("")'# first line will be ignored 
    push,iraflines,'cd '+curdir 
    push,iraflines,'images' 
    push,iraflines,'imgeom' 
    push,iraflines,'imshift("@'+maskinfile+'","@'+maskoutfile+'",shifts_file="'+maskshiftsfile+'",'+               'interp_type="linear",boundary_typ="constant",constant=0)' 
    push,iraflines,'logout' 
    imshiftscript = curdir+'/'+mchbase+'.imshift' 
    WRITELINE,imshiftscript,iraflines 
    IRAF_RUN,imshiftscript,irafdir,out=out,/silent,error=iraferror 
    if len(iraferror) != 0: 
        printlog,logf,'ERROR in running IMSHIFT with IRAF_RUN' 
        printlog,logf,iraferror 
        error = iraferror 
        return 
     
    # Trim 
    if keyword_set(trimcomb): 
        printlog,logf,'Trimming masks' 
        xstart = trimsection[0]-1# should be same as xoff 
        ximport pdb; pdb.set_trace() = trimsection[1]-1 
        ystart = trimsection[2]-1# should be same as yoff 
        yimport pdb; pdb.set_trace() = trimsection[3]-1 
         
        for i in range(nfiles): 
            im = PHOTRED_READFILE(outmaskfiles[i],head) 
            sz = size(im) 
            newim = im[xstart:ximport pdb; pdb.set_trace(),ystart:yimport pdb; pdb.set_trace()] 
            # Add LTV1/LTV2 to the header 
            #  these are IRAF keywords to convert from logical to physical coords 
            ltv1 = sxpar(head,'LTV1',/silent)# 0 if not found 
            ltv2 = sxpar(head,'LTV2',/silent)# 0 if not found 
            sxaddpar,head,'LTV1',ltv1-xstart 
            sxaddpar,head,'LTV2',ltv2-ystart 
            MWRFITS,newim,outmaskfiles[i],head,/create,/silent 
     
    # Combining masks 
    printlog,logf,'Combining masks' 
    undefine,bpm 
    for i in range(nfiles): 
        im = PHOTRED_READFILE(outmaskfiles[i],head) 
        if i == 0: 
            bpm = im 
            whead = head 
        else: 
            bpm = bpm+im 
    #bpm = bpm/float(nfiles) 
     
    # masks have 0-bad, 1-good. 
    # anything with less than 1.0 is considered bad 
    # weight map, -1 is bad, +1 is good 
    #weightmap = -2.0*float(bpm lt nfiles) + 1. 
    #combweightfile = mchbase+'_comb.mask.fits' 
    #FITS_WRITE,combweightfile,weightmap,whead 
    # 
    # THIS IS NOW DONE BELOW AFTER THE IMAGE IS COMBINED 
    # DEPENDING ON IF THE IMAGES ARE SCALED OR NOT!!! 
     
     
    ############################################ 
    # STEP E: COMBINE IMAGES 
    printlog,logf,'' 
    printlog,logf,'Step E: IMCOMBINE' 
    printlog,logf,'-----------------' 
     
    # SCALE the images for combining 
    #------------------------------- 
    if not keyword_set(nocmbimscale): 
         
        # Put BPM mask names in file headers 
        #  these will be used by IMCOMBINE 
        for i in range(nfiles): 
            head = PHOTRED_READFILE(outfiles[i],/header) 
            sxaddpar,head,'BPM',outmaskfiles[i] 
            MODFITS,outfiles[i],0,head 
         
        # Combine the frames WITH scaling/offset/masking, for the bright stars 
        #printlog,logf,'Creating SCALED image' 
        combfile = mchbase+'_comb.fits' 
        os.remove(combfile,/allow 
        os.remove(mchbase+'_comb.bpm.pl',/allow 
        IRAF_IMCOMBINE,'@'+outfile,combfile,combine='average',reject='avsigclip',                 weight='@'+weightfile,rdnoise='!rdnoise',gain='!gain',                 irafdir=irafdir,error=imcombineerror2,scale='@'+scalefile,zero='@'+zerofile,                 masktype='badvalue',maskvalue=0,bpmasks=mchbase+'_comb.bpm' 
         
        if len(imcombineerror2) != 0: 
            printlog,logf,'ERROR in IRAF_IMCOMBINE' 
            printlog,logf,imcombineerror2 
            error = imcombineerror2 
            return 
         
        # Convert BPM mask from PL to FITS 
        os.remove(mchbase+'_comb.bpm.fits',/allow 
        undefine,lines 
        cd,current=curdir 
        push,lines,'print("")'# first line will be ignored 
        push,lines,'cd '+curdir 
        push,lines,'imcopy '+mchbase+'_comb.bpm.pl '+mchbase+'_comb.bpm.fits' 
        push,lines,'logout' 
        tempfile = mktemp('tiraf') 
        WRITELINE,tempfile,lines 
        IRAF_RUN,tempfile,irafdir,silent=silent,out=out,error=error 
         
        # Delete temporary scripts and PL file 
        os.remove([tempfile,mchbase+'_comb.bpm.pl'],/allow 
         
         
        # Fix the rdnoise and background/sky level and saturate 
        #  the bad pixels for DAOPHOT 
        #------------------------------------------------------ 
         
        # 10/02/12 
        # THE IMAGES ARE (1) ZERO-SUBTRACTED, (2) SCALED, AND (3) WEIGHT AVERAGED 
        # The algorithm is: 
        # 1.) add zero-level correction.  im = im+zero 
        # 2.) scale the images.  im = im*scale 
        # 3.) take weighted average.  combim=total(weight*im) 
        #      there is also clipping that takes place during the averaging 
        # The final RDNOISE is essentially: comb_rdnoise = sqrt(total((weights*rdnoise*scale)^2)) 
        # A gain that changes from frame to frame could be problematic, 
        # but this shouldn't happen since it's the same chip from the same night. 
         
        # IMCOMBINE wants rdnoise in electrons and gain in electrons/DN. 
        # DAOPHOT expects rdnoise in DN.  That's why mkopt converts 
        #  it with the gain.  So we are fine.  The header should have 
        #  rdnoise in ELECTRONS. 
         
        # page 63-65 of daophot2.pdf shows how you need to modify rdnoise/gain 
        # when averaging/summing frames. in observing/mosaic/. 
         
        # Load the IMCOMBINE output combined file and BPM 
        combim = PHOTRED_READFILE(combfile,combhead) 
        badmask = PHOTRED_READFILE(mchbase+'_comb.bpm.fits',maskhead)# 0-good, 1-bad 
         
         
        # Fix the gain 
        # For N averaged frames gain(N)=N*gain(1) 
        # Leave the gain as is!  We are scaling everything to the reference 
        # and using its gain.  It's nearly impossible to figure out the real 
        # gain since we are scaling the images and then taking a weighted 
        # average with outlier rejection.  Find a gain that properly 
        # describes/follows Poisson errors for the final combined image is 
        # difficult/impossible.  But that's okay.  This is just for source 
        # detection and DAOPHOT FIND just cares about the noise in the 
        # background.  We just need to ensure that the sky and rdnoise 
        # are correct. 
         
        # Fix the rdnoise 
        # The final RDNOISE is essentially: comb_rdnoise = sqrt(total((weights*rdnoise*scale)^2)) 
        rdnoisearr = fltarr(nfiles) 
        for i in range(nfiles): 
            rdnoisearr[i] = PHOTRED_GETRDNOISE(base[i]+'.fits') 
        #  the "scales" array here is actually 1/scales used by IMCOMBINE. 
        rdnoise = sqrt(np.sum((weights*rdnoisearr/scales)**2)) 
        rdnoise = rdnoise > 0.01# must be >=0.01 or it will be 0.00 in the opt file and daophot will crash 
        dummy = PHOTRED_GETRDNOISE(combfile,keyword=rdnoisekey)# get keyword 
        sxaddpar,combhead,rdnoisekey,rdnoise 
         
        # Fix the sky 
        # DAOPHOT FIND computes the random error per pixel in ADU as 
        # noise = sqrt( sky level/gain + rdnoise^2) 
        # So it assumes that the noise in the background is sqrt(sky/gain) 
        # in ADU.  We need to set the sky level so this is correct. 
        # The final noise should be 
        # final sky noise = sqrt(total((weights*scale*sqrt(sky/gain))^2)) 
        # So the final sky level should be 
        # final sky = total((weights*scale*sqrt(sky/gain))^2)*gain 
        gain = PHOTRED_GETGAIN(combfile,keyword=gainkey) 
        comb_sky = np.sum((weights*sqrt((sky>0)/gain)/scales)**2)*gain 
        # the "scales" array here is actually 1/scale 
        combim += float(comb_sky)# keep it float 
         
         
        # set the maximum to a "reasonable" level 
        # Rescale the image and increase the gain 
        if max(combim) > 50000: 
            rescale = 50000./max(combim) 
            combim = combim*rescale 
            sxaddpar,combhead,gainkey,gain/rescale 
            # rdnoise does NOT get modified since it's in electrons 
            # we just need to modify the gain which takes you from ADUs to electrons 
         
        maskdatalevel = max(combim) + 10000# set "bad" data level above the highest "good" value 
        combim2 = combim*(1-badmask) + maskdatalevel*badmask# set bad pixels to maskdatalevel 
        MWRFITS,combim2,combfile,combhead,/create# fits_write can create an empty PDU 
         
        # Create the weight map for Sextractor using the BPM output by IMCOMBINE 
        #  bad only if bad in ALL images 
        weightmap = -2.0*float(badmask == 1) + 1.0 
        combweightfile = mchbase+'_comb.mask.fits' 
        MWRFITS,weightmap,combweightfile,whead,/create 
         
        # NO SCALING of the images for combining 
        #--------------------------------------- 
    else: 
         
        combfile = mchbase+'_comb.fits' 
        os.remove(combfile,/allow 
        IRAF_IMCOMBINE,'@'+outfile,combfile,combine='average',reject='avsigclip',                 weight='@'+weightfile,rdnoise='!rdnoise',gain='!gain',                 irafdir=irafdir,error=imcombineerror2 
         
        if len(imcombineerror2) != 0: 
            printlog,logf,'ERROR in IRAF_IMCOMBINE' 
            printlog,logf,imcombineerror2 
            error = imcombineerror2 
            return 
         
        # Fix the rdnoise and background/sky level and saturate 
        #  the bad pixels for DAOPHOT 
        #------------------------------------------------------ 
         
        # See the explanations for all these steps above!! 
         
        # Load the IMCOMBINE output combined file and BPM 
        combim = PHOTRED_READFILE(combfile,combhead) 
        badmask = PHOTRED_READFILE(mchbase+'_comb.bpm.fits',maskhead)# 0-good, 1-bad 
         
        # Fix the rdnoise 
        # The final RDNOISE is essentially: comb_rdnoise = sqrt(total((weights*rdnoise)^2)) 
        rdnoisearr = fltarr(nfiles) 
        for i in range(nfiles): 
            rdnoisearr[i] = PHOTRED_GETRDNOISE(base[i]+'.fits') 
        rdnoise = sqrt(np.sum((weights*rdnoisearr)**2)) 
        rdnoise = rdnoise > 0.01# must be >=0.01 or it will be 0.00 in the opt file and daophot will crash 
        dummy = PHOTRED_GETRDNOISE(combfile,keyword=rdnoisekey)# get keyword 
        sxaddpar,combhead,rdnoisekey,rdnoise 
         
        # Fix the sky 
        # So the final sky level should be 
        # final sky = total((weights*scale*sqrt(sky/gain))^2)*gain 
        gain = PHOTRED_GETGAIN(combfile,keyword=gainkey) 
        comb_sky = np.sum((weights*sqrt(sky/gain))**2)*gain 
        # the "scales" array here is actually 1/scale 
        combim += comb_sky 
         
         
        # set the maximum to a "reasonable" level 
        # Rescale the image and increase the gain 
        if max(combim) > 50000: 
            rescale = 50000./max(combim) 
            combim = combim*rescale 
            sxaddpar,combhead,gainkey,gain/rescale 
            # rdnoise does NOT get modified since it's in electrons 
            # we just need to modify the gain which takes you from ADUs to electrons 
         
         
        # Making Sextractor "weight" map file 
        #------------------------------------ 
        # masks have 0-bad, 1-good. 
        # anything with less than 1.0 is considered bad 
        # weight map, -1 is bad, +1 is good 
        # "bpm" is the SUM of the bad pixel masks 
        # consider a pixel bad that is bad in ANY image 
        weightmap = -2.0*float(bpm < nfiles) + 1. 
        combweightfile = mchbase+'_comb.mask.fits' 
        FITS_WRITE,combweightfile,weightmap,whead 
         
        #--------------------------------------------- 
        # SATURATE BAD pixels in the COMBINED IMAGE 
        # DAOPHOT needs to have the bad pixels "saturated", 
        # SExtractor will know which pixels are bad from the "weight" map. 
        # 
        # We could skip the fiximage.pro step but we still need the 
        # individual bpm masks and setting the bad pixels to the background 
        # probably helps in the IMALIGN/IMCOMBINE steps. 
        printlog,logf,'' 
        printlog,logf,'"Saturating" bad pixels in the COMBINED image' 
        printlog,logf,'' 
         
         
        badmask = float(weightmap < 0.5) 
        maskdatalevel = max(combim) + 10000# set "bad" data level above the highest "good" value 
        combim2 = combim*(1.0-badmask) + maskdatalevel*badmask# set bad pixels to 100,000 
        FITS_WRITE,combfile,combim2,combhead 
         
# no scaling of images for combining 
     
    # Add TILETYPE to the combined image 
    combhead = PHOTRED_READFILE(combfile,/header) 
    sxaddpar,combhead,'AFTILTYP','ORIG' 
    MODFITS,combfile,0,combhead 
     
    # Delete the shifted images 
    READLINE,outfile,shiftedfiles 
    os.remove(shiftedfiles,/allow,/quiet 
     
    # Delete mask files 
    os.remove([maskfiles,outmaskfiles],/allow 
    os.remove([maskinfile,maskoutfile,maskshiftsfile,imshiftscript],/allow 
     
    # Copy the original MCH file to COMB.MCH 
    FILE_COPY,file,mchdir+'/'+mchbase+'.comb.mch',/allow,/over 
     
    # Using CMN.LST of reference frame if it exists 
    if os.path.exists(mchbase+'.cmn.lst') and keyword_set(usecmn): 
        print('Using reference image COMMON SOURCE file' 
        FILE_COPY,mchbase+'.cmn.lst',mchbase+'_comb.cmn.lst',/over,/allow 
     
    # CD back to the original directory 
    cd,curdir 
     
    BOMB: 
     
    if keyword_set(stp) : 
        import pdb; pdb.set_trace() 
     
 
