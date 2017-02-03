pro complete_stargalaxy

; Run the ASTs through the star galaxy separation

field = 'Field100'
dir = '/datalab/users/dnidever/smash/cp/red/photred/addfakes/'+field+'/'
comb = mrdfits(dir+field+'_complete.fits.gz',1)
; only use recovered ASTs
gdobj = where(comb.ast eq 1 and comb.recovered eq 1,ngdobj)
obj = comb[gdobj]

; Copied from papers/smash_laf/stargalaxy_separation.pro

;; Select the stars
;; prob can remove bright stars, require g>20
;gdstars = where(abs(obj.sharp) lt 1 and obj.chi lt 2 and obj.prob gt 0.2 and $
;                obj.ndet gt 5 and obj.depthflag gt 1,ngdstars)
ndet = obj.ndetu+obj.ndetg+obj.ndetr+obj.ndeti+obj.ndetz
add_tag,obj,'ndet',0L,obj
obj.ndet = ndet
gdstars = where(abs(obj.sharp) lt 1 and obj.chi lt 2 and obj.prob gt 0.2 and $
                obj.ndet gt 5,ngdstars)
obj = obj[gdstars]

;; Morphology cuts
smashred_morphcuts,obj,ind1,nsig=3
obj1 = obj[ind1]
  
;; Color-color cuts
smashred_2cdcuts,obj1,ind2
obj2 = obj1[ind2]
; this cuts out some stars with g>24 on Field56

; YEP, THAT CUT OUT ABOUT 25% OF THE STARS

stop

; Figure out the completeness                                                                                                                      
gdall = where(comb.ast eq 1,ngdall)
gdrecover = where(comb.ast eq 1 and comb.recovered eq 1,ngdrecover)
dx = 0.2
dy = 0.4
xr = [-1,3.5]
yr = [17.0,27.0]
hess,comb[gdall].ast_g-comb[gdall].ast_i,comb[gdall].ast_g,dum,imall,dx=dx,dy=dy,xr=xr,yr=yr,xarr=xarr,yarr=yarr,/noplot
hess,comb[gdrecover].ast_g-comb[gdrecover].ast_i,comb[gdrecover].ast_g,dum,imrec,dx=dx,dy=dy,xr=xr,yr=yr,xarr=xarr,yarr=yarr,/noplot
;displayc,float(imrec)/(imall>1),xarr,yarr,/yflip,xtit='g-i',ytit='g',tit='Completeness'                                                           

; Make some figures                                                                                                                                
setdisp
;loadcol,3                                                                                                                                         
!p.font = 0
; Input ASTS                                                                                                                                       
file = dir+'plots/'+field+'_input'
ps_open,file,/color,thick=4,/encap
device,/inches,xsize=8.5,ysize=10.5
displayc,imall,xarr,yarr,/yflip,xtit='g-i',ytit='g',tit='Input ASTs for '+field,charsize=1.3
ps_close
ps2png,file+'.eps',/eps
spawn,['epstopdf',file+'.eps'],/noshell
; Recovered                                                                                                                                        
file = dir+'plots/'+field+'_recovered'
ps_open,file,/color,thick=4,/encap
device,/inches,xsize=8.5,ysize=10.5
displayc,imrec,xarr,yarr,/yflip,xtit='g-i',ytit='g',tit='Recovered ASTs for '+field,charsize=1.3
ps_close
ps2png,file+'.eps',/eps
spawn,['epstopdf',file+'.eps'],/noshell
; Completeness                                                                                                                                     
file = dir+'plots/'+field+'_completeness'
ps_open,file,/color,thick=4,/encap
device,/inches,xsize=8.5,ysize=10.5
displayc,float(imrec)/(imall>1),xarr,yarr,/yflip,xtit='g-i',ytit='g',tit='Completeness for '+field,charsize=1.3
ps_close
ps2png,file+'.eps',/eps
spawn,['epstopdf',file+'.eps'],/noshell
; Combine the figures                                                                                                                              
pdffiles = dir+'/plots/'+field+'_'+['input','recovered','completeness']+'.pdf'
spawn,'gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile='+dir+'/plots/'+field+'_complete.pdf '+strjoin(pdffiles,' ')


stop

end
