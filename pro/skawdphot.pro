;=====================================================================
;
;   SKAWDPHOT.PRO = Spanky's Kick-Ass Widget Driven Photometry
;
;   Version 1.4  12 July 1999
;
;   This is an updated version of SRM's photom2.for code.  It takes
;   data from either fillphot2.cl (SRM) or PHOTINFO.PRO (CP update)
;   and performs a fit of the observed photometry to the standard
;   photometry.  The user can interactively improve the fit using
;   IDL widget driven options, such as star deletion, eliminating terms,
;   etc.
;
;   CALLING SEQUENCE:  skawdphot,ddfile,sameair=sameair,samezero=samezero
;
;   ddfile = string defined as the name of the DATA.DAT file of photometric
;            data (as produced by fillphot2.cl or PHOTINFO.PRO
;   sameair = keyword which determines whether to solve for 1 airmass 
;             term for all nights, or for a different airmass term for
;             each night.  (/sameair=same for all, sameair=-1 different)
;   samezero = keyword which determines whether to solve for 1 zero
;              point term for all nights, or for a diff. zero for each
;              night.  (/samezero=same for all, samezero=-1 different)
;
;   History:  1989:  photom2.for written by Majewski
;             7/21/97:  rewritten and tested as skawdphot.pro by CP
;             9/29/97:  added if/then/else statement to take into
;                       account frames with only 1 star when calculating
;                       mean residual for that frame
;             12/12/97: changed definition of nights to take into account
;                       runs with nonsequential nights and/or that don't
;                       start with night 1
;             5/15/98:  changed mean vs. frame plot to show weighted 
;                       mean, not straight mean
;             5/18/98:  tweaked some small things (fixed plot range in 
;                       mean vs. frame plot, eliminated a variable that
;                       was never used)
;             7/7/99:   beginning to implement some changes suggested by
;                       Steve and Jamie, added delete all stars per frame,
;                       it now remembers which frames were selected previously
;                       for the offset subtraction, and there are now sliders
;                       allowing one to zoom in on the residual plots.
;
;=====================================================================

Pro  CurrentFit
;=====================================================================
;   
;    CURRENTFIT.PRO
;
;    Should print the current fit (matrix a) and std error (sigma) of
;    the current fit to the message window.
;
;=====================================================================

common skawdcommon,  r,w,winvert,a,sigsq,x,starname,realmag,realcol,secz,$
                     instmag,night,num,num2,UT,wt,vconst,Asigma,term,numterms,$
                     endframe,frame,mstars,m,i,j,k,Sv,whichterms,currentx,$
		     goodterms,residj,xstar,flag,applyme,samestar,$
                     invstat,solutionsig,same,momentrun,momentnight,sdrun,$
		     sdnight,numframes,momentframe,sdframe,mnights,frameindex,$
		     meanopt,subnights,meansuboff,sdsuboff,subframes,kill,$
                     numsubnights,numsubframes,framenames,numsamestars,$
		     sortstar,uniqstar,wtflag,plotflag,nowplaying,deathstar,$
                     tagged,printme,fpernight,picdiff,diff,frnight,outfile,$
                     same2,goodframes

common widgecommon,  base,options,main,plotme,togselect,togdone,vsselect,vsdone,$
                     offselect,offdone,ncselect,ncdone,fcselect,fcdone,fcopts,$
                     pfselect,pfdone,messages,pfopts,pssselect,pssdone,pssopts,$
                     pstar,pframe,dsdone,dsrestore,dsbase,showbase,smsigma,smdone,$
                     smdeleted,smclick,dssigma,killframe,fc2select,fc2done,dsframe,$
                     yholdmin,yholdmax,vsmax,vsmin
                     


;
;  The first iteration of CurrentFit (which is called when the program starts up)
;  Prints out some basic info, #stars, #nights, etc.  flag is a variable which is
;  initialized to 1 when the program is first called, but is set to 0 afterwards.
;
if (flag eq 1) then begin 		     
   if (same2 eq 1) then begin
      widget_control, messages, /append, set_value=$
      'This run assumes same zero point term and same airmass term for all nights'
   endif
   if (same2 ne 1 and same eq 1) then begin
      widget_control, messages, /append, set_value=$
      'This run solves for different zero point term but assumes same airmass term for all nights'
   endif
   if (same2 ne 1 and same ne 1) then begin
      widget_control, messages, /append, set_value=$
      'This run solves for different zero point and airmass term for each night'
      widget_control, messages, /append, set_value=$
      'WARNING: Solution will be poor for nights with little airmass coverage'
   endif
   firststring='Calculation proceeding with '+strcompress(string(j))+' stars from'$
   +strcompress(string(m))+' nights.'
   widget_control, messages, /append, set_value=firststring
   secondstring=strarr(m)
   for nak=0,m-1 do begin
      secondstring(nak)='Night '+strcompress(string(mnights(nak)))+' has'+$
      strcompress(string(mstars(nak)))+' stars'
;      widget_control, messages, /append, set_value=secondstring(nak)
   endfor
widget_control,messages,/append,set_value=secondstring
widget_control,messages,/append,set_value=' '
endif
flag=0
;
;  In PhotAlg the matrix inversion is done and the status of the inversion is
;  stored in a variable called invstat.  This is 0 for a successful inversion
;  and 1 or 2 if there is a problem.  The following prints out the status of
;  the most recent matrix inversion.
;
if (invstat eq 0) then begin
   widget_control, messages, /append, set_value='Inversion completed successfully'
endif 
if (invstat eq 1) then begin
   widget_control, messages, /append, set_value=$
   'Singular matrix, inversion invalid, solution unsuccessful'
endif 
if (invstat eq 2) then begin
      widget_control, messages, /append, set_value=$
      'WARNING: small pivot element used in gaussian elimination, solution inaccurate'
endif
;
;   Below we format and print out the values for the current fit.
;
solutionsig=strcompress(string(format='(E10.3)',solutionsig))
currentsig='Error of Solution, SIGMA = ' + solutionsig
widget_control, messages, /append, set_value=currentsig
widget_control, messages, /append, set_value='Coefficients and errors:'
solution=string(format='(F10.4)',a)
errorsol=string(format='(E10.3)',Asigma)
for zzz=0,i-1 do begin
   holder=strcompress(string(format='(I2)',zzz+1))
   currentval='a('+holder+') = '+solution(zzz)+' sigma = '+errorsol(zzz)$
      +'  '+term(whichterms(zzz))
   widget_control, messages, /append, set_value=currentval
endfor
widget_control,messages,/append,set_value=' '
;
;
end


Pro TermTogglerEvent_0,ev
;======================================================================
;
;   TERMTOGGLEREVENT_0.PRO
;
;   This is the event handler which accepts input from the nonexclusive
;   menu of terms and sets up the input matrix for PhotAlg based on
;   the terms selected by the user.  This particular code will be used
;   only if the user is using different airmass terms for all nights
;   and different zero points for all nights.
;
;======================================================================

common skawdcommon
common widgecommon

widget_control,ev.id,get_uvalue=widgetname

;
;  applyme is a 3-element array which holds the user's choice of which
;  optional terms to include out of color, airmass*color, and color^2.
;
applyme=intarr(3)

;  
;  When 'done' is selected, the current values of the 3 choices are
;  fixed, and the number of terms and which terms used in the solution
;  are updated.  Finally PhotAlg is called to update the solution using
;  only those terms chosen.
;
case widgetname of
   'togselect' : begin
                 end
     
    'togdone'  :  begin
                     widget_control,togselect,get_value=applyme
                     widget_control,ev.top,/destroy
                     whichhold=where(applyme eq 1)
                     numhold=n_elements(whichhold)
                     if (whichhold(0) eq -1) then begin
                        whichterms=indgen(2*m)
                     endif else begin
                        whichterms=intarr((2*m)+numhold)
                        whichterms(0:((2*m)-1))=indgen((2*m))
                        whichterms((2*m):(((2*m)-1)+numhold))=whichhold(0:(numhold-1))+(2*m)
                     endelse
 		     numterms=n_elements(whichterms)
		     goodterms=applyme
                     widget_control,main,/hourglass
                     PhotAlg
                  end
endcase

end


Pro TermTogglerEvent_1,ev
;======================================================================
;
;   TERMTOGGLEREVENT_1.PRO
;
;   This is the event handler which accepts input from the nonexclusive
;   menu of terms and sets up the input matrix for PhotAlg based on
;   the terms selected by the user.  This particular code will be used
;   only if the user is keeping the same airmass term for all nights,
;   but using different zero points for all nights
;
;======================================================================

common skawdcommon
common widgecommon

widget_control,ev.id,get_uvalue=widgetname

;  
;  see comments in _0 above
;
applyme=intarr(3)

case widgetname of
   'togselect' : begin
                 end
     
    'togdone'  :  begin
                     widget_control,togselect,get_value=applyme
                     widget_control,ev.top,/destroy
                     whichhold=where(applyme eq 1)
                     numhold=n_elements(whichhold)
                     if (whichhold(0) eq -1) then begin
                        whichterms=indgen(m+1)
                     endif else begin
                        whichterms=intarr(m+1+numhold)
                        whichterms(0:m)=indgen(m+1)
                        whichterms((m+1):(m+numhold))=whichhold(0:(numhold-1))+m+1
                     endelse
 		     numterms=n_elements(whichterms)
		     goodterms=applyme
                     widget_control,main,/hourglass
                     PhotAlg
                  end
endcase

end


Pro TermTogglerEvent_2,ev
;======================================================================
;
;   TERMTOGGLEREVENT_2.PRO
;
;   This is the event handler which accepts input from the nonexclusive
;   menu of terms and sets up the input matrix for PhotAlg based on
;   the terms selected by the user.  This particular code will be used
;   only if the user is keeping the same airmass term for all nights
;   and the same zero point for all nights
;
;======================================================================

common skawdcommon
common widgecommon

widget_control,ev.id,get_uvalue=widgetname

;
;  see comments in _0 above.
;
applyme=intarr(3)

case widgetname of
   'togselect' : begin
                 end
     
    'togdone'  :  begin
                     widget_control,togselect,get_value=applyme
                     widget_control,ev.top,/destroy
                     whichhold=where(applyme eq 1)
                     numhold=n_elements(whichhold)
                     if (whichhold(0) eq -1) then begin
                        whichterms=indgen(2)
                     endif else begin
                        whichterms=intarr(2+numhold)
                        whichterms(0:1)=indgen(2)
                        whichterms(2:(1+numhold))=whichhold(0:(numhold-1))+2
                     endelse
 		     numterms=n_elements(whichterms)
		     goodterms=applyme
                     widget_control,main,/hourglass
                     PhotAlg
                  end
endcase

end


Pro  TermToggler
;======================================================================
;
;   TERMTOGGLER.PRO
;
;   This brings up a menu widget which allows you to select
;   which of the optional terms of the solution you would like to include
;   in your current solution.  It then sources the PhotAlg code to 
;   calculate a new solution using the terms currently selected.
;
;======================================================================

common skawdcommon
common widgecommon

;
;  the if statements allows only 1 instance of TermToggler at a time.
;
if XRegistered('TermToggler') then return
;
;  Build a nonexclusive menu widget which allows selection of the optional
;  terms in the solution.
;
togglebase = widget_base(title='Select Terms',/column,group_leader=base)
togopts=['Color','Airmass * Color','Color^2']
togselect = cw_bgroup(togglebase,togopts,/column,/nonexclusive,$
                      uvalue='togselect',set_value=goodterms)
togdone=widget_button(togglebase,uvalue='togdone',value='Apply Changes')
;
widget_control,togglebase,/realize
;
;   How the solution is altered depends on the solution type the user selects
;   initially (same zero/air, diff zero/same air, diff zero/diff air).  This
;   is done in the Event Handler itself, so which one called depends on the
;   values of same and same2.
;
if (same2 eq 1) then begin
   xmanager,'TermToggler',togglebase,event_handler='TermTogglerEvent_2'
endif
if (same2 ne 1 and same eq 1) then begin
   xmanager,'TermToggler',togglebase,event_handler='TermTogglerEvent_1'
endif 
if (same2 ne 1 and same ne 1) then begin
   xmanager,'TermToggler',togglebase,event_handler='TermTogglerEvent_0'
endif

end


Pro PlotFrameEvent,ev
;======================================================================
;
;    PLOTFRAMEEVENT.PRO
;  
;    Event Handler for PlotFrame, exclusive menu of frame names for 
;    plotting.
;
;======================================================================
common skawdcommon
common widgecommon

widget_control,ev.id,get_uvalue=widgetname
;
;  nowplaying tells StarFind which type of data coordinates to use.
;  (i.e. for each different type of plot, there are different values plotted,
;  and therefore different coordinate systems apply to each).
;
nowplaying=4
;
;   When 'done' is selected, the stars for the frame chosen in the menu
;   are plotted in the draw widget.  Any deleted stars are shown with
;   a different symbol (done as oplot to change a + to a *)
;
case widgetname of
   'pfselect' :  begin
		 end

   'pfdone'   :  begin
                 widget_control,pfselect,get_value=pframe
                 widget_control,ev.top,/destroy
                 temp=where(frame eq frame(endframe(pframe)))
                 dead=where(wtflag(temp) eq 0)
                 tempname=pfopts(pframe)
                 temptitle='Residuals vs Star id for '+tempname
                 x1=min(samestar(temp))-1
                 x2=max(samestar(temp))+1
                 y1=(min(residj(temp))-(abs(0.25*(min(residj(temp))))))/yholdmin
                 y2=(max(residj(temp))+(abs(0.25*(max(residj(temp))))))/yholdmax
                 plot, samestar[temp], residj[temp], psym=1, xtitle='Star ID',$
                 ytitle='Residuals',title=temptitle,yrange=[y1,y2],xrange=[x1,x2]
                 if (dead(0) ne -1) then begin
                    oplot,samestar(temp(dead)),residj(temp(dead)), psym=7
                 endif
                 end
endcase

end

               
Pro PlotFrame
;======================================================================
;
;    PLOTFRAME.PRO
;
;    Here we have another widget menu, this time when you choose the
;    option to plot frame & star id in PlotResid, this asks you which
;    frame to plot against.
;
;======================================================================

common skawdcommon
common widgecommon
;
;  allow only one plotframe menu at a time with XRegistered.
;
if XRegistered('PlotFrame') then return
;
;   Build the menu widget of framenames using the framenames array 
;   initialized in the main program (SkawdPhot.pro).
;
pfbase = widget_base(title='Which Frame',/column,group_leader=vsbase)
pfopts = framenames
pfselect = cw_bgroup(pfbase,pfopts,/column,/exclusive,/scroll,$
                      x_scroll_size=120,y_scroll_size=150,$
                      uvalue='pfselect',set_value=0)
pfdone = widget_button(pfbase,uvalue='pfdone',value='All Done')
;
widget_control,pfbase,/realize
;
xmanager,'PlotFrame',pfbase,event_handler='PlotFrameEvent'

end


Pro PlotSameStarEvent,ev
;======================================================================
;
;    PLOTSAMESTAREVENT.PRO
;  
;    Event Handler for PlotSameStar, exclusive menu of star names for 
;    plotting.
;
;======================================================================
common skawdcommon
common widgecommon

widget_control,ev.id,get_uvalue=widgetname
;
nowplaying=5
;
;   When 'done' is selected, all stars matching the StarID chosen by
;   the user in the pssselect menu are plotted against there residual
;   in the draw widget, overplotted are deleted stars with a diff. symbol.
;
case widgetname of
   'pssselect' :  begin
        	  end

   'pssdone'  :  begin
                 widget_control,pssselect,get_value=pstar
                 widget_control,ev.top,/destroy
                 me=samestar(sortstar(uniqstar(pstar)))
                 temp=where(samestar eq me)
                 dead=where(wtflag(temp) eq 0)
                 tempname=pssopts(pstar)
                 temptitle='Residuals vs Star id for '+tempname
                 x1=me-1
                 x2=me+1
                 y1=(min(residj(temp))-(abs(0.25*(min(residj(temp))))))/yholdmin
                 y2=(max(residj(temp))+(abs(0.25*(max(residj(temp))))))/yholdmax
                 plot, samestar(temp), residj(temp), psym=1, xtitle='Star ID',$
                 ytitle='Residuals',title=temptitle,xrange=[x1,x2],yrange=[y1,y2]
                 if (dead(0) ne -1) then begin
                    oplot,samestar(temp(dead)),residj(temp(dead)), psym=7
                 endif
                 end
endcase

end
               

Pro PlotSameStar
;======================================================================
;
;    PLOTSAMESTAR.PRO
;
;    Here we have another widget menu, this time when you choose the
;    option to plot same star id in PlotResid, this asks you which
;    star to plot against.
;
;======================================================================

common skawdcommon
common widgecommon
;
;  only one PSS menu at a time, please.  (using XRegistered)
;
if XRegistered('PlotSameStar') then return
;
;  Build widget menu using by adding 'Star' to list of StarIDs initialized
;  in main pro (SkawdPhot.pro).
;
pssbase = widget_base(title='Which Star',/column,group_leader=vsbase)
pssopts = strarr(numsamestars)
for abc=0,numsamestars-1 do begin
   pssopts(abc)='Star '+strcompress(samestar(sortstar(uniqstar(abc))))
endfor
pssselect = cw_bgroup(pssbase,pssopts,/column,/exclusive,/scroll,$
                      x_scroll_size=120,y_scroll_size=150,$
                      uvalue='pssselect',set_value=0)
pssdone = widget_button(pssbase,uvalue='pssdone',value='All Done')

widget_control,pssbase,/realize

xmanager,'PlotSameStar',pssbase,event_handler='PlotSameStarEvent'

end


Pro PlotResidEvent,ev
;======================================================================
;
;   PLOTRESIDEVENT.PRO
;
;   This is the event handler which accepts input from the exclusive
;   menu of possible x-axes to plot the residuals versus.  
;
;======================================================================

common skawdcommon
common widgecommon
;
;  For all plots, the yrange is determined by adding/subtracting
;  1/4 of the max/min to the max/min of the range.
;
; y1=min(residj)-(abs(0.25*(min(residj))))/yholdmin
; y2=max(residj)+(abs(0.25*(max(residj))))/yholdmax
;
yholdmax=1 & yholdmin=1
;
;  we want to overplot with a new symbol all stars that have been deleted
;  dead is the index of these stars.
;
dead=where(wtflag eq 0)
;
widget_control,ev.id,get_uvalue=widgetname
;
;   When 'done' is selected, either the chosen option for the x-axis
;   is plotted vs. residuals, or in the case of frame or samestar, the
;   menu of these options is pulled up for input.
;
case widgetname of
   'vsmax'    : begin
                end

   'vsmin'    : begin
                end

   'vsselect' : begin
                end
     
    'vsdone'  :  begin
                     widget_control,vsmax,get_value=yholdmax
                     widget_control,vsmin,get_value=yholdmin
                     y1=(min(residj)-(abs(0.25*(min(residj)))))/yholdmin
                     y2=(max(residj)+(abs(0.25*(max(residj)))))/yholdmax
                     widget_control,vsselect,get_value=ordinate
                     widget_control,plotme,get_value=index
                     wset,index
                     widget_control,ev.top,/destroy
                     case ordinate of
                        0 : begin
                               nowplaying=0
                               plot, secz, residj, psym=1, xtitle='Airmass',$
                               ytitle='Residuals',yrange=[y1,y2]
                               if (dead(0) ne -1) then begin
                                  oplot,secz(dead),residj(dead),psym=7
                               endif
                            end
			1 : begin
                               nowplaying=1
                               plot, realcol, residj, psym=1, xtitle='Color',$
                               ytitle='Residuals',yrange=[y1,y2]
                               if (dead(0) ne -1) then begin
                                  oplot,realcol(dead),residj(dead),psym=7
                               endif
                            end
 			2 : begin
                               nowplaying=2
                               plot, realmag, residj, psym=1, xtitle='Real Magnitude',$
                               ytitle='Residuals',yrange=[y1,y2]
                               if (dead(0) ne -1) then begin
                                  oplot,realmag(dead),residj(dead),psym=7
                               endif
                            end
			3 : begin
                               nowplaying=3
                               plot, samestar, residj, psym=1, xtitle='Star ID',$
                               ytitle='Residuals',yrange=[y1,y2]
                               if (dead(0) ne -1) then begin
                                  oplot,samestar(dead),residj(dead),psym=7
                               endif
                            end
                        4 : PlotFrame
                        5 : PlotSameStar
        		6 : begin
                                nowplaying=-1
     				y3=(min(momentnight)-(abs(0.25*(min(momentnight)))))/yholdmin
				y4=(max(momentnight)+(abs(0.25*(max(momentnight)))))/yholdmax
                                x1=min(mnights)-1
				x2=max(mnights)+1
 				plot, mnights, momentnight, psym=1, xtitle='Night',$
			        ytitle='Mean Residual',yrange=[y3,y4],xrange=[x1,x2]
			    end
			7 : begin
                               nowplaying=-1
                              ; y5=min(momentframe)-(abs(0.25*(min(momentframe))))
                              ; y6=max(momentframe)+(abs(0.25*(max(momentframe))))
                              ; plot, frameindex, momentframe, psym=1, xtitle='Frame #',$
                                y5=(min(picdiff)-(abs(0.25*(min(picdiff)))))/yholdmin
                                y6=(max(picdiff)+(abs(0.25*(max(picdiff)))))/yholdmax
                                plot, frameindex, picdiff, psym=1, xtitle='Frame #',$
			        ytitle='Mean Residual',yrange=[y5,y6]
                            end
                     endcase
                  end

endcase

end



Pro PlotResid
;======================================================================
;
;   PLOTRESID.PRO
;
;   This brings up yet another menu widget that asks you to choose what
;   you would like to plot the residuals against.  Then it plots the
;   residual of each star vs. what you chose in the draw widget.
;
;=====================================================================

common skawdcommon
common widgecommon
;
;  Allow only 1 plotresid menu at a time with XRegistered.
;
if XRegistered('PlotResid') then return
;
;  If you try and choose the 'click on star' option in delete stars
;  or showme, it tests to see if there is a plot currently displayed
;  or not.  This is done with the plotflag flag variable (0=yes there
;  is a plot, 1=no plot)
;
plotflag=0
;
;  Build menu widget of possible x-axes for plotting residuals vs.
;
vsbase = widget_base(title='Which X-Axis',/column,group_leader=base,xsize=250)
vsmax  = widget_slider(vsbase,uvalue='vsmax',value=yholdmax,title=$
                       'Y-axis max value zoom factor',minimum=1)
vsmin  = widget_slider(vsbase,uvalue='vsmin',value=yholdmin,title=$
                       'Y-axis min value zoom factor',minimum=1)
vsopts=['Airmass','Color','Real Magnitude','Star ID','Frame & Star ID',$
        'Same Star ID','Mean vs Night','Mean vs Frame']
vsselect = cw_bgroup(vsbase,vsopts,/column,/exclusive,uvalue='vsselect',$
                     set_value=0)
vsdone = widget_button(vsbase,uvalue='vsdone',value='Plot Me')

widget_control,vsbase,/realize

xmanager,'PlotResid',vsbase,event_handler='PlotResidEvent'

end


PRO NightChooserEvent,ev
;=====================================================================
;
;   NIGHTCHOOSEREVENT.PRO
;
;   This is the event handler which accepts input from the nonexclusive
;   menu of nights and lets OffSet know which nights to include in the
;   mean residual calculation.
;
;=====================================================================

common skawdcommon
common widgecommon

widget_control,ev.id,get_uvalue=widgetname
;
;  nightme is an array which holds the user's selection of which
;  nights to include in the OffSet calculation...  1 for include 0 for not.
;
nightme=intarr(m)
;
;   When 'done' is selected, the index of nights is passed to OffSetSubtract
;   which only includes the selected ones in the zero point calculation.
;
case widgetname of
   'ncselect' : begin
		end

   'ncdone'   : begin
		  widget_control,ncselect,get_value=nightme
                  widget_control,ev.top,/destroy
                  subnights=where(nightme eq 1)
                  numsubnights=n_elements(subnights)
                  widget_control,main,/hourglass
                  OffSetSubtract
                end
endcase

end


PRO NightChooser
;=====================================================================
;
;   NIGHTCHOOSER.PRO
;
;   Guess what??  Another menu widget!!  This one allows you to choose
;   which nights you would like to include in the calculation of the
;   mean which is then subtracted in OffSet.
;
;=====================================================================

common skawdcommon
common widgecommon
;
;  only 1 NC menu allowed using XRegistered.
;
if XRegistered('NightChooser') then return
;
;   OffSetSubtract uses a different solution for each of the 3 choices
;   all nights, certain nights, certain frames.  Meanopt is the flag
;   variable which tells OSS which to use.
;
meanopt=1
;  
;   goodnights is an array of which nights have been selected in the
;   menu (1=selected, 0=not).  Below it is initialized to all 0s.
;
goodnights=intarr(m)
goodnights(*)=0
;
;   Build the widget menu of nights by adding the string 'Night'+the array
;   of night names initialized in SkawdPhot.pro.
;
ncbase = widget_base(title='Which Nights',/column,group_leader=offbase,xsize=150)
ncopts = strarr(m)
for cyndi=0,m-1 do begin
   ncopts(cyndi)='Night'+strcompress(mnights(cyndi))
endfor
ncselect = cw_bgroup(ncbase,ncopts,/column,/nonexclusive,uvalue='ncselect',$
                      set_value=goodnights)
ncdone = widget_button(ncbase,uvalue='ncdone',value='All Done')

widget_control,ncbase,/realize

xmanager,'NightChooser',ncbase,event_handler='NightChooserEvent'

end


PRO FrameChooserEvent,ev
;=====================================================================
;
;   FRAMECHOOSEREVENT.PRO
;
;   This is the event handler which accepts input from the nonexclusive
;   menu of frames and lets OffSet know which frames to include in the
;   mean residual calculation.
;
;=====================================================================

common skawdcommon
common widgecommon

widget_control,ev.id,get_uvalue=widgetname
;
;   frameme is an array which holds which frames to use in the offset
;   zeropoint calculation in OffSetSubtract.  For each frame 1=include,
;   0=not.
;
frameme=intarr(numframes)  &  frameme=goodframes
;
;   When done is selected, the current value of frameme is retrieved
;   and passed to OSS.
;
case widgetname of
   'fcselect' : begin
		end

   'fcdone'   : begin
		  widget_control,fcselect,get_value=frameme
                  widget_control,ev.top,/destroy
                  subframes=where(frameme eq 1)
                  goodframes=frameme
                  numsubframes=n_elements(subframes)
                  widget_control,main,/hourglass
                  OffSetSubtract
                end
endcase

end

PRO FrameChooserEvent2,ev
;=====================================================================
;
;   FRAMECHOOSEREVENT2.PRO
;
;   This is the event handler which accepts input from the exclusive
;   menu of frames and lets DeleteStars know which frame to delete 
;   all stars from
;
;=====================================================================

common skawdcommon
common widgecommon

widget_control,ev.id,get_uvalue=widgetname
;
;
case widgetname of
   'fc2select' : begin
               	 end

   'fc2done'   : begin
		  widget_control,fc2select,get_value=killframe
                  widget_control,ev.top,/destroy
                  temp=where(frame eq frame[endframe[killframe]])
                  wtflag[temp]=0
                  fc2string='All stars from '+frame[endframe[killframe]]+$
                           ' have been deleted.'
                  widget_control,messages,/append,set_value=fc2string
                end
endcase

end

PRO FrameChooser
;=====================================================================
;
;   FRAMECHOOSER.PRO
;
;   ARRRGGHHH!  Another menu widget!!  This one allows you to choose
;   which frames you would like to include in the calculation of the
;   mean which is then subtracted in OffSet.
;
;=====================================================================

common skawdcommon
common widgecommon
;
;  Only 1 FC menu allowd using XRegistered.
;
if XRegistered('FrameChooser') then return
;
;   meanopt is the flag for OffSetSubtract telling it which option to
;   use for the solution (Whole run, some nights, some frames)
;
meanopt=2
;
;   goodframes is an array which holds 1 for selected, 0 for not for
;   each frame in the menu.  It is initialized to 0 for all on startup
;   of FC.
;
; goodframes=intarr(numframes)
; goodframes(*)=0
;
;    Build the menu widget of framenames using the framenames variable
;    which is initialized in SkawdPhot.pro.
;
fcbase = widget_base(title='Which Frames',/column,group_leader=offbase)
fcopts = framenames
fcselect = cw_bgroup(fcbase,fcopts,/column,/nonexclusive,/scroll,$
                      x_scroll_size=120,y_scroll_size=150,$
                      uvalue='fcselect',set_value=goodframes)
fcdone = widget_button(fcbase,uvalue='fcdone',value='All Done')

widget_control,fcbase,/realize

xmanager,'FrameChooser',fcbase,event_handler='FrameChooserEvent'

end

PRO FrameChooser2
;=====================================================================
;
;   FRAMECHOOSER2.PRO
;
;   This time we have an exclusive menu of frames to choose which 
;   one to have all of its stars deleted.
;
;=====================================================================

common skawdcommon
common widgecommon
;
;  Only 1 FC menu allowd using XRegistered.
;
if XRegistered('FrameChooser2') then return

;
;    Build the menu widget of framenames using the framenames variable
;    which is initialized in SkawdPhot.pro.
;
fc2base = widget_base(title='Which Frame',/column,group_leader=dsbase)
fc2opts = framenames
fc2select = cw_bgroup(fc2base,fc2opts,/column,/exclusive,/scroll,$
                      x_scroll_size=120,y_scroll_size=150,$
                      uvalue='fc2select',set_value=0)
fc2done = widget_button(fc2base,uvalue='fc2done',value='All Done')

widget_control,fc2base,/realize

xmanager,'FrameChooser2',fc2base,event_handler='FrameChooserEvent2'

end

Pro OffSetEvent,ev
;======================================================================
;
;   OFFSETEVENT.PRO
;
;   This is the event handler which accepts input from the exclusive
;   menu of possible stars to use for the mean calculation and then
;   offset subtraction
;
;======================================================================

common skawdcommon
common widgecommon

widget_control,ev.id,get_uvalue=widgetname
;
;  When 'done' is selected, the meanopt flag is set based on which
;  choice was made by the user, and if Night or Frame is chosen, 
;  the NightChooser or FrameChooser menus are called to allow the
;  user to select which nights or frames to use in the offset zero
;  point calculation.  Then OffSetSubtract is called to do the calculation.
;
case widgetname of
   'offselect' : begin
                 end
     
    'offdone'  :  begin
                     widget_control,offselect,get_value=whichmean
                     case whichmean of
                        0 : begin
                              widget_control,ev.top,/destroy
		              meanopt=0
                              widget_control,main,/hourglass
                              OffSetSubtract
                            end
			1 : begin
                              NightChooser
                              widget_control,ev.top,/destroy
                            end
 			2 : begin
                              FrameChooser
                              widget_control,ev.top,/destroy
                            end
                     endcase
                  end
endcase

end


PRO  OffSet
;======================================================================
;
;  OFFSET.PRO
;
;  This brings up yes, **another** menu widget that asks you to choose
;  either the whole run, certain nights, or certain frames to use to
;  determine the mean residual value, which is then subtracted from the
;  instrumental magnitude of all stars.
;
;======================================================================

common skawdcommon
common widgecommon
;
;  Only one instance of the OS menu at a time using the XRegistered.
;
if XRegistered('OffSet') then return
;
;   Build menu widget of 3 options, all stars, certain nights, or certain
;   frames.  
;
offbase = widget_base(title='Which Mean',/column,group_leader=base)
offopts = ['All Stars','Certain Nights','Certain Frames']
offselect = cw_bgroup(offbase,offopts,/column,/exclusive,uvalue='offselect',$
                      set_value=0)
offdone = widget_button(offbase,uvalue='offdone',value='All Done')

widget_control,offbase,/realize

xmanager,'OffSet',offbase,event_handler='OffSetEvent'

end


Pro  OffSetSubtract
;======================================================================
;
;   OFFSETSUBTRACT.PRO
;
;   Once the user chooses which method to use in calculating the mean
;   residual, this program is called which performs the mean calculation
;   and subtraction.
;
;======================================================================

common skawdcommon
common widgecommon
;
;
if (meanopt eq 1) then begin
   totalframes=0
   for huey=0,numsubnights-1 do begin
      totalframes=totalframes+fpernight(subnights(huey))
   endfor
   avghold=0.00
   for ghi=0,numsubnights-1 do begin
      holder=where(frnight eq mnights(subnights(ghi)))
      avghold=avghold+total(picdiff(holder))
   endfor
   refresid=avghold/totalframes
   picdiff=picdiff-refresid
endif
if (meanopt eq 2) then begin
   avghold=(picdiff(subframes))
   refresid=(total(avghold)/numsubframes)
   picdiff=picdiff-refresid
endif
;
;   picdiff (calc. in PhotAlg) is the weighted average residual per frame.  If the user
;   chooses to use the whole run, this is subtracted from diff (v-V)
;   for all stars in each frame to improve the fit.  If the user chooses
;   the night or frame option, then a reference "zero point" avg residual
;   is determined from the frames specified.  The difference between
;   the wt. avg. residual for each frame - this zero point is then subtracted
;   from all stars in each frame in this case.
;
widget_control,messages,/append,set_value=$
   'CCD Offset Correction being applied:'
for piphi=0,numframes-1 do begin
   frameindices=where(frame eq frame(endframe(piphi)))
   diff(frameindices)=diff(frameindices)-picdiff(piphi)
   picdiffstr=string(format='(E12.3)',picdiff)
   submsg='Subtracting '+picdiffstr(piphi)+$
          ' from all stars in '+framenames(piphi)
   widget_control,messages,/append,set_value=submsg
endfor

widget_control,messages,/append,set_value=' '
;
;  After the subtraction is performed, PhotAlg is called to calculate the
;  current solution.
;
PhotAlg

end


PRO  DeleteStarsEvent,ev
;====================================================================
;
;   DELETESTARSEVENT.PRO
;
;   This is the event handler for delete stars.  It sets up the plotme
;   widget to accept button click events, and then sends the coords
;   of the click to StarFind to find the star which matches the click.
;   It sets the kill flag to 1, so when StarFind finds the star, it
;   has its flag set to 0.
;
;====================================================================

common skawdcommon
common widgecommon

widget_control,ev.id,get_uvalue=widgetname

case widgetname of
   'dssigma'  :  begin
                    threesig=3.*sdrun
                    diffmean=abs(residj-momentrun[0])
                    killme=where(diffmean gt threesig,tsct)
                    if (tsct ne 0) then wtflag[killme]=0
                    widget_control,messages,/append,set_value=$
                    'All stars with residual > 3 sigma have been deleted.'
                 end

   'dsframe'  :  begin
                    FrameChooser2
                 end

;
;  by clicking on dsrestore, the wtflag of each star is reset to 1, 
;  thereby making all stars again available when calculating the
;  solution.
;
   'dsrestore' :  begin
                    wtflag[*]=1
                    widget_control,messages,/append,set_value=$
                    'All deleted stars have been restored.'
                  end
;
;  by clicking on dsundel, it sets kill to 2, which means when StarFind
;  determines which star has been clicked on it sets wtflag to 1, which
;  restores the star to being available to the solution if it had previously
;  been deleted.
;
   'dsundel'   :  begin
                    kill=2
                    widget_control,messages,/append,set_value=$
                    'Click on any deleted star to undelete it'
                  end
;
;   dsdel allows you to return to deleting after you have been undeleting
;   (the default is deletion, so this button does not need to be clicked
;   before starting to delete.
;
   'dsdel'     :  begin
                    widget_control,messages,/append,set_value=$
                    'Click on any star to delete it'
                    kill=1
                  end
;
;  When 'done' is clicked, events generated by clicking on the draw windo
;  are no longer accepted.  Also, PhotAlg is called to recalculate the 
;  solution using the stars that have not been deleted.
;
   'dsdone'    :  begin
                    widget_control,ev.top,/destroy
                    kill=0
                    widget_control,plotme,draw_button_events=0
                    widget_control,messages,/append,set_value='Done Deleting'
                    widget_control,messages,/append,set_value=' '
                    widget_control,main,/hourglass
                    PhotAlg
                  end

   'plotme'    :  begin
                    if (ev.type eq 0) then begin
                       if (kill ne 2) then kill=1
                       coo=[ev.x,ev.y]
                       deathstar=convert_coord(coo,/device,/to_data)
                       StarFind
                    endif
                  end
endcase


end



PRO  DeleteStars
;=====================================================================
;
;   DELETESTARS.PRO
;
;   This program will allow the user to click on a star which is then
;   removed from the solution.  This sets up the widgets for DSEvent.
;   
;======================================================================

common skawdcommon
common widgecommon
;
;  Only 1 instance of DeleteStars allowed.
;
if XRegistered('DeleteStars') then return
;
;  Also, we don't want DeleteStars and ShowMe running at the same time,
;  since they give conflicting button events, so we'll use XRegistered
;  to prevent this, too.
;
if XRegistered('ShowMe') then return
;
widget_control,plotme,/draw_button_events
;
;  Default when began is to delete stars (kill=1), user can choose
;  to undelete by clicking on dsundel.
;
kill=1
;
;  If there is no current plot, calls PlotResid to allow user to
;  plot stars, which can then be deleted by clicking if necessary.
;
if (plotflag eq 1) then begin
   widget_control,messages,/append,set_value=$
   'Please Plot Residuals before attempting to Delete Stars'
   PlotResid
endif
;
;  build widget menu with the following options...
;
dsbase=widget_base(title='Deleting Stars',/column,group_leader=base)
dssigma=widget_button(dsbase,uvalue='dssigma',value=$
                      'Delete all stars with residuals > 3 sigma')
dsframe=widget_button(dsbase,uvalue='dsframe',value=$
                      'Delete all stars from one frame')
dsrestore=widget_button(dsbase,uvalue='dsrestore',value=$
                        'Restore All Deleted Stars')
dsundel=widget_button(dsbase,uvalue='dsundel',value=$
                      'Switch to Undeletion')
dsdel=widget_button(dsbase,uvalue='dsdel',value=$
                      'Switch to Deletion')
dsdone=widget_button(dsbase,uvalue='dsdone',value='Done Deleting')

widget_control,dsbase,/realize
widget_control,messages,/append,set_value='Click on any star to delete it'

xmanager,'DeleteStars',dsbase,event_handler='DeleteStarsEvent'
xmanager,'DeleteStars',plotme,event_handler='DeleteStarsEvent'

end

PRO  StarFind
;======================================================================
;
;   STARFIND.PRO
;
;   This program is called by either DeleteStars or ShowMe.  Those two
;   programs accept click events from the plotme widget and send the
;   coordinates here.  This code narrows the range of possible values
;   of the x-axis of the plot from which the click came from by finding
;   the minimum value of (x-click)-(all-x-values).  Then, it does the same
;   for the residual (y-axis) and finds where the possible x choices match
;   the y choice.  This algorithm is therefore more sensitive to erroneous
;   clicks in x, than in y.
;
;======================================================================

common skawdcommon
common widgecommon
;
;  You can't delete stars or find stars when no stars are plotted!
;  (nowplaying=-1 means mean vs. night or mean vs. frame is plotted.
;
if (nowplaying eq -1) then begin
   widget_control,messages,/append,set_value=$
   'You can not use this plot option for star deletion or star info.'
endif else begin
;
;  findit is the absolute value of the x-coordinate of the click event
;  minus the x-coordinate of all stars.
;
   if (nowplaying eq 0) then begin
      findit=abs(deathstar[0]-secz[*])
      findtemp=secz
   endif
   if (nowplaying eq 1) then begin
      findit=abs(deathstar[0]-realcol[*])
      findtemp=realcol
   endif
   if (nowplaying eq 2) then begin
      findit=abs(deathstar[0]-realmag[*])
      findtemp=realmag
   endif
   if (nowplaying eq 3) then begin
      findit=abs(deathstar[0]-samestar[*])
      findtemp=samestar
   endif
;
;  In the Frame & Same star plot options, we can narrow the range
;  of possible stars to find by clicking by only using what's displayed
;  (frame or star id) rather than the whole array of star values.
;
   if (nowplaying eq 4) then begin
      helper=where(frame eq frame[endframe[pframe]])
      findit=dblarr(j) & findit[*]=1.0d4
      findit[helper]=abs(deathstar[0]-samestar[helper])
      findtemp=samestar
   endif
   if (nowplaying eq 5) then begin
      helper=where(samestar eq samestar[sortstar[uniqstar[pstar]]])
      findit=dblarr(j) & findit[*]=1.0d4
      findit[helper]=abs(deathstar[0]-samestar[helper])
      findtemp=samestar
   endif
;
;  min1 is the minimum of the difference of all stars from the click pos'n.
;  Since this will probably be multi-valued (for instance many stars have the
;  same airmass) we then use the minimum of (yclick)-(residual) of all stars
;  to narrow down which star is exactly it.
;
   min1=where(findit eq min(findit))
   find=abs(deathstar[1]-residj[min1])
   min2=where(find eq min(find))
   tagged=min1[min2]
   killmsg='Star '+strcompress(string(samestar[tagged]))+' on Frame '+$
              strmid(strcompress(frame[tagged]),0,6)
   gotit=fix(wtflag[tagged[0]])
;
;   Print out the star has been deleted, has already been deleted (if you
;   try and delete a previously deleted star) has been restored (undelete 
;   option in DeleteStars) or has already been restored (if you try and
;   restore a previously restored star or restore a star that was never
;   deleted to begin with.)
;
   if (kill eq 1) then begin
      if (gotit eq 0) then begin
         killmsg=killmsg+' has already been deleted'
         widget_control,messages,/append,set_value=killmsg
      endif else begin
      killmsg=killmsg+' has been deleted'
      wtflag[tagged]=0
      widget_control,messages,/append,set_value=killmsg
      endelse
   endif 
   if (kill eq 2) then begin
      if (gotit eq 1) then begin
         killmsg=killmsg+' has not been deleted, or has already been restored'
         widget_control,messages,/append,set_value=killmsg
      endif else begin
      killmsg=killmsg+' has been restored'
      wtflag[tagged]=1
      widget_control,messages,/append,set_value=killmsg
      endelse
   endif
endelse

end


PRO  StarPrint
;====================================================================
;
;    STARPRINT.PRO
;
;    This is called by ShowMe.  It prints information on the star(s)
;    selected by the user to the messages widget.
;
;====================================================================

common skawdcommon
common widgecommon
;
;  This just determines for each star selected by ShowMe whether it
;  has or has not been deleted.  This is included below in the info
;  which is printed out.
;
if (printme[0] ne -1) then begin
delmsg=strarr(n_elements(printme))
for mno=0,n_elements(printme)-1 do begin
   if (wtflag[printme[mno]] eq 0) then delmsg[mno]='deleted' else $
        delmsg[mno]='not deleted'
endfor
endif
;
;  Below we print the info listed in starhead to the messages window for
;  each star selected.  If none fit the criterion (none >3sig, no deleted stars,
;  etc.) than this is printed.
;
starmsg=strarr(n_elements(printme))
starhead='   ID     Night   Residual     MDiff   Color     Air       Wt       Status'
if (printme[0] eq -1) then begin
   starmsg='No stars fit the criterion selected'
endif else begin
   for efg=0,n_elements(printme)-1 do begin
      starmsg[efg]=strcompress(starname[printme[efg]])+'   '+strcompress(string(night[printme[efg]]))+$
                   '     '+string(format='(F7.4)',residj[printme[efg]])+$
                   '     '+string(format='(F6.3)',diff[printme[efg]])+$
                   '  '+string(format='(F6.3)',realcol[printme[efg]])+$
                   '   '+string(format='(F6.3)',secz[printme[efg]])+$
                   '   '+string(format='(F8.5)',wt[printme[efg]])+$
                   '   '+strcompress(delmsg[efg])
   endfor
endelse

widget_control,messages,/append,set_value=starhead
widget_control,messages,/append,set_value=starmsg
widget_control,messages,/append,set_value=' '

end

PRO  ShowMeEvent,ev
;====================================================================
;
;   SHOWMEEVENT.PRO
;
;   This is the event handler for ShowMe.  it accepts input from a
;   list of button widgets.  Then it calls StarPrint to print out
;   info on the stars selected.
;
;====================================================================

common skawdcommon
common widgecommon

widget_control,ev.id,get_uvalue=widgetname
;
;   You can select stars that have residj >3 sig from the mean of
;   all stars, those which have been deleted, or by clicking.  Kill
;   is set to 0 for clicked on stars to prevent them from being
;   deleted by StarFind.
;
case widgetname of
   'smsigma' :  begin
                    widget_control,messages,/append,set_value=$
                    'All stars with residuals > 3 sigma:'
                    threesig=3.*sdrun
                    diffmean=abs(residj-momentrun[0])
                    printme=where(diffmean gt threesig)
                    StarPrint
                end

   'smdeleted' : begin
                    widget_control,messages,/append,set_value=$
                    'All stars which have been deleted from the solution:'
                    printme=where(wtflag eq 0)
                    StarPrint
                 end

   'smclick'   :  begin
                     widget_control,messages,/append,set_value=$
                     'The star you clicked on:'
                     widget_control,plotme,/draw_button_events
                     if (plotflag eq 1) then begin
                        widget_control,messages,/append,set_value=$
                        'Please Plot Residuals before attempting to Click on Stars'
                        PlotResid
                     endif
                     kill=0
                  end

   'plotme'   :  begin
                    kill=0
                    if (ev.type eq 0) then begin
                    coo=[ev.x,ev.y]
                    deathstar=convert_coord(coo,/device,/to_data)
                    StarFind
                       if (nowplaying ne -1) then begin
                          printme=tagged
                          StarPrint
                       endif
                    endif
                 end


   'smdone'    :  begin
                    widget_control,ev.top,/destroy
                    widget_control,messages,/append,set_value=' '
                    widget_control,plotme,draw_button_events=0
                  end

endcase

end


PRO  ShowMe
;======================================================================
;
;   SHOWME.PRO
;
;   This is a widget creation program which sets up buttons that allow
;   the user to select certain stars, which then have their info printed
;   to the messages widget.
;
;======================================================================

common skawdcommon
common widgecommon
;
;   Allow only 1 SM button list.
;
if XRegistered('ShowMe') then return
;   
;  To prevent button click conflicts, limit the use of DS at the same
;  time as SM, too.
;
if XRegistered('DeleteStars') then return
;
widget_control,plotme,/draw_button_events
;
;  build several button widgets with the following options...
;
showbase=widget_base(title='Star Information',/column,group_leader=base)
smsigma=widget_button(showbase,uvalue='smsigma',value=$
                        'Show all stars with resid > 3 sigma')
smdeleted=widget_button(showbase,uvalue='smdeleted',value=$
                        'Show all deleted stars')
smclick=widget_button(showbase,uvalue='smclick',value=$
                        'Click on a star')
smdone=widget_button(showbase,uvalue='smdone',value='Done')

widget_control,showbase,/realize

xmanager,'ShowMe',showbase,event_handler='ShowMeEvent'
xmanager,'ShowMe',plotme,event_handler='ShowMeEvent'

end
 

PRO  MainHandler, ev
;======================================================================
;
;   MAINHANDLER.pro
;
;   This is the event handler for the main widget which creates the
;   menu of options, plotting area, and message window.  
;
;======================================================================

common  skawdcommon
common widgecommon

widget_control,ev.id
;
;  The main window has a row of buttons on top.  Each one calls a different
;  program as listed below.  When quit is selected, all info printed to the
;  messages window is dumped to the file named ddfile+'.out'.
;
case ev.value of
     0 : TermToggler
     1 : CurrentFit
     2 : PlotResid
     3 : DeleteStars
     4 : OffSet
     5 : ShowMe
     6 : begin
          widget_control,messages,get_value=printfile
          widget_control,ev.top,/destroy
          get_lun,unit
          openw,unit,outfile
          for jkl=0,n_elements(printfile)-1 do begin
             printf,unit,printfile(jkl)
          endfor
          close,unit
          free_lun,unit
         end
endcase

end


PRO  PhotAlg
;======================================================================
;
;   PHOTALG.PRO
;
;   This program contains the algorithm for performing the matrix 
;   inversion solution of the transformation equation using the 
;   technique in Harris et al. 1981 (PASP 83, 507)
;
;======================================================================

common skawdcommon
common widgecommon
;
;   The initial fit is done using all of the terms of the solution.  
;   However, subsequent fits may toggle one or more terms on and off. 
;   Here we will define a new array which contains only those columns 
;   currently toggled on by the user.
;
;   x is the initialarray of terms and is never altered, so terms can be 
;   toggled on and off at will, and only currentx needs to be updated.
;
currentx=dblarr(numterms,j)
for yyy=0,numterms-1 do begin
   currentx[yyy,*]=x[whichterms[yyy],*]
endfor
i=numterms & k=i
;
;   The following set of nested loops fills the weighting array which
;   is inverted to give you the transformation equation constants.
;
;   The weighting array is the wt. of the star (calculated in photinfo.pro
;   or in fillphot2.cl) times a flag which is 1 if the star is ok, and 0 if
;   it is deleted (and hence removed from the solution) * the terms of the
;   observed quantities array.
;
w=dblarr(i,k)
wxx=dblarr(j)
sum=double(0.0)
for alpha=0,i-1 do begin
   for beta=0,k-1 do begin
      for gamma=0,j-1 do begin
         wxx[gamma]=wtflag[gamma]*wt[gamma]*currentx[alpha,gamma]*currentx[beta,gamma]
      endfor
   sum=total(wxx)
   w[alpha,beta]=sum
   sum=double(0.0)
   endfor
endfor
;
;  Now we need to fill the r-array, which is the array of observed quantities
;  to be fit.
;
r=dblarr(k)
wxvV=dblarr(j)
sum2=double(0.0)
for delta=0,k-1 do begin
   for epsilon=0,j-1 do begin
      wxvV[epsilon]=wtflag[epsilon]*wt[epsilon]*currentx[delta,epsilon]*$
                    diff[epsilon]
   endfor
   sum2=total(wxvV)
   r[delta]=sum2
   sum2=double(0.0)
endfor
;
;   Ok.  The first solution will be a=W^(-1)*r.  So let's do it.
;   invstat is a status variable which flags with 1 or 2 if the 
;   inversion is unsuccessful.
;
a=dblarr(k)
winvert=dblarr(i,k)
winvert=invert(w,invstat,/double)
a=winvert#r
;
;   Now we compute the residual per star using the just calculated fit, a.
;   xstar is the array of solution terms*the observed quantities, and should
;   therefore equal the standard magnitude.  The difference between the fit
;   and the standard is the residual.
;
xstar=dblarr(numterms,j)
if (same2 eq 1) then begin
   for jon=0,numterms-1 do begin
      xstar[jon,*]=a[jon]*currentx[jon,*]
   endfor
endif
if (same2 ne 1 and same eq 1) then begin
   for jon=0,m-1 do begin
      xstar[jon,*]=a[jon]*currentx[jon,*]
   endfor
   for lambda=0,(numterms-m-1) do begin
      xstar[(lambda+m),*]=a[lambda+m]*currentx[(lambda+m),*]
   endfor
endif
if (same2 ne 1 and same ne 1) then begin
   for jon=0,m-1 do begin
      xstar[jon,*]=a[jon]*currentx[jon,*]
      xstar[(jon+m),*]=a[jon+m]*currentx[(jon+m),*]
   endfor
   for lambda=0,((numterms)-(2*m)-1) do begin
      xstar[(lambda+(2*m)),*]=a[lambda+(2*m)]*currentx[(lambda+(2*m)),*]
   endfor
endif
;
;
Ajstar=dblarr(j)
for nu=0,numterms-1 do begin
   Ajstar[*]=Ajstar[*]+xstar[nu,*]
endfor
residj=diff
residj[*]=residj[*]-Ajstar[*]
;
;   Now we compute the standard error matrix defined as 
;   sigsquared=W^(-1)*sigsquared_v.  where ss_v=residuals of fit.
;
Sv=double(0.0)
vconst=double(0.0)
solutionsig=double(0.0)
wvVA=dblarr(j)
for zeta=0,j-1 do begin
   wvVA[zeta]=wtflag[zeta]*wt[zeta]*((diff[zeta]-Ajstar[zeta])^(2))
endfor
Sv=total(wvVA,/double)
;
;  vconst is the variance of the entire solution, so sqrt(vconst) is
;  the sigma of the solution.
;
vconst=Sv/(j-numterms)
solutionsig=sqrt(vconst)
sigsq=dblarr(numterms)
Asigma=dblarr(numterms)
for omicron=0,numterms-1 do begin
   sigsq[omicron]=winvert[omicron,omicron]*vconst
endfor
Asigma=sqrt(sigsq)
;
;    We want to compute the mean residual for each frame, each night, and
;    the entire run for use in plotting residuals and deleting >3sig stars.
;
momentrun=moment(residj,/double,sdev=sdrun)
momentnight=dblarr(m)
sdnight=dblarr(m)
for iota=0,m-1 do begin
    nighthold=where(night eq mnights[iota])
    yarhold=(residj[nighthold]*wt[nighthold]*wtflag[nighthold])
    momenthold=moment(yarhold,/double,sdev=sdhold)
    momentnight[iota]=momenthold[0]
    sdnight[iota]=sdhold
endfor
momentframe=dblarr(numframes)
sdframe=dblarr(numframes)
for theta=0,numframes-1 do begin
    frameindices=where(frame eq frame[endframe[theta]])
    if (n_elements(frameindices) eq 1) then begin
       momentframe[theta]=residj[frameindices[0]]
       sdframe[theta]=99.99
    endif else begin
       momentholder=moment(residj[frameindices],/double,sdev=sdholder)
       momentframe[theta]=momentholder[0]
       sdframe[theta]=sdholder
    endelse
endfor
;
;   We want to calculate a weighted average residual per frame for use in
;   The offset correction
;
picdiff=dblarr(numframes)
for tau=0,numframes-1 do begin
   pdhold=where(frame eq frame[endframe[tau]])
   pd3hold=total(wt[pdhold]*wtflag[pdhold]*residj[pdhold])
   pd4hold=total(wt[pdhold]*wtflag[pdhold])
   if (pd3hold eq 0.0 or pd4hold eq 0.0) then begin
      picdiff[tau]=0.0
   endif else begin
      picdiff[tau]=pd3hold/pd4hold
   endelse
endfor
;
;    Now that the solution is finished, the event handler should be called
;    so that the user can begin to improve the solution.  
;    If this is the first iteration of PhotAlg, the initial fit is printed 
;    with CurrentFit.  Subsequent iterations require the user to click on 
;    CurrentFit to show the current solution.
;
if (flag eq 1) then begin
   CurrentFit
   xmanager,'PhotAlg',base,event_handler='MainHandler'
endif

end


PRO  SkawdPhot,ddfile,SAMEAIR=sameair,SAMEZERO=samezero
;=====================================================================
;
;   SKAWDPHOT.PRO = Spanky's Kick-Ass Widget Driven Photometry
;
;   This is the main program.  It reads in the data from the datafile
;   initializes many variables, builds the main widget window, and then
;   sources PhotAlg to determine the first solution.
;
;=====================================================================

common  skawdcommon
common  widgecommon
;
;   All messages printed to the text widget will be printed out at the
;   conclusion of the program.  So the first thing we do is create a
;   filename to store this in by adding .OUT to the name of the ddfile.
;
outfile=ddfile+'.out'
;
;   For multi-night solutions, you have the option of using the same
;   airmass term for all nights, or solving for the airmass term 
;   separately for each night.  This option is specified as the keyword
;   sameair in the call to skawdphot.  The default is to use the same
;   airmass term for all nights.
;
if (keyword_set(sameair) eq 0) then sameair=1
same=sameair
;
;   For multi-night solutions, you also have the option of using the
;   same zero point for all nights, or solving the zero point separately
;   for each night.  This option is specified as the keyword 
;   samezero in the call to skawdphot.  The default is to use different
;   zeros for all nights.  Also, if you choose samezero, by default you
;   want sameair, so sameair is set to 1
;
if (keyword_set(samezero) eq 0) then samezero=-1
same2=samezero
if (same2 eq 1) then same=1
;
;   As usual, we start by reading in the data from an ascii file.
;   the name of the file is included as a call to the program, and
;   is stored in the string 'ddfile'.
;
;
readcol,ddfile,starname,realmag,realcol,secz,instmag,night,num,UT,wt,$
    format='a,d,d,d,d,i,i,a,d'
peach=uniq(night)
mnights=night[peach]  
m=n_elements(mnights)
j=n_elements(starname)
num2=indgen(j)
;
;   When we delete stars, we want this to be reversible.  So let's introduce
;   a flag variable which is initially 1 for all stars.  It will be multiplied
;   by the weight, so a deleted star will have its flag = 0, and will therefore
;   drop out of the calculation.
wtflag=intarr(j)
wtflag[*]=1
;
;   Next we'll separate the frame names from the star numbers so we can
;   later plot things vs. stars from individual frames.  We will also determine
;   where each frame ends so we can separate the stars in each frame.
;   Let's also find stars of the same name so we can use that, too.
;
dashpos=strpos(starname,'-')
frame=strarr(j)
samestar=strarr(j)
for phi=0,j-1 do begin
    frame[phi]=strmid(starname[phi],0,dashpos[phi])
    samestar[phi]=strmid(starname[phi],(dashpos[phi]+1),4)
endfor
endframe=uniq(frame)
numframes=n_elements(endframe)
frameindex=indgen(numframes)
frameindex=frameindex+1
framenames=strarr(numframes)
frnight=intarr(numframes)
for bcd=0,numframes-1 do begin
   framenames[bcd]='Frame '+strcompress(frame[endframe[bcd]])
   frnight[bcd]=night[endframe[bcd]]
endfor
sortstar=sort(samestar)
uniqstar=uniq(samestar[sortstar])
numsamestars=n_elements(uniqstar)
;
;   Here we determine the number of stars in the datafile, and define
;   the x vector (zeropoint,airmass,color,airmass*color,color squared)
;   for each of the j stars.  
;
;   for multiple night, simultaneous solutions, you have the choice
;   of using the same airmass term, or different airmass terms for each
;   night.  The huge if then else structure fills x depending on how
;   the keyword sameair is set.  If sameair is true, it fills x using
;   a different zeropoint for each night, if false it fills x using a
;   different zero and air term for each night.  Also, you have the
;   choice of multiple zero points or not.
;
mstars=intarr(m)
fpernight=intarr(m)
for iii=0,m-1 do begin
   jj=where(night eq mnights[iii])
   mstars[iii]=n_elements(jj)
   jjhold=uniq(frame[jj])
   fpernight[iii]=n_elements(jjhold)
endfor
if (same2 eq 1) then begin
   numterms=5
   x=dblarr(5,j)
   for lou=0,j-1 do begin
      x[0,lou]=1
      x[1,lou]=secz[lou]
      x[2,lou]=realcol[lou]
      x[3,lou]=(secz[lou]*realcol[lou])
      x[4,lou]=((realcol[lou])^(2))
   endfor
;
;
   term=strarr(numterms)
   term[0]='zero point'
   term[1]='airmass'
   term[2]='color'
   term[3]='airmass*color'
   term[4]='color^2'
endif
if (same2 ne 1 and same eq 1) then begin
   numterms=m+4
   x=dblarr(numterms,j)
   for lou=0,j-1 do begin
      for joe=0,m-1 do begin
         if (night[lou] eq mnights[joe]) then begin
            x[joe,lou]=1
         endif else begin
            x[joe,lou]=0
         endelse
      endfor
      x[m,lou]=secz[lou]
      x[(m+1),lou]=realcol[lou]
      x[(m+2),lou]=(secz[lou]*realcol[lou])
      x[(m+3),lou]=((realcol[lou])^(2))
   endfor
;
;
   term=strarr(numterms)
   for roy=0,m-1 do begin
      term[roy]='Night '+strcompress(string(mnights[roy]))+' zero point'
   endfor
   term[m]='airmass'
   term[m+1]='color'
   term[m+2]='airmass*color'
   term[m+3]='color^2'
;
;
endif
if (same2 ne 1 and same ne 1) then begin
   numterms=(2*m)+3
   x=dblarr(numterms,j)
   for lou=0,j-1 do begin
      for joe=0,m-1 do begin
         if (night[lou] eq mnights[joe]) then begin
            x[joe,lou]=1
            x[(joe+m),lou]=secz[lou]
         endif else begin
            x[joe,lou]=0
            x[(joe+m),lou]=0
         endelse
      endfor
      x[(2*m),lou]=realcol[lou]
      x[((2*m)+1),lou]=(secz[lou]*realcol[lou])
      x[((2*m)+2),lou]=((realcol[lou])^(2))
   endfor
;
;
   term=strarr(numterms)
   for roy=0,m-1 do begin
      term[roy]='Night '+strcompress(string(mnights[roy]))+' zero point'
      term[roy+m]='Night '+strcompress(string(mnights[roy]))+' airmass'
   endfor
   term[2*m]='color'
   term[(2*m)+1]='airmass*color'
   term[(2*m)+2]='color^2'
endif
;
;  The initial fit performed upon startup of the program will include
;  all terms.  So we initialize the number of terms, and the i and k
;  indices to reflect this.  
;
flag=1
plotflag=1
whichterms=indgen(numterms)
goodterms=intarr(3) & goodterms[*]=1
;
;   Jamie and Steve want the program to remember which frames/nights you
;   selected in the offset correction so that each time you don't have to
;   redo the choice.  So I need to initialize these arrays, as well.
;   I'm also initializing the y-axis residual plot zoom factor.
;
goodframes=intarr(numframes) & goodframes[*]=1
yholdmin=1 & yholdmax=1
;
;   We will want to compute the residuals by subtracting the fitted magnitude
;   from the known standard magnitude.  
;
diff=dblarr(j)
diff=instmag-realmag
residj=dblarr(j)
;
;   
;===========================================================================
;  WIDGET CREATION
;
;  Below the menu, plot area, and message area widgets are created.
;
;
base = widget_base(title='Interactive Photometric Calibration',/column)
;
options=['Toggle Terms','Show Current Fit','Plot Residuals','Delete Stars',$
               'CCD offset correction','Star Info','Quit']
main = cw_bgroup(base,options,/row,uvalue='main')
;
;
;
plotme = widget_draw(base,xsize=720,ysize=256,/button_events,$
                     uvalue='plotme')
messages = widget_text(base,ysize=25,/scroll,/wrap,value='Initial fit:',$
                       uvalue='messages')
widget_control,base,/realize
widget_control,plotme,get_value=index
wset,index
;
widget_control,plotme,draw_button_events=0
widget_control,main,/hourglass
;
;
PhotAlg

end
