pro ia_trim,xshift,yshift,xsize,ysize,trimsection,vignette

; Found this in the IMCENTROID source code:
; /iraf/iraf/pkg/images/immatch/src/listmatch/t_imctroid.x

;# IA_TRIM -- Compute the trim section.
;
;procedure ia_trim (cp)
;
;pointer cp                      #I center structure pointer
;
;real    xlo, xhi, ylo, yhi, xmin, ymin
;int     ixlo, ixhi, iylo, iyhi, ixlonew, ixhinew, iylonew, iyhinew, i
;int     vxlo, vxhi, vylo, vyhi          # vignetted versions
;bool    firsttime
;
;begin

nimages = n_elements(xshift)

        firsttime = 1
        for i=0,nimages-1 do begin

        ;firsttime = true
        ;do i = 1, NIMAGES(cp) {

            ;if (IS_INDEFR(XSHIFT(cp,i)) || IS_INDEFR(YSHIFT(cp,i)))
            ;    next

            ;# Compute limits.
            ;xlo = 1. + XSHIFT(cp,i)
            ;ylo = 1. + YSHIFT(cp,i)
            ;xhi = XSIZE(cp,i) + XSHIFT(cp,i)
            ;yhi = YSIZE(cp,i) + YSHIFT(cp,i)
            xlo = 1.0 + xshift[i]
            ylo = 1.0 + yshift[i]
            xhi = xsize[i] + xshift[i]
            yhi = ysize[i] + yshift[i]

            ;ixlonew = int (xlo)
            ;if (xlo > ixlonew)                  # round up
            ;    ixlonew = ixlonew + 1
            ixlonew = ceil(xlo)

            ;ixhinew = int (xhi)
            ;if (xhi < ixhinew)                  # round down
            ;    ixhinew = ixhinew - 1
            ixhinew = floor(xhi)

            ;iylonew = int (ylo)                 # round up
            ;if (ylo > iylonew)
            ;    iylonew = iylonew + 1
            iylonew = ceil(ylo)

            ;iyhinew = int (yhi)                 # round down
            ;if (yhi < iyhinew)
            ;    iyhinew = iyhinew - 1
            iyhinew = floor(yhi)

            if (firsttime) then begin
                ixlo = ixlonew
                ixhi = ixhinew
                iylo = iylonew
                iyhi = iyhinew

                ;xmin = XSIZE(cp,i)
                ;ymin = YSIZE(cp,i)
                xmin = xsize[i]
                ymin = ysize[i]

                firsttime = 0
             endif else begin
                ;ixlo = max (ixlo, ixlonew)
                ;ixhi = min (ixhi, ixhinew)
                ;iylo = max (iylo, iylonew)
                ;iyhi = min (iyhi, iyhinew)
                ixlo = max([ixlo,ixlonew])
                ixhi = min([ixhi,ixhinew])
                iylo = max([iylo,iylonew])
                iyhi = min([iyhi,iyhinew])

                ;xmin = min (XSIZE(cp,i), xmin)
                ;ymin = min (YSIZE(cp,i), ymin)
                xmin = min([xsize[i], xmin])
                ymin = min([ysize[i], ymin])
            endelse
        endfor

        ;# Don't bother to complain.
        ;if (firsttime)
        ;    return

        ;call printf ("\n")

        ;# Vignetting is possible downstream since imshift and other tasks
        ;# preserve the size of the input image.

        ;vxlo = max (1, min (ixlo, int(xmin)))
        ;vxhi = max (1, min (ixhi, int(xmin)))
        ;vylo = max (1, min (iylo, int(ymin)))
        ;vyhi = max (1, min (iyhi, int(ymin)))
        vxlo = max([1, min([ixlo, fix(xmin)]) ])
        vxhi = max([1, min([ixhi, fix(xmin)]) ])
        vylo = max([1, min([iylo, fix(ymin)]) ])
        vyhi = max([1, min([iyhi, fix(ymin)]) ])

        ;if (vxlo != ixlo || vxhi != ixhi || vylo != iylo || vyhi != iyhi) {
        ;    call eprintf ("#Vignette_Section = [%d:%d,%d:%d]\n")
        ;        call pargi (vxlo)
        ;        call pargi (vxhi)
        ;        call pargi (vylo)
        ;        call pargi (vyhi)
        ;}
        if (vxlo ne ixlo or vxhi ne ixhi or vylo ne iylo or vyhi ne iyhi) then begin
          vignette = [vxlo,vxhi,vylo,vyhi]
          print,'Vignette_Section = [',strtrim(vxlo,2),':',strtrim(vxhi,2),',',strtrim(vylo,2),':',strtrim(vyhi,2),']'
        endif

        ;# Output the trim section.
        ;call printf ("#Trim_Section = [%d:%d,%d:%d]\n")
        ;    call pargi (ixlo)
        ;    call pargi (ixhi)
        ;    call pargi (iylo)
        ;    call pargi (iyhi)
        trimsection = [ixlo,ixhi,iylo,iyhi]
        print,'Trim_Section = [',strtrim(ixlo,2),':',strtrim(ixhi,2),',',strtrim(iylo,2),':',strtrim(iyhi,2),']'

       ; call flush (STDOUT)
end
