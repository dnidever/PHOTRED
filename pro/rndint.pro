function rndint

; returns a random integer
;
; WARNING!! Currently if RNDINT is called by multiple
; IDL sessons at exactly the same time then it can
; output the *SAME* "random" digits.  Not sure how
; to get around that.  The "seed" is initialized by
; the current time.

;time = systime(1)
;seed = time*1000.-floor(time*1000.)
;seed = seed*100
rnd = randomu(seed,/uniform)

integer = floor(rnd*10.)
if integer eq 10 then integer = 9

;stop

return,integer

end
