;-------------------------------------------------------------
;+
; NAME:
;       STRESS
; PURPOSE:
;       String edit by sub-string. Precede, Follow, Delete, Replace.
; CATEGORY:
; CALLING SEQUENCE:
;       new = stress(old,cmd,n,oldss,newss,ned)
; INPUTS:
;       old = string to edit.                               in
;       cmd = edit command:                                 in
;         'P' = precede.
;         'F' = follow.
;         'D' = delete.
;         'R' = replace.
;       n = occurrence number to process (0 = all).         in
;       oldss = reference substring.                        in
;       oldss may have any of the following forms:
;         1. s	  a single substring.
;         2. s...    start at substring s, end at end of string.
;         3. ...e    from start of string to substring e.
;         4. s...e   from subs s to subs e.
;         5. ...     entire string.
;       newss = substring to add. Not needed for 'D'        in
; KEYWORD PARAMETERS:
; OUTPUTS:
;       ned = number of occurrences actually changed.       out
;       new = resulting string after editing.               out
; COMMON BLOCKS:
; NOTES:
; MODIFICATION HISTORY:
;       Written by R. Sterner, 6 Jan, 1985.
;       RES --- 23 May, 1988 fixed a bug in SSTYP = 2.
;       Converted to SUN 13 Aug, 1989 --- R. Sterner. (FOR loop change).
;       Johns Hopkins University Applied Physics Laboratory.
;
; Copyright (C) 1985, Johns Hopkins University/Applied Physics Laboratory
; This software may be used, copied, or redistributed as long as it is not
; sold and this copyright notice is reproduced on each copy made.  This
; routine is provided as is without any express or implied warranties
; whatsoever.  Other limitations apply as described in the file disclaimer.txt.
;-
;-------------------------------------------------------------
 
	FUNCTION STRESS,STRNG,CMDX,N,OLD,NEW,NED, help = h
 
	if (n_params(0) lt 3) or keyword_set(h) then begin
	  print,' String edit by sub-string. Precede, Follow, Delete, Replace.'
	  print,' new = stress(old,cmd,n,oldss,newss,ned)
	  print,'   old = string to edit.                               in'
	  print,'   cmd = edit command:                                 in'
	  print,"     'P' = precede.
	  print,"     'F' = follow.
	  print,"     'D' = delete.
	  print,"     'R' = replace.
	  print,'   n = occurrence number to process (0 = all).         in'
	  print,'   oldss = reference substring.                        in'
	  print,'   oldss may have any of the following forms:
	  print,'     1. s	  a single substring.
	  print,'     2. s...    start at substring s, end at end of string.
	  print,'     3. ...e    from start of string to substring e.
	  print,'     4. s...e   from subs s to subs e.
	  print,'     5. ...     entire string.
	  print,"   newss = substring to add. Not needed for 'D'        in"
	  print,'   ned = number of occurrences actually changed.       out'
	  print,'   new = resulting string after editing.               out'
	  return, -1
	endif
 
	CMD = STRUPCASE(CMDX)
	PDOT = STRPOS(OLD,'...')
	SSL = STRLEN(OLD)
	SSTYP = 0
	POS1 = -1
	POS2 = -1
	RSTR = STRNG
	IF (PDOT EQ -1) THEN SSTYP = 1
;	IF ((PDOT>0) EQ SSL-3) THEN SSTYP = 2
	IF (PDOT GT 0) AND (PDOT EQ SSL-3) THEN SSTYP = 2
	IF (PDOT EQ 0) AND (SSL GT 3) THEN SSTYP = 3
	IF (PDOT GT 0) AND (PDOT LT SSL-3) THEN SSTYP = 4
	IF (PDOT EQ 0) AND (SSL EQ 3) THEN SSTYP = 5
	NED = 0		; number of occurrences actually changed.
 
 
	CASE SSTYP OF
1:	  BEGIN
	    S = OLD
	    E = ''
	  END
2:	  BEGIN
	    S = STRSUB(OLD,0,SSL-4)
	    E = ''
    	  END
3:  	  BEGIN
	    S = ''
	    E = STRSUB(OLD,3,SSL-1)
	  END
4:  	  BEGIN
	    S = STRSUB(OLD,0,PDOT-1)
	    E = STRSUB(OLD,PDOT+3,SSL-1)
	  END
5:  	  BEGIN
	    S = ''
	    E = ''
	  END
ELSE: 	  PRINT, 'ERROR IN SSTYP'
	ENDCASE
 
 
;---------------  Find substring # N start  ---------------
	POS = -1
	nfor = n>1
LOOP:
	FOR I = 1, nfor DO BEGIN
	  POS = POS + 1
	  CASE SSTYP OF
    1:      POS = STRPOS(RSTR,S,POS)
    2:      POS = STRPOS(RSTR,S,POS)
    4:      POS = STRPOS(RSTR,S,POS)
    3:      POS = STRPOS(RSTR,E,POS)
    5:      POS = 0
	  ENDCASE
  	  IF POS LT 0 THEN GOTO, DONE
	ENDFOR
 
;----------  Find substring # N END  ----------------
    	CASE SSTYP OF
1:  	  BEGIN
	    POS1 = POS
	    POS2 = POS + STRLEN(S) - 1
	  END
2:  	  BEGIN
	    POS1 = POS
	    POS2 = STRLEN(RSTR) - 1
	  END
3:  	  BEGIN
	    POS1 = 0
	    POS2 = POS + STRLEN(E) - 1
	  END
4:  	  BEGIN
	    POS1 = POS
	    POS2 = STRPOS(RSTR,E,POS+1)
	    IF (POS2 LT 0) THEN GOTO, DONE
	    POS2 = POS2 + STRLEN(E) - 1
	  END
5:  	  BEGIN
	    POS1 = 0
	    POS2 = STRLEN(RSTR) - 1
	  END
	ENDCASE
 
;------------  edit string  --------------
    	CASE CMD OF
'P':  	  BEGIN
	    RSTR = STREP(RSTR,CMD,POS1,NEW)
	    POS = POS + STRLEN(NEW)
	  END
'F':  	  BEGIN
	    RSTR = STREP(RSTR,CMD,POS2,NEW)
	    POS = POS + STRLEN(NEW)
	  END
'R':  	  BEGIN
	    RSTR = STREP(RSTR,'D',POS1,POS2-POS1+1)
	    IF (POS1 GT 0) THEN $
	      RSTR = STREP(RSTR,'F',POS1-1,NEW)
	    IF (POS1 EQ 0) THEN $
	      RSTR = STREP(RSTR,'P',0,NEW)
	    POS = POS + STRLEN(NEW) - 1
	  END
'D':  	  BEGIN
	    RSTR = STREP(RSTR,CMD,POS1,POS2-POS1+1)
	    POS = POS - 1
	  END
ELSE: 	  BEGIN
	    PRINT, 'Error in cmd'
	    RETURN,RSTR
	  END
ENDCASE
 
	NED = NED + 1
	IF SSTYP EQ 5 THEN RETURN,RSTR
	IF N EQ 0 THEN GOTO, LOOP
 
DONE:
	RETURN, RSTR
	END
