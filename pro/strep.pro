;-------------------------------------------------------------
;+
; NAME:
;       STREP
; PURPOSE:
;       Edit a string by position. Precede, Follow, Replace, Delete.
; CATEGORY:
; CALLING SEQUENCE:
;       newstring = strep(string,cmd,p,ss,[iflg])
; INPUTS:
;       string = string to edit.                               in
;       cmd = edit command:                                    in
;         'P' = precede position p with substring ss.
;         'F' = follow position p with substring ss.
;         'R' = replace text starting at position p
;               with text from substring ss.
;         'D' = delete N characters starting at
;               position p.  The calling sequence for
;               this command is slightly different:
;               IFLG = STREP(string,'D',p,n,[iflg])
;               Where n = number of characters to delete.
;       p = character position to use.                         in
;           0 = first char.  Any number larger
;           than the string length = last char.
;       ss = substring to use.  For 'D' command                in
;            n is used instead of ss.
; KEYWORD PARAMETERS:
; OUTPUTS:
;       iflg = 0 for a successful edit,                        out
;       iflg = -1 for an error and no change to string.
;       newstring = edited string.                             out
; COMMON BLOCKS:
; NOTES:
; MODIFICATION HISTORY:
;       Written by R. Sterner, 27 Dec, 1984.
;       Converted to SUN 13 Aug, 1989 --- R. Sterner.
;       Johns Hopkins University Applied Physics Laboratory.
;
; Copyright (C) 1984, Johns Hopkins University/Applied Physics Laboratory
; This software may be used, copied, or redistributed as long as it is not
; sold and this copyright notice is reproduced on each copy made.  This
; routine is provided as is without any express or implied warranties
; whatsoever.  Other limitations apply as described in the file disclaimer.txt.
;-
;-------------------------------------------------------------
 
	FUNCTION strep,S,CMD,IP,SS,IFLG, help = h 
 
	IF (N_PARAMS(0) LT 4) or keyword_set(h) THEN BEGIN
	  PRINT,' Edit a string by position. Precede, Follow, Replace, Delete.'
	  print,' newstring = strep(string,cmd,p,ss,[iflg])
	  print,'   string = string to edit.                               in'
	  print,'   cmd = edit command:                                    in'
	  print,"     'P' = precede position p with substring ss.
	  print,"     'F' = follow position p with substring ss.
	  print,"     'R' = replace text starting at position p
	  print,'           with text from substring ss.
	  print,"     'D' = delete N characters starting at
	  print,'           position p.  The calling sequence for
	  print,'           this command is slightly different:
	  print,"           IFLG = STREP(string,'D',p,n,[iflg])
	  print,'           Where n = number of characters to delete.
	  print,'   p = character position to use.                         in'
	  print,'       0 = first char.  Any number larger
	  print,'       than the string length = last char.
	  print,"   ss = substring to use.  For 'D' command                in"
	  print,'        n is used instead of ss.
	  print,'   iflg = 0 for a successful edit,                        out'
	  print,'   iflg = -1 for an error and no change to string.
	  print,'   newstring = edited string.                             out'
	  RETURN, -1
	ENDIF
 
	N = STRLEN(S) - 1
	IF N LT 0 THEN $
	  IF CMD EQ 'D' THEN RETURN, '' ELSE RETURN, SS
	P = IP>0<N
	BS = BYTE(S)
	IF (CMD NE 'D') THEN  BSS = BYTE(SS)
	JFLG = 0
 
	CASE CMD OF
  'P':	BEGIN
	IF P EQ 0 THEN T = [BSS,BS] ELSE $
		T = [BS(0:P-1),BSS,BS(P:*)]
	END
  'F':	BEGIN
	IF P EQ N THEN BEGIN
	  T = [BS,BSS]
	ENDIF ELSE BEGIN
	  IF IP GE 0 THEN BEGIN
	    T = [BS(0:P),BSS,BS(P+1:*)]
	  ENDIF ELSE BEGIN
	    T = [BSS, BS]
	  ENDELSE
	ENDELSE
	END
  'R':  BEGIN
	NSS = STRLEN(SS)
	IF NSS GT N - P + 1 THEN $
		T = [BS(0:P-1),BSS] $
	ELSE BEGIN $
		T = BS
		T(P) = BSS
	ENDELSE
	END
  'D':	BEGIN
	ND = SS < (N + 1)
	IF P EQ 0 THEN $
	  IF ND EQ N+1 THEN T = '' ELSE $
	    T = [BS(ND:*)]
	IF P GT 0 THEN $
	  IF ND+P-1 GE N THEN $
	    T = [BS(0:P-1)] ELSE $
	    T = [BS(0:P-1),BS(ND+P:*)]
	END
  ELSE: JFLG = -1
  ENDCASE
 
	IF JFLG NE 0 THEN T = S
	IF JFLG EQ 0 THEN T = STRING(T)
	IF N_PARAMS(0) EQ 5 THEN IFLG = JFLG
 
 
	RETURN,T
 
	END
