; $Id: running_gdl.pro,v 1.1 2010/05/06 10:18:17 lmugnier Exp $
; Copyright (C) 2010 Laurent Mugnier.

; This program is free software: you can redistribute it and/or modify
; it under the terms of the GNU General Public License as published by
; the Free Software Foundation, either version 3 of the License, or
; (at your option) any later version.

; This program is distributed in the hope that it will be useful,
; but WITHOUT ANY WARRANTY; without even the implied warranty of
; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
; GNU General Public License for more details.

; You should have received a copy of the GNU General Public License
; along with this program.  If not, see <http://www.gnu.org/licenses/>.

FUNCTION RUNNING_GDL, VERBOSE = verbose, VERSION = version, HELP = help

;+
; NAME:
;	RUNNING_GDL - check whether we are running GDL or IDL
;
; CATEGORY:
;	Help Routines
;
; CALLING SEQUENCE:
;	result = RUNNING_GDL([/VERBOSE] [, /VERSION] [, /HELP])
;
; PURPOSE:
;   This function checks whether we are running GDL or IDL.
;   it returns 1 if running GDL, 0 if running IDL.
;
;   It may be useful for instance to change your !path depending on the
;   running engine or to make routines handle gracefully things not yet
;   implemented in GDL.
;
; POSITIONAL PARAMETERS:
;   none.
;
; KEYWORD PARAMETERS:
;	VERBOSE  : (optional input) print a message telling whether we are
;	           running GDL or IDL.
;
;   /VERSION : (optional input) prints version number before execution.
;
;   /HELP    : (optional input) prints the documentation and exits.
;
; AUTHORS:
;   $Author: lmugnier $
;
; RESTRICTIONS:
;   This program is copyright (c) Laurent Mugnier, 2010.
;   This program is free software: you can redistribute it and/or modify
;   it under the terms of the GNU General Public License as published by
;   the Free Software Foundation, either version 3 of the License, or
;   (at your option) any later version.
;
; EXAMPLE:
;   You may want to add some librairies to your path only if you are running
;   GDL. To do this, use:
;	if running_gdl() then $
;	    !path = !path + ':' + expand_path('+~/gdl/lib')
;
; SEE ALSO:
;   DEFSYSV (called by this routine)
;
; HISTORY:
;   $Log: running_gdl.pro,v $
;   Revision 1.1  2010/05/06 10:18:17  lmugnier
;   Removed call to 'routine_courante()'
;
;   Revision 1.0  2010/05/06 10:14:18  lmugnier
;   Initial revision
;
;-

on_error,2
IF keyword_set(version) THEN $
    printf, -2, '% RUNNING_GDL: $Revision: 1.1 $, $Date: 2010/05/06 10:18:17 $'

IF (n_params() NE 0) OR keyword_set(help) THEN BEGIN
    message, 'Help required or incorrect syntax. Documentation:', /INFO
    doc_library, 'running_gdl'
    retall
ENDIF

DEFSYSV, '!GDL', EXISTS = WeAreRunningGDL

IF keyword_set(verbose) THEN $
    IF (WeAreRunningGDL EQ 1) THEN $
        print, 'We are running GDL, compatible with IDL version ', $
               !version.release $
    ELSE $
        print, 'We are running IDL version ', !version.release

return, WeAreRunningGDL

END
