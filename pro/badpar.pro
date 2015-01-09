;+
; NAME:
;     badpar
; PURPOSE: (one line)
;     Validate an input parameter against valid entries.
; DESCRIPTION:
;
;     This is a general parameter checking function for validating input
;     quantities in other procedures and functions.  This routine will
;     generate an error message indicating what is wrong with the item.
;
;     Example of use:
;
;     pro foo,array
;     if badpar(array,[4,5],2,CALLER='foo') then return
;        .
;        . code for foo .
;        .
;     end
;
;
;     This would cause an immediate return to the routine that called foo
;     with an error message if the input was not either floating or double
;     and 2 dimensional.
;
;     As of IDL v3.0, these are the recognized type codes (see 1-218 in
;        reference guide).
;
;        Type
;        Code     Data Type
;        ----    -----------------------------
;          0      Undefined
;          1      Byte
;          2      Integer
;          3      Longword integer
;          4      Floating point
;          5      Double-precision floating
;          6      Complex floating
;          7      String
;          8      Structure
;          9      Double-precision complex
;         10      Pointer
;         11      Object reference
;         12      Unsigned integer
;         13      Unsigned Longword integer
;         14      64-bit integer
;         15      Unsigned 64-bit integer
;
; CATEGORY:
;  Utility
; CALLING SEQUENCE:
;     val = badpar(param,goodtype,goodrank)
; INPUTS:
;     param    - IDL variable to validate.
;     goodtype - Scalar or vector of type codes that are valid.
;     goodrank - Scalar or vector of valid ranks.
; OPTIONAL INPUT PARAMETERS:
; KEYWORD PARAMETERS:
;     CALLER   - String identifying the calling routine.
;     DEFAULT  - Value to return in param if undefined and undefined allowed.
;     DIMEN    - Dimensions of variable.
;     NPTS     - Total number of elements in variable.
;     RANK     - Rank of variable.
;     TYPE     - Type of variable.
; OUTPUTS:
;     Return value is true if the parameter is bad.  False if good.
; COMMON BLOCKS:
; SIDE EFFECTS:
; RESTRICTIONS:
; PROCEDURE:
; MODIFICATION HISTORY:
;     3/24/93 - Written by Marc W. Buie, Lowell Observatory.
;     4/27/93 - MWB, added TYPE and DEFAULT keywords.
;     2000/11/1, MWB, added new types for IDL v5.4
;-
function badpar,param,goodtype,goodrank, $
            CALLER=caller, DEFAULT=default, DIMEN=dimen, $
            NPTS=npts, RANK=rank, TYPE=type

   errval = 0

   sz = size(param)

   rank = sz[0]
   type = sz[rank+1]
   npts = sz[rank+2]

   err1=''
   err2=''

   if rank eq 0 then dimen=0 else dimen = sz[1:rank]

   z=where(goodtype eq type, count)
   if count eq 0 then begin
      case type of
         0 :    err1 = 'Undefined variable is not allowed.'
         1 :    err1 = 'Byte variable type is not allowed.'
         2 :    err1 = 'Integer variable type is not allowed.'
         3 :    err1 = 'Longword integer variable type is not allowed.'
         4 :    err1 = 'Floating point variable type is not allowed.'
         5 :    err1 = 'Double-precision floating point variable type is not allowed.'
         6 :    err1 = 'Complex floating point variable type is not allowed.'
         7 :    err1 = 'String variable type is not allowed.'
         8 :    err1 = 'Structure variable type is not allowed.'
         9 :    err1 = 'Double-precision complex type is not allowed.'
         10:    err1 = 'Pointer type is not allowed.'
         11:    err1 = 'Object reference type is not allowed.'
         12:    err1 = 'Unsigned integer type is not allowed.'
         13:    err1 = 'Unsigned longword integer type is not allowed.'
         14:    err1 = '64-bit integer type is not allowed.'
         15:    err1 = 'Unsigned 64-bit integer type is not allowed.'
         else : err1 = 'Unrecognized variable type code.  Impossible!'
      endcase
      errval=1
   endif

   if type ne 0 then begin
      z=where(goodrank eq rank, count)
      if count eq 0 then begin
         case rank of
            0 :    err2 = 'Scalar variables are not allowed.'
            1 :    err2 = 'Vector variables are not allowed.'
            2 :    err2 = '2-D variables are not allowed.'
            3 :    err2 = '3-D variables are not allowed.'
            4 :    err2 = '4-D variables are not allowed.'
            5 :    err2 = '5-D variables are not allowed.'
            6 :    err2 = '6-D variables are not allowed.'
            7 :    err2 = '7-D variables are not allowed.'
            8 :    err2 = '8-D variables are not allowed.'
            else : err2 = 'Unrecognized variable rank.  Impossible!'
         endcase
         errval=1
      endif
   endif

   if errval then begin
      if not keyword_set(caller) then caller = ''
      print,caller,'Illegal variable encountered.'
      if err1 ne '' then print,err1
      if err2 ne '' then print,err2
      return,errval
   endif

   if type eq 0 then begin
      szd = size(default)
      if szd[szd[0]+1] ne 0 then begin
         param = default
         sz    = size(param)
         rank  = sz[0]
         type  = sz[rank+1]
         npts  = sz[rank+2]
      endif
   endif

   return,errval

end
