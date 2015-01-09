function wcstnxcor, xi, eta, tnx

    ;compute normalized values for xi and eta used in the 
    ;chebyshev and legendre polynomial functions:

    ; see http://iraf.noao.edu/projects/ccdmosaic/tnx.html for details
    
    xin = (2.d * xi - (tnx.ximax + tnx.ximin)) / (tnx.ximax - tnx.ximin)
    etan = (2.d * eta - (tnx.etamax + tnx.etamin)) / (tnx.etamax - tnx.etamin)

    ncoor = n_elements(xi)

    ; Calculating the polynomial terms Pmn
    CASE (tnx.fun_type) of

    ; Chebyshev
    1: begin

         ; Pmn = Pm(xin) * Pn(etan)
         ;
	 ; P0(xin) = 1.0
	 ; P1(xin) = xin
	 ; Pm+1(xin) = 2.0 * xin * Pm(xin) - Pm-1(xin) 
         ;
	 ; P0(etan) = 1.0
	 ; P1(etan) = etan
	 ; Pn+1(etan) = 2.0 * etan * Pn(etan) - Pn-1(etan)

         pm = dblarr(4,ncoor)
         pm[0,*] = 1.
         pm[1,*] = xin
         for m=1,2 do pm[m+1,*] = 2.0d * xin * Pm[m,*] - Pm[m-1,*] 
         pn = dblarr(4,ncoor)
         pn[0,*] = 1.
         pn[1,*] = etan
         for n=1,2 do pn[n+1,*] = 2.0d * etan * Pn[n,*] - Pn[n-1,*] 

         sum = dblarr(ncoor)
         for i = 0, 9 do $
           sum = sum + tnx.c[i] * Pm[tnx.m[i],*] * Pn[tnx.n[i],*]

       end ; chebyshev

    ; Legendre
    2: begin

         ; Pmn = Pm(xin) * Pn(etan)     (legendre)
         ; 
         ; P0(xin) = 1.0
         ; P1(xin) = xin
         ; Pm+1(xin) = ((2m+1) * xin * Pm(xin) - m * Pm-1(xin))/ (m+1)   
         ; 
         ; P0(etan) = 1.0
         ; P1(etan) = etan
         ; Pn+1(etan) = ((2n+1) * etan * Pn(etan) - n * Pn-1(etan))/ (n+1)

         pm = dblarr(4,ncoor)
         pm[0,*] = 1.
         pm[1,*] = xin
         for m=1,2 do pm[m+1,*] = ((2.0d*m+1) * xin * Pm[m,*] - m * Pm[m-1,*])/ (m+1) 
         pn = dblarr(4,ncoor)
         pn[0,*] = 1.
         pn[1,*] = etan
         for n=1,2 do pn[n+1,*] = ((2.0d*m+1) * etan * Pn[n,*] - n * Pn[n-1,*])/ (n+1) 

         sum = dblarr(ncoor)
         for i = 0, 9 do $
           sum = sum + tnx.c[i] * Pm[tnx.m[i],*] * Pn[tnx.n[i],*]

       end ; legendre


    ; Polynomial
    3: begin

         ; Pmn = xi ** m * eta ** n    (polynomial)

         sum = dblarr(ncoor)
         for i = 0, 9 do $
           sum = sum + tnx.c[i] * (xi^tnx.m[i]) * (eta^tnx.n[i])

       end ; polynomial

    ENDCASE
   
    return, sum

end
