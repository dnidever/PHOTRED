function parsetnx, sin

    s = sin

    ; all we are interested here is the bit between quotes
    q1 =  strpos(s, '"')
    q2 =  strpos(s, '"', q1+1)
    scor = strsplit(strmid(s,q1+1,q2-q1-1)," ", /extract)
    
    fun_type = fix(scor[0])
    xiorder = fix(scor[1])
    etaorder = fix(scor[2])
    cross_type = fix(scor[3])
    ximin = double(scor[4])
    ximax = double(scor[5])
    etamin = double(scor[6])
    etamax = double(scor[7])
    
    m = intarr(10)
    n = intarr(10)
    c = dblarr(10)
    i = 0
    for nn = 0, 3 do begin
        for mm = 0, 3-nn do begin
            m[i] = mm
            n[i] = nn
            c[i] = double(scor[8+i])
            i = i + 1
        end
    end

    tnx = {fun_type:fun_type, xiorder:xiorder, etaorder:etaorder,$
           cross_type:cross_type, ximin:ximin, ximax:ximax, $
           etamin:etamin, etamax:etamax, m:m, n:n, c:c}

    return, tnx

end

