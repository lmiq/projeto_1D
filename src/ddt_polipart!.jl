
function ddt_polipart!(ddt, c, p, t)
    l = p[1]::Vector{Int64}; 
    k0 = p[2]::Float64; alpha = p[3]::Float64; omga = p[4]::Float64; uref = p[5]::Float64; I = p[6]::Float64; dr = p[7]::Float64; kapp = p[8]::Float64;
    i = get_tmp(p[9], c); i0 = get_tmp(p[10], c); uend = get_tmp(p[11], c);
    Asup = p[12]::Vector{Float64};
    A = p[13]::Vector{Float64};
    
    Solve_BV!(i, i0, uend, Asup, c, l, k0, alpha, omga, uref, I);
    ddt_regsol_4!(ddt, c, l, omga, kapp, i, dr, A);

    return nothing
end