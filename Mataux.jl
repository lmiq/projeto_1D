using LinearAlgebra

function Mataux(l, dr)
    Mvec = [3/4 for ii = 1:l[end]];
    Muvec, Mlvec = [[1/8 for ii = 1:l[end]-1] for jj = 1:2]; Muvec[end] = 1/4; Mlvec[1] = 1/4;
    Vvec = Array{Float64}(undef, l[end]);
    A = Array{Float64}(undef, l[end]);
    Asup = Array{Float64}(undef, length(l));

    init = 0;
    for jj = eachindex(l)
        if jj < length(l)
            Muvec[l[jj]] = 0; Muvec[l[jj]-1] = 1/4; 
            Mlvec[l[jj]] = 0; Mlvec[l[jj]+1] = 1/4;
        end

        nv = l[jj]- init;
        Vvec[init+1] = 4*pi*(dr^3/24); A[init+1] = 4*pi*(dr/2)^2;
        r = dr;
        for ii = 2:nv-1
            Vvec[init+ii] = 4*pi*(r^2*dr + dr^3/12); A[init+ii] = 4*pi*(r + dr/2)^2;
            r = ii*dr;
        end
        Vvec[init+nv] = 4*pi*(r^3/3 - (r - dr/2)^3/3); Asup[jj] = A[init+nv] = 4*pi*r^2;
        
        init = l[jj]
    end

    Md = [3/4 for ii = 1:l[end]]; Ml = Mu = [1/8 for ii = 1:l[end]-1]; Ml[1] = 1/4; Mu[end] = 1/4;
    for ii = 2:length(l)-1
        Mu[l[ii] - 1] = 1/4; Mu[l[ii]] = 0;
        Ml[l[ii-1]+1] = 1/4; Ml[l[ii-1]] = 0;
    end

    M = Tridiagonal(Mlvec, Mvec, Muvec); V = Diagonal(Vvec);
    MV = M*V; Vtot = sum(Vvec);

    return MV, Vtot, A, Asup;
end