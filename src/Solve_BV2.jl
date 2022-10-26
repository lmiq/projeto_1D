function Solve_BV2(i, i0, uend, Asup::Vector{Float64}, c, l::Vector{Int64}, k0::Float64, alpha::Float64, omga::Float64, uref::Float64, I::Float64)
    for ii = eachindex(l)
        i0[ii] = currtroc(c[l[ii]], k0, alpha, omga);
        uend[ii] = poteq(c[l[ii]], omga, uref);
    end

    phiprev = uref; 
    phipos = IdI(i0, uend, Asup, alpha, phiprev, I);

    while abs(phipos - phiprev) > 1e-6 + 1e-6*abs(phiprev)
        phiprev = phipos;
        phipos = IdI(i0, uend, Asup, alpha, phiprev, I);
    end

    for ii = eachindex(l)
        eta = phipos - uend[ii];
        i[ii] = i0[ii]*(exp.(-alpha*(eta)) - exp.((1-alpha)*(eta)));
    end

    return phipos
end