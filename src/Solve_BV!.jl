function Solve_BV!(i, i0, uend, Asup::Vector{Float64}, c, l::Vector{Int64}, k0::Float64, alpha::Float64, omga::Float64, uref::Float64, I::Float64)
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

    return nothing
end

function currtroc(cloc, k0::Float64, alpha::Float64, omga::Float64)
    return k0*(ativ(cloc, omga))^alpha/(1/(1-cloc))
end

function ativ(cloc, omga::Float64)
    return (cloc/(1-cloc))*exp(omga*(1-2*cloc))
end

function poteq(cloc, omga::Float64, uref::Float64)
    return uref + log(cloc/(1-cloc)) + omga*(1 - 2*cloc)
end

function IdI(i0, uend, Asup::Vector{Float64}, alpha::Float64, phi, I::Float64)
    f = -I;
    df = 0.;
    for ii = eachindex(i0)
        eta = phi - uend[ii];
        f = f + Asup[ii]*i0[ii]*(exp(-alpha*eta) - exp((1-alpha)*eta));
        df = df + Asup[ii]*i0[ii]*exp(-alpha*eta)*(alpha*(exp(eta)-1)-exp(eta));
    end

    return phi - f/df;
end