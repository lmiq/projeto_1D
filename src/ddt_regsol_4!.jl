function ddt_regsol_4!(ddt, c, l::Vector{Int64}, omga::Float64, kapp::Float64, i, dr::Float64, A::Vector{Float64})
    init = 1;
    for jj = eachindex(l)
        fim = l[jj];
        ddt_part!(ddt, c, init, fim, jj, l, omga, kapp, i, dr, A);
        init = l[jj] + 1;
    end

    return nothing
end

function ddt_part!(ddt, c, init::Int64, fim::Int64, jj::Int64, l::Vector{Int64}, omga::Float64, kapp::Float64, i, dr::Float64, A::Vector{Float64})
    s = log(c[init] / (1 - c[init]))
    H = omga * (1 - 2 * c[init])
    pg = kapp*(6*(c[init+1] - c[init])/dr^2)
    u = s + H - pg
    f1 = 0.
    @inbounds for ii = 1:fim-init-1
        i2 = init + ii
        i3 = i2 + 1
        im3 = i2 - 1
        s = log(c[i2] / (1 - c[i2]))
        H = omga*(1 - 2*c[i2])
        pg = kapp*((c[i3] - c[im3])/(ii*dr*dr) + (c[i3] + c[im3] - 2*c[i2])/dr^2)
        du = (s + H - pg - u) / dr
        cm = 2/(inv(c[im3]) + inv(c[i2]))
        f2 = A[im3]*(-cm * (1-cm) * du)
        ddt[im3] = (f1 - f2)
        f1 = f2
        u = s + H - pg
    end
    s = log(c[fim] / (1 - c[fim]))
    H = omga * (1 - 2 * c[fim])
    pg = kapp * (2*c[fim-1] - 2*c[fim]) / dr^2
    du = (s + H - pg - u) / dr
    cm = 2/(inv(c[fim-1]) + inv(c[fim]))
    f2 = A[fim-1]*(-cm * (1-cm) * du)
    ddt[fim-1] = (f1 - f2)
    ddt[fim] = (f2 - A[fim]*i[jj])
    return nothing
end