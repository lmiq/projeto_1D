
function salvasimul(sol, sufix, p, e, Kb, T, tref)
    l = p[1]::Vector{Int64}; 
    k0 = p[2]::Float64; alpha = p[3]::Float64; omga = p[4]::Float64; uref = p[5]::Float64 ; I = p[6]::Float64; dr = p[7]::Float64; kapp = p[8]::Float64;
    i, i0, uend = [Array{Float64}(undef, length(l)) for ii = 1:3]; 
    Asup = p[12]::Vector{Float64};
    A = p[13]::Vector{Float64};

    arquivou = open("u"*sufix*".txt", "w");
    arquivot = open("t"*sufix*".txt", "w");
    arquivol = open("l"*sufix*".txt", "w");
    arquivov = open("v"*sufix*".txt", "w");

    for ii in eachindex(sol.t)
        phi = Solve_BV2(i, i0, uend, Asup, sol.u[ii], l, k0, alpha, omga, uref, I);
        v = string(-1*phi*Kb*T/e);
        write(arquivov, v * "\n");

        auxt = string(sol.t[ii]*tref/3600);
        write(arquivot, auxt * "\n");

        for jj in eachindex(sol.u[ii])
            auxu = string(sol.u[ii][jj];)
            write(arquivou, auxu * " ")
        end
        write(arquivou, "\n")
    end

    for ii in eachindex(l)
        auxl = string(l[ii]);
        write(arquivol, auxl * " ");
    end

    close(arquivou);
    close(arquivot);
    close(arquivol);
    close(arquivov);
end