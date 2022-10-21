import Pkg; Pkg.activate("projeto_1D");
include("config.jl"); include("salvasimul.jl");

function Run_polipart(npart, Lref, stdev, disc)
    prob, p, e, Kb, T, tref, trun = config(npart, Lref, stdev, disc);

    @time sol = solve(prob, Rodas4(), saveat = trun/200, abstol = 1e-6, reltol = 1e-6, progress = true);

    #Salvando a solução
    salvasimul(sol, "", p, e, Kb, T, tref)

    return sol
end