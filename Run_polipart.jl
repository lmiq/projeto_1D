import Pkg; Pkg.activate("projeto_1D");

using DifferentialEquations
include("config.jl"); include("salvasimul.jl")

function Run_polipart(npart, Lref, stdev, disc)
    #Função principar para as simulações
    #npart: número de partículas simuladas
    #Lref: raio de partícula de médio
    #stdev: desvio padrão da distribuição de raios das partículas
    #disc: comprimento de discretização espacial interna das partíuclas 

    prob, p, e, Kb, T, tref, trun = config(npart, Lref, stdev, disc);
    

    #@time sol = solve(prob, Rodas4(), saveat = trun/200, abstol = 1e-6, reltol = 1e-6, progress = true);
    @time sol = solve(prob, FBDF(), saveat = trun/200, abstol = 1e-6, reltol = 1e-6, progress = true);
    

    #Salvando a solução
    salvasimul(sol, "", p, e, Kb, T, tref)


    return sol
end