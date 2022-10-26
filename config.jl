include("geradist.jl"); include("ddt_polipart!.jl"); include("Mataux.jl"); include("Configsimul.jl"); include("teste.jl"); include("jacpattern.jl");

function config(npart, Lref, stdev, disc)
    #Gera distribuição de partículas
    l, dr = geradist(npart, Lref, stdev, disc)

    #matrizes auxiliares
    MV, Vtot, A, Asup = Mataux(l, dr);

    #Configuração
    c0, trun, p, e, Kb, T, tref = Configsimul(l, Vtot, A, Asup, dr, Lref)
    #teste(1, ddt_polipart!, c0, p);
    jc = jacpattern(l);

    #resolvendo o sistema 1
    func = ODEFunction(ddt_polipart!, mass_matrix = MV, jac_prototype = jc);
    prob = ODEProblem(func, c0, (0, trun), p);

    return prob, p, e, Kb, T, tref, trun
end