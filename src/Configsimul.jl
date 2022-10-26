function Configsimul(l, Vtot, A, Asup, dr, Lref)
    #Parâmetros e adimensionalização
    D = 4e-16; e = 1.6e-19; rhos = 1.35e28; Kb = 1.381e-23; T = 298; alpha = 0.5;
    tref = Lref^2/D; k0 = 3.6*Lref/(rhos*e*D); omga = 1.43e-20/(Kb*T); kapp = 1e-10/(rhos*Kb*T*Lref^2); uref = -1*1.55*e/(Kb*T);

    #Fluxo imposto e tempo de (des)carga
    cr = 0.1*(tref/3600); trun = 0.98/(cr); I = -1*Vtot*cr; 

    #Variáveis auxiliares
    i, i0, uend = [Array{Float64}(undef, length(l)) for ii = 1:3];

    #Condição inicial e parâmetros
    c0 = [0.01 for ii = 1:l[end]];
    p = (l, k0, alpha, omga, uref, I, dr, kapp, dualcache(i, 12), dualcache(i0, 12), dualcache(uend, 12), Asup, A);

    return c0, trun, p, e, Kb, T, tref
end