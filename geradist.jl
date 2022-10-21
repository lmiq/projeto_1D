using Distributions

function geradist(npart, Lref, stdev, disc)
    l = Int.(round.(rand(myLogNormal(Lref/disc, stdev/disc), npart))) .+ 1; dr = disc / Lref;
    for ii = length(l):-1:2
        l[ii] = sum(l[1:ii])
    end

    return l, dr
end

function myLogNormal(m,std)
    γ = 1+std^2/m^2
    μ = log(m/sqrt(γ))
    σ = sqrt(log(γ))

    return LogNormal(μ,σ)
end