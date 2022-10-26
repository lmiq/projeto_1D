module Projeto1D

    using PreallocationTools
    using DifferentialEquations
    using Distributions
    using SparseArrays
    using LinearAlgebra

    include("geradist.jl")
    include("ddt_polipart!.jl")
    include("Mataux.jl")
    include("Configsimul.jl")
    include("jacpattern.jl")
    include("./config.jl")
    include("./salvasimul.jl")
    include("Solve_BV!.jl");
    include("Solve_BV2.jl");
    include("ddt_regsol_4!.jl");
    include("./Run_polipart.jl")

    function main()
        Run_polipart(5, 80e-9, 8e-9, 0.5e-9)
    end

end # module