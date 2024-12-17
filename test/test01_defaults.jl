import EcotoxSystems as EctSys
import EcotoxSystems: defaultparams

@testset "Default parameters" begin  

    p = defaultparams
    p.glb.t_max = 56.
    p.spc.Z = Dirac(1.)
    global sim = ODE_simulator(p, reltol = 1e-3)


    @test isapprox(maximum(sim.H), p.spc.H_p, rtol = 0.1) 
    @test isapprox(maximum(sim.S), EctSys.calc_S_max(p.spc), rtol = 0.1)
end;

@time EcotoxSystems.ODE_simulator(EcotoxSystems.defaultparams);

using BenchmarkTools
@benchmark EctSys.ODE_simulator(EctSys.defaultparams)

VSCodeServer.@profview_allocs [EctSys.ODE_simulator(EctSys.defaultparams) for _ in 1:100]


#=
Basic test of @replicates macro
=#

import EcotoxSystems: @replicates
using Chain

@testset "@replicates macro" begin
    p = EctSys.params()
    p.spc.Z = Truncated(Normal(1., 0.1), 0, Inf)
    sim = @replicates EctSys.ODE_simulator(p) 10

    plt = @df sim plot(
        plot(:t, :S, group = :replicate, color = 1),
        plot(:t, :H, group = :replicate, color = 1)
    )

    display(plt)

    cvs = @chain sim begin # compute coefficients of variation in final values
        groupby(:replicate)
        combine(_) do df
            return DataFrame(
                S_max = maximum(df.S),
                H_max = maximum(df.H),
                R_max = maximum(df.R)
            )
        end
        (
            S = std(_.S_max) / mean(_.S_max),
            H = std(_.H_max) / mean(_.H_max),
            R = std(_.R_max) / mean(_.R_max)
        )
    end

    @test cvs.S > 0.025 # test for plausible coefficients of variation in the final values
    @test cvs.H > 0.025
    @test cvs.R > 0.025
end;
