begin # setting up a debkiss instance
    debkiss = SimplifiedEnergyBudget() |> instantiate
    debkiss.parameters.glb.t_max = 56.  
end

@testset "Default debkiss parameters" begin
    p = debkiss.parameters  
    p.spc.Z = Dirac(1.)  
    @time global sim = simulate(debkiss; reltol = 1e-3)
    @test isapprox(maximum(sim.H), p.spc.H_p, rtol = 0.1) 
    @test isapprox(maximum(sim.S), EcotoxSystems.calc_S_max(p.spc), rtol = 0.1)
end;

#@time sim = EcotoxSystems.ODE_simulator(p; alg = Rosenbrock23());
#using BenchmarkTools; 
#@benchmark EcotoxSystems.ODE_simulator(p; alg = Tsit5())
#VSCodeServer.@profview_allocs [EcotoxSystems.ODE_simulator(p) for _ in 1:100];

#=
Basic test of @replicates macro
=#

import EcotoxSystems: @replicates
using Chain

@testset "@replicates macro" begin
    debkiss.parameters.spc.Z = Truncated(Normal(1., 0.1), 0, Inf)
    sim = @replicates simulate(debkiss) 10

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

#=
Basic test of replicates function.
=#
@testset "replicates function" begin
    
    debkiss.parameters.spc.Z = Truncated(Normal(1., 0.1), 0, Inf)

    function simulate(p; kwargs...)
        return EcotoxSystems.ODE_simulator(
            debkiss.parameters;
            model = debkiss.complete_derivatives!, 
            statevars_init = debkiss.initialize_all_statevars,
            gen_ind_params = debkiss.generate_individual_params, 
            kwargs...
        )
    end

    sim = replicates(simulate, debkiss.parameters, 10)

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