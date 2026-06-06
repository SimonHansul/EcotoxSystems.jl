include("test00_setup.jl")

using EcotoxSystems.DEBkiss

debkiss = FullDEBkiss()
sim = DEBkiss.sim_all(debkiss.parameters)

@df sim plot(
    plot(:t, :S), 
    plot(:t, :H)
)



begin # setting up a debkiss instance
    debkiss = FullDEBkiss() |> instantiate
    debkiss.parameters.glb.t_max = 56.  
end

@testset "Default debkiss parameters" begin
    p = debkiss.parameters  
    p.spc.Z = Dirac(1.)  
    @time global sim = simulate(debkiss; reltol = 1e-3)
    @test isapprox(maximum(sim.H), p.spc.H_p, rtol = 0.1) 
    @test isapprox(maximum(sim.S), DEBkiss.calc_S_max(p.spc), rtol = 0.1)
end;

#@time sim = simulate(debkiss); # bencmark around 2.5 ms is decent...✅
#using BenchmarkTools; 
#@benchmark simulate(debkiss)

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
    sim = replicates(p->(DEBkiss.sim_all(p)), debkiss.parameters, 10)

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