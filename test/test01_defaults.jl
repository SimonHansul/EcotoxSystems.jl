@testset "Default parameters" begin  

    global p = EcotoxSystems.params()
    p.glb.t_max = 56.
    p.spc.Z = Dirac(1.)
    @time global sim = EcotoxSystems.ODE_simulator(p, reltol = 1e-3)

    @test isapprox(maximum(sim.H), p.spc.H_p, rtol = 0.1) 
    @test isapprox(maximum(sim.S), EcotoxSystems.calc_S_max(p.spc), rtol = 0.1)
end;

#=
 inital allocations: 50k (>1 mb)
   changed signature of sig and clipneg from Real to Float64
   did not change much 
   adding @. to du.ind.D_z, despite already looping explicitly?
       does not work
   move TK to separate inlined function with Float64 signature?
       maybe improved things a little (50k -> 43k allocations)
   did the same with caluclation of y_T
       some more, small improvement (43k -> 40k)
   same with f_X
       allocations down to 38k...(little less than 1MB)
   changed TKTD_mix_IA to apply product to y_j in inner loop
  in TKTD_mix_IA, removed call to reduce in favor of updating y_j and h_z within the loop
       allocations down from 50k to 5k! (≈ 300 kB)
       runtime down to 1-2 ms 
       no more suspicious in ODE part for now
  switching from callbacks to sig made it slightly worse
        9k allocations, 500 kB
=#

@time sim = EcotoxSystems.ODE_simulator(p);
#VSCodeServer.@profview_allocs [EcotoxSystems.ODE_simulator(p) for _ in 1:100];

#=
Basic test of @replicates macro
=#

import EcotoxSystems: @replicates
using Chain

@testset "@replicates macro" begin
    p = EcotoxSystems.params()
    p.spc.Z = Truncated(Normal(1., 0.1), 0, Inf)
    sim = @replicates EcotoxSystems.ODE_simulator(p) 10

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


@testset "replicates function" begin
    
    p = EcotoxSystems.params()
    p.spc.Z = Truncated(Normal(1., 0.1), 0, Inf)
    sim = replicates(EcotoxSystems.ODE_simulator, p, 10)

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