@testset "Default parameters" begin  

    global p = EcotoxSystems.params()
    p.glb.t_max = 56.
    p.spc.Z = Dirac(1.)
    @time global sim = EcotoxSystems.ODE_simulator(p, reltol = 1e-3)

    @test isapprox(maximum(sim.H), p.spc.H_p, rtol = 0.1) 
    @test isapprox(maximum(sim.S), EcotoxSystems.calc_S_max(p.spc), rtol = 0.1)
end;

# inital allocations: 50k (>1 mb)
#   changed signature of sig and clipneg from Real to Float64
#   did not change much 
#   adding @. to du.ind.D_z, despite already looping explicitly?
#       does not work
#   move TK to separate inlined function with Float64 signature?
#       maybe improved things a little (50k -> 43k allocations)
#   did the same with caluclation of y_T
#       some more, small improvement (43k -> 40k)
#   same with f_X
#       allocations down to 38k...(little less than 1MB)
#   changed TKTD_mix_IA to apply product to y_j in inner loop

@time EcotoxSystems.ODE_simulator(p);
VSCodeServer.@profview_allocs [EcotoxSystems.ODE_simulator(p) for _ in 1:10];

Base.summarysize(sim)


using BenchmarkTools
x = 1
@benchmarks x = EcotoxSystems.sig(
    1.,
    0.,
    0.,
    1.
);

using BenchmarkTools
let m = [1, 2, 3], x = 0
    @info "Benchmarking slicing"
    @benchmark $x + $m[1]
end

let m = [1, 2, 3], x = 0
    @info "Benchmarking version with @view"
    @benchmark @. $x + @view($m[1])
end


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