#=
Testing the default parameters

- is the maximum maturity equal to maturity at puberty?
- does the maximum structure match with what is calculated from parameters?
- does the mass balance check out?
=#

@testset begin  
    p = defaultparams
    p.glb.t_max = 56.
    p.spc.Z = Dirac(1.)
    global sim = DEBODE.simulator(p, reltol = 1e-3)


    @test isapprox(maximum(sim.H), p.spc.H_p_0, rtol = 0.1) # test for maximum maturity
    @test isapprox(maximum(sim.S), DEBODE.calc_S_max(p.spc), rtol = 0.1)
end;


#=
Basic test of @replicates macro
=#

@testset begin
    p = defaultparams
    p.spc.Z = Truncated(Normal(1., 0.1), 0, Inf)
    sim = @replicates DEBODE.simulator(p) 10

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

    @test cvs.S > 0.05 # test for plausible coefficients of variation in the final values
    @test cvs.H > 0.05
    @test cvs.R > 0.05
end;
