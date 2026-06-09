# ======================================== #
# Simulate temp effects
# ======================================== #

include("test00_setup.jl")

using EcotoxSystems.DEBkiss
import EcotoxSystems: relative_response
using EcotoxSystems
using DataFramesMeta

@testset "Temperature" begin # effect of food input
    norm(x) = x ./ sum(x)
    # prepare the plot
    plt = plot(
        layout = grid(1,3, widths = norm([2, 1, 1])),
        leg = false, 
        title = ["Growth" "Reproduction" "Food density"], 
        leftmargin = 5mm, bottommargin = 6mm, 
        size = (1200,350), 
        xlabel = "Time (d)"
        )

    debkiss = FullDEBkiss() 
    p = debkiss.parameters
    p.glb.T = 293.15
    p.glb.t_max = 56.
    p.glb.dX_in = 1e10 

    p.spc.K_X = 12e3

    global sim = DataFrame()

    # iterate over temperatures
    let T_degC = 10. # temperature in °C
        for _ in 1:4
            T_degC += 5

            p.glb.T = T_degC + 273.15
  
            # generate the predidction
            @time sim_i = simulate(debkiss; reltol = 1e-4, alg = Tsit5())

            # plot the trajectories
            @df sim_i plot!(plt, :t, :S, ylabel = "S", subplot = 1, leg = :outertopleft, label = "T = $(T_degC)") 
            @df sim_i plot!(plt, :t, :R, ylabel = "R", subplot = 2)
            @df sim_i plot!(plt, :t, :X ./p.glb.V_patch, ylabel = "[X]", subplot = 3, 
                yscale = :log10
                )

            sim_i[!,:T_degC] .= T_degC
            append!(sim, sim_i)
        end
        hline!(
            plt, 
            [DEBkiss.calc_S_max(p.spc)], 
            linestyle = :dash, 
            color = "gray", 
            subplot = 1, 
            label = "S_max"
            )
        display(plt)
    end

    # verify that we have simulated a case without food competition - this would change the validty of the following tests

    fX = DEBkiss.f_X.(sim.X, p.glb.V_patch, p.spc.K_X)
    @test  minimum(fX) .> 0.99
    
    # verify that higher temperature == faster growth
    perfect_rankcorr = combine(groupby(sim, :T_degC)) do df 
        
        df_t = @subset(df, )
        a = sort(df, :T_degC)[:,[:T_degC,:S]]
        b = sort(df, :S)[:,[:T_degC,:S]]

        return a == b
    end

    @test unique(perfect_rankcorr.x1) == [true]
end