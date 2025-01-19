using Plots.Measures
using StatsBase

@testset"food availability" begin 
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
    
    sim = DataFrame()
    p = EcotoxSystems.params()
    # iterate over nutrient input concentrations
    let dX_in = 4800.
        for i in 1:5
            dX_in /= 2
            # generate the predidction
            p.glb.dX_in = dX_in
            p.glb.t_max = 56.
            p.spc.K_X = 12e3
            sim_i = EcotoxSystems.ODE_simulator(p, reltol = 1e-10)

            # plot the trajectories
            @df sim_i plot!(plt, :t, :S, ylabel = "S", subplot = 1, leg = :outertopleft, label = "dX_in = $(dX_in)") 
            @df sim_i plot!(plt, :t, :R, ylabel = "R", subplot = 2)
            @df sim_i plot!(plt, :t, :X ./ p.glb.V_patch, ylabel = "[X]", subplot = 3, 
                yscale = :log10
                )

            sim_i[!,:dX_in] .= dX_in 
            append!(sim, sim_i)
        end
        hline!(plt, [EcotoxSystems.calc_S_max(p.spc)], linestyle = :dash, color = "gray", subplot = 1, label = "S_max")
        display(plt)
    end

    # checking that maximum size increases strictly monotonically with increasing food availability
    rankcor_size = combine(groupby(sim, :dX_in), :S => maximum) |> x -> corspearman(x.dX_in, x.S_maximum)

    @test rankcor_size == 1 
end
