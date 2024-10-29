using Plots.Measures
using StatsBase

@testset begin # effect of food input
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
    p = DEB.params()
    # iterate over nutrient input concentrations
    let Xdot_in = 4800.
        for i in 1:5
            Xdot_in /= 2
            # generate the predidction
            p.glb.Xdot_in = Xdot_in
            p.glb.t_max = 56.
            p.spc.K_X = 12e3
            sim_i = DEB.simulator(p, reltol = 1e-10)

            # plot the trajectories
            @df sim_i plot!(plt, :t, :S, ylabel = "S", subplot = 1, leg = :outertopleft, label = "Xdot_in = $(Xdot_in)") 
            @df sim_i plot!(plt, :t, :R, ylabel = "R", subplot = 2)
            @df sim_i plot!(plt, :t, :X_p ./ p.glb.V_patch, ylabel = "[X_p]", subplot = 3, 
                yscale = :log10
                )

            sim_i[!,:Xdot_in] .= Xdot_in 
            append!(sim, sim_i)
        end
        hline!(plt, [DEB.calc_S_max(p.spc)], linestyle = :dash, color = "gray", subplot = 1, label = "S_max")
        display(plt)
    end

    # checking that maximum size increases strictly monotonically with increasing food availability
    rankcor_size = combine(groupby(sim, :Xdot_in), :S => maximum) |> x -> corspearman(x.Xdot_in, x.S_maximum)

    @test rankcor_size == 1 
end
