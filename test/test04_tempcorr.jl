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

    p = DEB.params()
    p.glb.T = 293.15
    p.glb.t_max = 56.
    p.spc.K_X = 12e3
    
    sim = DataFrame()
    # iterate over nutrient input concentrations
    let T_C = 10. # temperature in °C
        for _ in 1:4
            T_C += 5

            p.glb.T = T_C + 273.15

            # generate the predidction
            sim_i = DEB.simulator(p)

            # plot the trajectories
            @df sim_i plot!(plt, :t, :S, ylabel = "S", subplot = 1, leg = :outertopleft, label = "T = $(T_C)") 
            @df sim_i plot!(plt, :t, :R, ylabel = "R", subplot = 2)
            @df sim_i plot!(plt, :t, :X_p ./ p.glb.V_patch, ylabel = "[X_p]", subplot = 3, 
                yscale = :log10
                )

            sim_i[!,:T_C] .= T_C
            append!(sim, sim_i)
        end
        hline!(plt, [DEB.calc_S_max(p.spc)], linestyle = :dash, color = "gray", subplot = 1, label = "S_max")
        display(plt)
    end

    rankcor = combine(groupby(sim, :T_C), :S => maximum) |> x -> corspearman(x.T_C, x.S_maximum)
    
    @test rankcor == -1 
end


