@testset"food availability" begin 

    debkiss = SimplifiedEnergyBudget() |> instantiate
    p = debkiss.parameters

    norm(x) = x ./ sum(x)
    # prepare the plot
    plt = plot(
        layout = grid(1,4, widths = norm([2, 1, 1, 1])),
        leg = false, 
        title = ["Growth" "Reproduction" "Food density" "Functional response"], 
        leftmargin = 5mm, bottommargin = 6mm, 
        size = (1500,350), 
        xlabel = "Time (d)"
        )
    
    global sim = DataFrame()
    p = EcotoxSystems.debkiss_defaultparams |> deepcopy
    # iterate over nutrient input concentrations
    let dX_in = 4800.
        for i in 1:5
            dX_in /= 2
            # generate the predidction
            p.glb.dX_in = dX_in
            p.glb.t_max = 56.
            p.spc.K_X = 12e3
            @time sim_i = simulate(debkiss, saveat = 1);

            @test nrow(sim_i) == 57

            # plot the trajectories
            @df sim_i plot!(plt, :t, :S, ylabel = "S", subplot = 1, leg = :outertopleft, label = "dX_in = $(dX_in)") 
            @df sim_i plot!(plt, :t, :R, ylabel = "R", subplot = 2)
            @df sim_i plot!(plt, :t, :X ./ p.glb.V_patch, ylabel = "[X]", subplot = 3, 
                yscale = :log10
                )
            @df sim_i plot!(plt, :t, :f_X, ylabel = "f", subplot = 4)

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
