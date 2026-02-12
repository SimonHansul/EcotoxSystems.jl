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

    debkiss = SimplifiedEnergyBudget() |> instantiate
    p = debkiss.parameters
    p.glb.T = 293.15
    p.glb.t_max = 56.
    p.glb.dX_in = 1e10 

    p.spc.K_X = 12e3

    sim = DataFrame()

    # iterate over nutrient input concentrations
    let T_degC = 10. # temperature in °C
        for _ in 1:4
            T_degC += 5

            p.glb.T = T_degC + 273.15
  
            # generate the predidction
            sim_i = simulate(debkiss)

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
            [EcotoxSystems.calc_S_max(p.spc)], 
            linestyle = :dash, 
            color = "gray", 
            subplot = 1, 
            label = "S_max"
            )
        display(plt)
    end

    # check that we have simulated a case without food competition - this would change the validty of the following tests

    cv_X_sum = 0. # coefficient of variation in food abundance X across temp treatments, summed over time points

    for t in unique(sim.t)
        sim_rand = @subset(sim, :t .== t)
        cv_X = sim_rand.X |> x -> std(x)/mean(x)
        cv_X_sum += cv_X
    end

    cv_X_mean = cv_X_sum/(length(unique(sim.t))) # average CVs

    @test cv_X_mean < 1e-3 # verify that CV is very small
    
    # check that higher temperature == faster growth

    rankcor = combine(groupby(sim, :T_degC), :S => maximum) |> 
    x -> corspearman(x.T_degC, x.S_maximum)
    
    @test rankcor == 1 

    # check that final mass is barely affected

    sim_end = @subset(sim, :t .== maximum(:t))
    cv_S = sim_end.S |> x -> std(x) / mean(x)

    @test cv_S < 0.01

    
end


