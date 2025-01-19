#=
Simulate single stressors with different PMoAs
=#

using DataFramesMeta
import EcotoxSystems: exposure, relative_response

begin
    p = EcotoxSystems.params()
    
    p.glb.t_max = 42.
    p.glb.dX_in = 2400.
    p.glb.k_V = 0.

    p.spc.k_D_z .= 0
    p.spc.e_z .= 1.
    p.spc.b_z .= 0.1

    let C_Wvec =  hcat([0], round.(10 .^ range(log10(1.01), log10(10.), length = 5), sigdigits = 2)...)' |> Matrix
        global sims = DataFrame()
        pmoas = ["G", "M", "A", "R"]
        for (j,pmoa) in enumerate(pmoas)

            p.spc.k_D_z .= 0.
            p.spc.k_D_z[1,j] = 1.
            p.spc.k_D_h[1] = 1.

            # using the exposure function to iterate over treatments
            sim_j = exposure(EcotoxSystems.ODE_simulator, p, C_Wvec)
            sim_j[!,:pmoa] .= pmoa
            append!(sims, sim_j)
        end 
        
        global sims = relative_response(sims, [:S, :R], :C_W_1; groupby_vars = [:t, :pmoa])

        sort!(sims, :t)
        
        rankcor = combine(groupby(sims, :pmoa)) do sim_j
            @chain sim_j begin
                combine(groupby(_, :C_W_1), :y_R => last) 
                corspearman(_.C_W_1, _.y_R_last)
                return DataFrame(r = _)
            end 
        end

        plt = plot(
            layout = (2,4), title = hcat(pmoas...), 
            ylim = (0, 1.01),
            xlabel = "t", 
            ylabel = "relative response", 
            size = (800,450), 
            bottommargin = 5mm, leftmargin = 2.5mm, lw = 2,
            leg = [false false false true false false false false],
            legtitle = "C_W / e", legendtitlefontsize = 10
            )
            
        for (j,pmoa) in enumerate(pmoas)
            @df @subset(sims, :pmoa .== pmoa) plot!(plt, :t, :y_S, group = :C_W_1, subplot = j, label = hcat(unique(:C_W_1)...))
            @df @subset(sims, :pmoa .== pmoa) plot!(plt, :t, :y_R, group = :C_W_1, subplot = 4+j, label = hcat(unique(:C_W_1)...))
        end
        
        display(plt)
        @test unique(rankcor.r .<= 0.9) == [true]
        @test minimum(sims.y_R) .< 0.5
    end
end