#=
Simulate single stressors with different PMoAs
=#

using Plots.Measures

using MechanisticEffectModels.ParamStructs
using MechanisticEffectModels.DEBODE

ENV["JULIA_DEBUG"] = Main

params = Params()
isolate_pmoas!(params, ["M"])
params.spc.k_D_M = [1.]
params.spc.e_M = [1.]
params.spc.b_M = [2.]


let C_Wvec =  vcat([0], round.(10 .^ range(log10(1.01), log10(10.), length = 5), sigdigits = 2))
    global sim = DataFrame()
    pmoas = ["G", "M", "A", "R"]
    for pmoa in pmoas
        for C_W in C_Wvec
            glb = GlobalParams(t_max = 42., C_W = [C_W])
            spc = SpeciesParams(
                kappa_0 = 0.538,
                k_D_G = [0.1], 
                k_D_M = [0.1], 
                k_D_A = [0.1], 
                k_D_R = [0.1], 
                k_D_h = [0.1], 
                e_G = [1.],
                e_M = [1.],
                e_A = [1.],
                e_R = [1.],
                e_h = [1.],
                b_G = [0.5],
                b_M = [0.5],
                b_A = [0.5],
                b_R = [0.5],
                b_h = [0.5]
                )
            theta = Params(glb = glb, spc = spc)
            isolate_pmoas!(theta.spc, [pmoa])
            sim_zj = DEBODE.simulator(theta)
            sim_zj[!,:C_W] .= C_W
            sim_zj[!,:pmoa] .= pmoa
            append!(sim, sim_zj)
        end # for C_W in ...
    end # for pmoa in ...
    
    sim = relative_response(sim, [:S, :R], :C_W; groupby_vars = [:t, :pmoa])

    sort!(sim, :t)
    
    rankcor = combine(groupby(sim, :pmoa)) do sim_j
        @chain sim_j begin
            combine(groupby(_, :C_W_1), :y_R => last) 
            corspearman(_.C_W_1, _.y_R_last)
            return DataFrame(r = _)
        end 
    end

    plt = plot(
        layout = (2,4), title = hcat(pmoas...), 
        ylim = (0, 1.01),
        xlabel = "t", ylabel = hcat([gridylabel("y_S", 1, 4), gridylabel("y_R", 1, 4)]...), 
        size = (800,450), 
        bottommargin = 5mm, leftmargin = 2.5mm, lw = 2,
        leg = [false false false true false false false false],
        legtitle = "C_W / e", legendtitlefontsize = 10
        )
        
    for (j,pmoa) in enumerate(pmoas)
        @df @subset(sim, :pmoa .== pmoa) groupedlineplot!(plt, :t, :y_S, :C_W, subplot = j, label = hcat(unique(:C_W)...))
        @df @subset(sim, :pmoa .== pmoa) groupedlineplot!(plt, :t, :y_R, :C_W, subplot = 4+j, label = hcat(unique(:C_W)...))
    end
    
    display(plt)
    @test unique(rankcor.r .<= 0.9) == [true]
end

