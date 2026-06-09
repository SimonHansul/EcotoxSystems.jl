#=
Simulate single stressors with different PMoAs
=#

include("test00_setup.jl")

using EcotoxSystems.DEBkiss
import EcotoxSystems: relative_response
using EcotoxSystems
using DataFramesMeta

@testset begin

    debkiss = FullDEBkiss()
    p = debkiss.parameters
    tktd_ref = deepcopy(debkiss.parameters.spc.TKTD)
    
    p.glb.t_max = 42.
    p.glb.dX_in = 2400.
    p.glb.k_V = 0.
    p.spc.H_p = 75.

    let Cvec =  vcat(0., round.(10 .^ range(log10(1.01), log10(10.), length = 5), sigdigits = 2)) 

        global sims = DataFrame()
        global pmoas = ["G", "M", "A", "R"]
        
        for (j,pmoa) in enumerate(pmoas)

            tktd = deepcopy(tktd_ref)
            DEBkiss.set_TKTD_param!(tktd, :k_D, 1, pmoa, 1.)
            DEBkiss.set_TKTD_param!(tktd, :e, 1, pmoa, 1.)
            DEBkiss.set_TKTD_param!(tktd, :b, 1, pmoa, 1.)

            p.spc.TKTD = tktd

            # using the exposure function to iterate over treatments
            sim_j = DEBkiss.simulate_constant_exposure(p, :C_W1 => Cvec; saveat = 1)
            sim_j[!,:pmoa] .= pmoa
            append!(sims, sim_j)
        end 
        
        sims.t = ceil.(sims.t)
        # FIXME: relative_response creates artefacts when used with life stage-concatenating simulators (like in DEBkiss.sim_all)
        sims = relative_response(sims, [:S, :R], :treatment_id; groupby_vars = [:t, :pmoa])
        sort!(sims, :t)
        
        rankcor = combine(groupby(sims, :pmoa)) do sim_j
            @chain sim_j begin
                combine(groupby(_, :C_W1), :y_R => last) 
                corspearman(_.C_W1, _.y_R_last)
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
            @df @subset(sims, :pmoa .== pmoa) plot!(plt, :t, :y_S, group = :treatment_id, subplot = j, label = hcat(Cvec...))
            @df @subset(sims, :pmoa .== pmoa) plot!(plt, :t, :y_R, group = :treatment_id, subplot = 4+j, label = hcat(Cvec...))
        end
        
        display(plt)

        sims_ad = @subset(sims, :is_adult .== 1)
        
        @test unique(rankcor.r .<= 0.9) == [true] # all pmoas affect reproduction
        @test minimum(sims_ad[sims_ad.pmoa.=="R",:].y_S) > 0.8 # no effect on growth for PMoA R

    end
end