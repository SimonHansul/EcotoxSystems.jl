#=
Simulate single stressors with different PMoAs
=#

using Pkg; Pkg.activate("test")

using Test
using Distributions

using OrdinaryDiffEq
using DataFrames, DataFramesMeta
using Plots, StatsPlots, Plots.Measures
default(leg = false)
using StatsBase

using Revise
@time using EcotoxSystems

import EcotoxSystems: exposure, relative_response
 
# NOTE: this will currently fail, because we assume that glb.C_W is always called as a function
@testset "Constant exposure" begin

    global p = EcotoxSystems.params()
    
    p.glb.t_max = 42.
    p.glb.dX_in = 2400.
    p.glb.k_V = 0.

    p.spc.H_p = 0.75

    p.spc.KD .= 0
    p.spc.E .= 1.
    p.spc.B .= 2.

    let C_Wvec =  vcat(0., round.(10 .^ range(log10(1.01), log10(10.), length = 5), sigdigits = 2)) #hcat([0], round.(10 .^ range(log10(1.01), log10(10.), length = 5), sigdigits = 2)...)' |> Matrix

        global sims = DataFrame()
        global pmoas = ["G", "M", "A", "R"]
        
        for (j,pmoa) in enumerate(pmoas)
            
            p.spc.KD .= 0.
            p.spc.KD[1,j] = 1.
            p.spc.KD_h[1] = 1.

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
        @test unique(rankcor.r .<= 0.9) == [true] # all pmoas affect reproduction
        @test minimum(sims.y_R) .< 0.5
        @test minimum(sims[sims.pmoa.=="R",:].y_S) > 0.99 # no effect on growth for PMoA R
    end
end

p.glb.C_W = 

using DataFrames 

exposure_profiles = DataFrame(
    treatment_id = repeat([1], 12),
    t     = Vector{Float64}([0, 2, 5, 7, 7,    9,    11,      14,   14,   16,   19,   21]), 
    C_W_1 = Vector{Float64}([0, 0, 0, 0, 0.1,  0.05, 0.025,   0.,   0.2,  0.1,  0.05,  0.025])
)

@df exposure_profiles scatter(:t, :C_W_1)

EcotoxSystems.preprocess_exposure_input!(p)

using Interpolations

xs = exposure_profiles.t
A = exposure_profiles[:,3]

Interpolations.deduplicate_knots!(xs)

linear_interpolator = linear_interpolation(xs, A);
linear_interpolator(22.)

linear_interpolator = linear_interpolation([0., Inf], [0., 0.])
linear_interpolator(rand())

let xvalues = 0:0.01:21
    plot!(xvalues, linear_interpolator.(xvalues))
end




@testset "Time-variable exposure" begin
    
    p.glb.C_W = 
    exposure_simulator()
    

end

