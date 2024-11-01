#using Pkg; Pkg.activate("test")
#using Plots, StatsPlots
#using Revise
@time using EcotoxSystems
@time import EcotoxSystems: params, ODE_simulator, IBM_simulator 
@time import EcotoxSystems: @replicates, DEBIndividual, treplicates

@time using DataFrames, DataFramesMeta

using DrWatson
using PackageCompiler

@time using Plots
default(leg = false)

using Distributions
using DataFramesMeta
using DataFrames

begin
    p = params()

    p.glb.dX_in = 1e3
    p.glb.k_V = 0.1
    p.glb.V_patch = 0.05
    p.glb.N0 = 1.
    p.glb.t_max = 21

    p.spc.Z = Truncated(Normal(1, 0.05), 0, Inf)
    p.spc.tau_R = Inf # don't simulate reproduction
    p.spc.f_Xthr = 0.9
    p.spc.H_p = 100.
    p.spc.a_max = Inf


    sim_ode = ODE_simulator(p);
    @time sim_ibm = IBM_simulator(p, dt = 1/24, showinfo = 7) 
end

# TODO: convert comparison between ODE and IBM simulator to a proper test
# there is some divergence but it looks o.k.

@df sim_ode plot(:t, :S, xlabel = "t", ylabel = "Resource abundance", leg = true, label = "ODE simulator")
@df sim_ibm.spc plot!(:t, :S, label = "IBM simulator")

@df sim_ode plot(
    :t, :X, 
    xlabel = "t", ylabel = "Resource abundance", 
    leg = true, label = "ODE simulator"
    )
@df sim_ibm.glb plot!(:t, :X, label = "IBM simulator")


@testset "Running IBM" begin
    # FIXME: lots of memory allocs in the ODE part
    # probably a whoopsie in how I used broadcasting in the TKTD model
    # (cf https://discourse.julialang.org/t/can-i-avoid-allocations-when-broadcasting-over-slices/102501/2)
    # keyword "broadcast fusion" https://bkamins.github.io/julialang/2023/03/31/broadcast.html

    # we can play a little with the parameters here, but we are not looking for a "correct" prediction!

    p.glb.dX_in = 30_000 #2e4 * 0.66 * 0.5
    p.glb.k_V = 0.1
    p.glb.V_patch = 0.5
    p.glb.N0 = 10
    p.glb.t_max = 56

    p.spc.Z = Truncated(Normal(1, 0.05), 0, Inf)
    p.spc.tau_R = 2.
    p.spc.eta_IA = 0.33
    
    p.spc.f_Xthr = 0.75
    p.spc.s_min = 0.25

    @time sim_ibm = treplicates(x -> IBM_simulator(p, saveat = 1, showinfo = 7), p, 4)
    
   
    plot(
        (@df sim_ibm.glb plot(:t, :N, group = :replicate)),
        #(@df sim_ibm.spc plot(:t, :S, group = :replicate .* :id)), 
        #(@df sim_ibm.spc plot(:t, :S, group = :replicate .* :id)), 
        (@df sim_ibm.spc plot(:t, :f_X, group = :replicate .* :id)), 
        leg = false
    )
end


# trying the simulation without TKTD to compare exec time

function simple!(du, u, p, t) 
    # turning off all possible toxic effects

    u.ind.y_j .= 1
    u.ind.h_z = 0.

    EcotoxSystems.DEBkiss!(du, u, p, t)
end

p.glb.t_max = 56.

@time sim_ibm = IBM_simulator(p, individual_ode! = simple!, showinfo = 7);
VSCodeServer.@profview IBM_simulator(p, individual_ode! = simple!)

popsize = combine(
    groupby(sim_ibm.spc, [:t]), 
    x -> (N = nrow(x), M = sum(x.S) + sum(x.R))
    )

@df popsize plot(:t, :N)

@df sim_ibm.spc plot(
    plot(:t, :I, group = :id),
    plot(:t, :S, group = :id), 
    plot(:t, :f_X, group = :id),
    leg = false
)

