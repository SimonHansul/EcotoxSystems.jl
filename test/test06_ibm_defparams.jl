using Pkg; Pkg.activate("test")
#using Plots, StatsPlots

using Plots
default(leg = false)

using Distributions
using DataFramesMeta
using DataFrames
using Test

using Revise
@time using EcotoxSystems
@time import EcotoxSystems: params, ODE_simulator, IBM_simulator 
@time import EcotoxSystems: @replicates, DEBIndividual, treplicates

using Plots, StatsPlots
begin
    p = params()

    p.glb.dX_in = 1e4
    p.glb.k_V = 0.1
    p.glb.V_patch = 0.05
    p.glb.N0 = 1.
    p.glb.t_max = 56

    p.spc.Z = Dirac(1) # Truncated(Normal(1, 0.05), 0, Inf)
    p.spc.tau_R = Inf # don't simulate reproduction for now

    p.spc.H_p = 100.
    p.spc.a_max = Inf

    @time sim_ode = ODE_simulator(p, saveat = 1);
    @time sim_ibm = IBM_simulator(p, dt = 1/24, showinfo = Inf);

    @df sim_ode plot(:t, :S, xlabel = "t", ylabel = "S", leg = true, label = "ODE simulator")
    @df sim_ibm.spc plot!(:t, :S, label = "IBM simulator")
end

# TODO: convert comparison between ODE and IBM simulator to a proper test

@testset "Running IBM" begin
    # FIXME: lots of memory allocs in the ODE part
    # maybe a problem in how I used broadcasting in the TKTD model, although removing TKTD part does not change so much
    # (cf https://discourse.julialang.org/t/can-i-avoid-allocations-when-broadcasting-over-slices/102501/2)
    # keyword "broadcast fusion" https://bkamins.github.io/julialang/2023/03/31/broadcast.html

    # we can play a little with the parameters here, but we are not looking for a "correct" prediction here!
    # the IBM_simulator just provides the infrastructure to combine the necessary components

    p.glb.dX_in = 50_000 #100_000
    p.glb.k_V = 0.1
    p.glb.V_patch = 0.5
    p.glb.N0 = 10
    p.glb.t_max = 56 #365

    p.spc.Z = Truncated(Normal(1, 0.1), 0, Inf)
    p.spc.tau_R = truncated(Normal(2., 0.2), 0, Inf)

    @time sim_ibm = treplicates(x -> IBM_simulator(p, saveat = 1, showinfo = 14), p, 1)
    #VSCodeServer.@profview_allocs global sim_ibm = treplicates(x -> IBM_simulator(p, saveat = 1, showinfo = 14), p, 1)
    
    plt = plot(
        (@df sim_ibm.glb plot(:t, :N, group = :replicate, ylabel = "N")),
        (@df sim_ibm.spc groupedlineplot(:t, :S, :cohort, ylabel = "S")), 
        (@df sim_ibm.spc groupedlineplot(:t, :H, :cohort, ylabel = "H")), 
        leg = false
    )

    display(plt)
end

