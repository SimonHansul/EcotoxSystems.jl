using Pkg; Pkg.activate("test")
using Plots, StatsPlots
using Revise
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

    #p.glb.dX_in = 100_000
    #p.glb.k_V = 0.1
    #p.glb.V_patch = 2.
    #p.glb.N0 = 10
    #p.glb.t_max = 63
    #p.spc.Z = Truncated(Normal(1, 0.05), 0, Inf)
    #p.spc.tau_R = 2.
    #p.spc.f_Xthr = 0.9
    #p.spc.H_p = 100.


    @time sim_ode = ODE_simulator(p);
end

@df sim_ode plot(:t, [:S, :R, :H, :X], layout = (1,4), size = (800,350))




# TODO: how can I test whether the food feedback is correct?

# v1
# simulate a single individual
# in output:
#   calculate dI from I
#   calculate dX from X_pF
# --> global vars are currently not returned!
#   - implement recording of global satates

# v2 
# include X in the ode-vs-ibm comparison

# trying to store individuals in Memory instead of Vector
begin
    # FIXME: lots of memory allocs in the ODE part
    # (cf https://discourse.julialang.org/t/can-i-avoid-allocations-when-broadcasting-over-slices/102501/2)
    # keyword "broadcast fusion" https://bkamins.github.io/julialang/2023/03/31/broadcast.html

    p.glb.dX_in = 100_000
    p.glb.k_V = 0.1
    p.glb.V_patch = 2.
    p.glb.N0 = 10
    p.glb.t_max = 63

    p.spc.Z = Truncated(Normal(1, 0.05), 0, Inf)
    p.spc.tau_R = 2.
    
    p.spc.f_Xthr = 0.9

    sim_ibm = IBM_simulator(p, saveat = 1, showinfo = 14)

    
    popsize = combine(
        groupby(sim_ibm, [:t, :replicate]), 
        x -> (N = nrow(x), M = sum(x.S) + sum(x.R))
        )

    plot(
        (@df popsize plot(:t, :N, group = :replicate)),
        (@df popsize plot(:t, :M, group = :replicate)),
        #(@df sim_ibm plot(:t, :S, group = :replicate .* :id)), 
        #(@df sim_ibm plot(:t, :f_X, group = :replicate .* :id)), 
        leg = false
    )

end

@df sim_ibm plot(
    plot(:t, :S, group = :id), 
    plot(:t, :f_X, group = :id),
    leg = false
)


(glb = (
    X = -21.667178163038983, 
    C_W = [1.0e-323]), 
ind = 
(
    embryo = 0.0, 
    juvenile = 0.0, 
    adult = 0.0, 
    X_emb = -30.50278570500914, 
    S = 3.6198141420927143, 
    S_max_hist = 0.0, 
    H = 3.8312451146890707, 
    R = 0.0, 
    f_X = 0.0,
    I_emb = 0.0, 
    I_p = 0.0, 
    I = 30.50278570500914, 
    A = 10.065919282653017, 
    M = 0.9007628157340829, 
    J = 0.8091436746139695, 
    D_z = [0.0 0.0 0.0 0.0], 
    D_h = [0.0], 
    y_z = [0.0 0.0 0.0 0.0], 
    y_j = [0.0, 0.0, 0.0, 0.0], 
    y_T = 0.0, 
    h_z = 0.0, 
    S_z = -0.0, 
    h_fX = 0.0, 
    id = 0.0, 
    cohort = 0.0, 
    age = 0.0, 
    cause_of_death = 0.0, 
    time_since_last_repro = 0.0,
    cum_offspring = 0.0
    ))


import EcotoxSystems: sig

smin = 0.5
f_Xthr = 0.25

x = 0.25

plot(x -> sig(x, f_Xthr, (1 - smin) * x / f_Xthr + smin, 1, beta = 1e3), xlim = (0,1), ylim = (0,1))