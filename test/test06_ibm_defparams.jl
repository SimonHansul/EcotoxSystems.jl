using Pkg; Pkg.activate("test")
using Plots, StatsPlots
using Revise
import EnergyBudgetDiffEqs: params, ODE_simulator, IBM_simulator

p = params()
p.spc.H_p = 100.
@time sim_ode = ODE_simulator(p);
@df sim_ode plot(:t, [:S, :R, :H], layout = (1,3))

using Distributions


# FIXME: there's no food

p.glb.Xdot_in = 1e10
p.glb.k_V
p.spc.Z = Truncated(Normal(1, 0.05), 0, Inf)
p.spc.tau_R = 2.
p.spc.eta_AR

sim_ibm = IBM_simulator(p)

@df sim_ibm plot(
    plot(:t, :S), 
    plot(:t, :f_X)
)


(glb = (
    X_p = -21.667178163038983, 
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


import EnergyBudgetDiffEqs: sig

smin = 0.5
f_Xthr = 0.25

x = 0.25

plot(x -> sig(x, f_Xthr, (1 - smin) * x / f_Xthr + smin, 1, beta = 1e3), xlim = (0,1), ylim = (0,1))