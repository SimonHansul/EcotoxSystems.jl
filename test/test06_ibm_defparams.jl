using Pkg; Pkg.activate("test")
using Plots, StatsPlots
using Revise
import EnergyBudgetDiffEqs: params, ODE_simulator, IBM_simulator

p = params()
p.spc.H_p = 100.
@time sim_ode = ODE_simulator(p);
@df sim_ode plot(:t, [:S, :R, :H], layout = (1,3))

using Distributions

p.glb.Xdot_in
p.glb.k_V
p.spc.Z = Truncated(Normal(1, 0.05), 0, Inf)

sim_ibm = IBM_simulator(p);

DEB.DEBODE_global!
DEB.DEBODE_individual!
DEB.default_individual_rules!

methods(DEB.IndividualBasedModel)

