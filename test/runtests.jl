using Pkg; Pkg.activate("test")
using Revise
@time using EnergyBudgetDiffEqs


using OrdinaryDiffEq
using DataFrames, ComponentArrays
using Plots, StatsPlots
default(leg = false)

p = EnergyBudgetDiffEqs.params()
p.glb.Xdot_in = 1e10

p.glb.T 
p.spc.T_ref
p.glb.t_max = 50.
p.spc.Idot_max_rel_emb

spc = p.spc 

@benchmark sim = EnergyBudgetDiffEqs.simulator(p)

@time EnergyBudgetDiffEqs.simulator(p, returntype=EnergyBudgetDiffEqs.odesol);

VSCodeServer.@profview_allocs EnergyBudgetDiffEqs.simulator(p, returntype=EnergyBudgetDiffEqs.odesol);

using BenchmarkTools

@df sim plot(
    plot(:t, :X_emb),
    plot(:t, :f_X, ylim = (-0.1, 1.1)),
    plot(:t, :S),
)

zprop = ind.Z .^ ind.propagate_zoom




spc[keys(spc.propagate_zoom)] .* 1 .^ spc.propagate_zoom

ind = ComponentVector(
    ind;
    zprop...
)