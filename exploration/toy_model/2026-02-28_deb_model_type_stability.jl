using Pkg; Pkg.activate("test")
include(joinpath(pwd(), "test", "test00_setup.jl"))

# WIP: fix circular reference
# TODO: fix larger maxiters needed
# TODO: fix bad benchmark


debkiss = SimplifiedEnergyBudget() |> instantiate
sim = simulate_ode(debkiss)
include(joinpath(pwd(), "test", "test01_defaults.jl"))


using BenchmarkTools

debkiss.parameters.spc.Z = 1.
simulate_ode(debkiss)
@benchmark simulate_ode(debkiss)