using Pkg; Pkg.activate("test")
include(joinpath(pwd(), "test", "test00_setup.jl"))



debkiss = SimplifiedEnergyBudget() |> instantiate
sim = simulate_ode(debkiss)
include(joinpath(pwd(), "test", "test01_defaults.jl"))


using BenchmarkTools

debkiss.parameters.spc.Z = 1.
simulate_ode(debkiss, alg = Rodas5())
@benchmark simulate_ode(debkiss, alg = Rodas5())


# we will probably need something like this when using Turing.jl efficientyl
# function getprob(ddeb)
# return ODEProblem

# function update_prob_params!(...)