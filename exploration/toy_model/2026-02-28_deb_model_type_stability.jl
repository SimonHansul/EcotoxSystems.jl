using Pkg; Pkg.activate("test")
include(joinpath(pwd(), "test", "test00_setup.jl"))

debkiss = SimplifiedEnergyBudget() |> instantiate
sim = simulate_ode(debkiss)
include(joinpath(pwd(), "test", "test01_defaults.jl"))

using BenchmarkTools

debkiss.parameters.spc.Z = 1.
sim = simulate_ode(debkiss, alg = Rodas5())

debkiss.parameters.spc.kappa

@df sim plot(
    plot(:t, :S), 
    plot(:t, :H), 
    plot(:t, :R), 
    plot(:t, [:is_embryo :is_adult])
)

debkiss.parameters.spc.H_p

@benchmark simulate_ode(debkiss, alg = Rodas5())


# we will probably need something like this when using Turing.jl efficientyl
# function getprob(ddeb)
# return ODEProblem

# function update_prob_params!(...)