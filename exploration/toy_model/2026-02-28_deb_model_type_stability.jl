using Pkg; Pkg.activate("test")
include(joinpath(pwd(), "test", "test00_setup.jl"))

sim = simulate_ode(debkiss)
include(joinpath(pwd(), "test", "test01_defaults.jl"))

using BenchmarkTools

debkiss = SimplifiedEnergyBudget() |> instantiate
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

VSCodeServer.@profview_allocs simulate_ode(debkiss, alg = Rodas5(), returntype = EcotoxSystems.odesol)
prob = simulate_ode(debkiss, alg = Rodas5(), returntype = EcotoxSystems.odeprob)
integrator = init(prob, Rodas5())

du = similar(integrator.u)
reset_du(du) = du .= 0.;
@code_warntype reset_du(du)

t = 0.

@code_warntype debkiss.complete_derivatives!(du, integrator.u, integrator.p, t)
@allocated debkiss.complete_derivatives!(du, integrator.u, integrator.p, t)


# we will probably need something like this when using Turing.jl efficientyl
# function getprob(ddeb)
# return ODEProblem

# function update_prob_params!(...)