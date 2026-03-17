# tests to check compatability with autodiff (specifically Turing / NUTS)

include("test00_setup.jl")

debkiss = SimplifiedEnergyBudget() |> instantiate

#@test debkiss.parameters |> eltype == Union{Real,Distribution} # the species-level parameters can be scalar or distributions
@test debkiss.parameters |> debkiss.generate_individual_params |> eltype |> x-> x <: Real # the individual-level parameters are all scalar

@time simulate(debkiss);

using ModelingToolkit

# TODO: use this info to deal with if/else in symbolis
#  https://docs.sciml.ai/ModelingToolkit/dev/basics/FAQ/#How-do-I-handle-if-statements-in-my-symbolic-forms?


debkiss.parameters

p_ind = debkiss.generate_individual_params(debkiss.parameters)
u0 = debkiss.initialize_all_statevars(p_ind)
tspan = (0,21)
p = p_ind
prob = ODEProblem(debkiss.complete_derivatives!, u0, tspan, p)

@mtkcompile sys = modelingtoolkitize(
    prob, 
    u_names = Symbol.(ComponentArrays.labels(u0)),
    )

prob_jac = ODEProblem(sys, u0, (0.,21.), jac = true)
sol = solve(prob_jac);
@time sol = solve(prob_jac);

hcat(sol.u...)

sol.t

prob_jac