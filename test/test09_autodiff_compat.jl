# tests to check compatability with autodiff (specifically Turing / NUTS)

include("test00_setup.jl")

debkiss = SimplifiedEnergyBudget() |> instantiate

debkiss.parameters |> eltype
debkiss.parameters |> generate_individual_params |> eltype