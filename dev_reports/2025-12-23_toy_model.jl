#=
Here I implement a toy model that is as simple as I can think of in order to

1. demonstrate basic usage of the package 
2. make sense of unexpected performance drops that I currently see for more complex models (2025-12-23)
=#

using Pkg; Pkg.activate("test")
include(joinpath("..", "test", "test00_setup.jl"))

using Glob


include(joinpath("toy_model", "derivatives.jl"))  
include(joinpath("toy_model", "statevars.jl"))  
include(joinpath("toy_model", "parameters.jl"))  


begin # having a lookg at the ODE output

    toy_parameters.glb.t_max = 30.
    toy_parameters.glb.N0 = 1.

    sim_ode = EcotoxSystems.ODE_simulator(
        toy_parameters; 
        model = toy_ODE!, 
        statevars_init = toy_statevars, 
        callback = adult_callback, 
        gen_ind_params = EcotoxSystems.generate_individual_params_nozoom
    )

    @df sim_ode plot(
        plot(:t, :L, ylabel = "length"), 
        plot(:t, :R, ylabel = "repro buffer"), 
        xlabel = "age"
    )

    toy_parameters.glb.N0 = 5

    sim_ode = EcotoxSystems.ODE_simulator(
        toy_parameters; 
        model = toy_ODE!, 
        statevars_init = toy_statevars, 
        callback = adult_callback, 
        gen_ind_params = EcotoxSystems.generate_individual_params_nozoom
    )

    @df sim_ode plot!(:t, :L, subplot = 1)
    @df sim_ode plot!(:t, :R, subplot = 2)

end