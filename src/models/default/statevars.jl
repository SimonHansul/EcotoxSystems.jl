# default_statevars.jl 
# functions to generate default state variables as ComponentArrays
# state variables are differentiated into global and individual-level
# this differentiation is essential for the functionality of the generic IBM schedule, 
# and to allow for extensions with mulitple types of agents

const X_EMB_INT_REL = 0.001 # for the default DEB model, this global constant determines the initial structural mass, relative to the mass of an egg

# some helper functions to construct mutable static vectors and matrices

constrmvec(x::AbstractMatrix; fillval::Float64 = 0.) = MVector{size(x)[1],Float64}(fill(fillval, size(x)[1]))
constrmvec(x::AbstractVector; fillval::Float64 = 0.) = MVector{length(x),Float64}(fill(fillval, length(x)))
constrmmat(x::AbstractMatrix; fillval::Float64 = 0.) = MMatrix{size(x)...,Float64}(fill(fillval, size(x)))
constrmmat(x::AbstractMatrix, dims::Int64; fillval::Float64 = 0.) = MMatrix{size(x)[dims],1,Float64}(fill(fillval, size(x)[dims]))

"""
    initialize_individual_statevars(
        p::ComponentVector; 
        id = 1., 
        cohort = 0.)::ComponentVector

This function defines the individual-level state variables and their initial values for the default model. 

`p` is a parameter vector containing all components which are relevant for this simulation (including global state varaiables).

`id` is an individual's unique identifier in the IBM simulation. 

`cohort`is the index of the cohort an indvidiual belongs to in the IBM simulation (i.e. initial generation is cohort 0, 
their offsrping is cohort 1, etc.)

To add state variables or overwrite the default initial value, one can define a wrapper around this function: 

```Julia
using ComponentArrays

custom_ind_statevars(p; kwargs...) = ComponentVector(
    EcotoxSystems.initialize_individual_statevars(); 
    new_statevar = 0.
    )
```

For IBM simulations, the new function `custom_ind_statevars` can be supplied as keyword-argument 
`init_individual_statevars` to `IBM_simulator`:

```Julia
sim = IBM_simulator(p; init_individual_statevars = custom_ind_statevars)
```

For the ODE, simulator, the state variables for all components are initialized simultaneously, 
and this function does not take any additional keyword-arguments: 

```Julia
custom_statevars(p) = ComponentVector( # state variables for some completely new model
    glb = initialize_global_statevars(p),
    spc = custom_ind_statevars(p)
)
```

Alternatively, the custom `ComponentVector` can of course be defined *from scratch* in the new initialization functions for state variables: 

```Julia
custom_statevars(p) = ComponentVector( # state variables for some completely new model
    glb = ComponentVector(t_max = 10., C_W = 0.),
    spc = ComponentVector(N = 1.)
)
```
"""
function initialize_individual_statevars(
    p::ComponentVector; 
    id = 1., 
    cohort = 0.)::ComponentVector
    
    ComponentVector(
        embryo = 1.,
        juvenile = 0.,
        adult = 0.,

        X_emb = p.spc.X_emb_int, # initial mass of vitellus
        S = p.spc.X_emb_int * X_EMB_INT_REL, # initial structure is a small fraction of initial reserve // mass of vitellus
        H = 0., # maturity
        R = 0., # reproduction buffer
        I = 0., # cumulative ingestion
        A = 0., # cumulative assimilation
        M = 0., # cumulative somatic maintenance
        J = 0., # cumulative maturity maintenance 
        f_X = 1., # scaled functional response 
        I_emb = 0., # cumulative ingestion from vitellus
        I_p = 0., # cumulative ingestion from external food resource

        D_z = constrmmat(p.spc.KD), # sublethal damage per stressor and PMoA
        D_h = constrmvec(p.spc.KD_h), # lethal damage per stressor

        y_T = 1.,

        y_z = constrmmat(p.spc.KD), # relative response per stressor and pmoa
        #y_j = constrmmat(p.spc.KD, 2, fillval = 1), # relative response per pmoa
        y_j = [0. 0. 0. 0.],
        h_z = 0., # hazard rate caused by chemical stressors
        S_z = 1., # chemical-related survival probability

        # these are curently only needed in the IBM version, 
        # but may find application in the pure-ODE implementation 

        S_max_hist = p.spc.X_emb_int * X_EMB_INT_REL, # initial reference structure
        id = id, 
        cohort = cohort,
        age = 0.,
        cause_of_death = 0.,
        time_since_last_repro = 0.,
        cum_repro = 0.,
    )
end


"""
    initialize_global_statevars(p::ComponentVector)

Function to initialize global state variables. 
In the default model, these are the resource abundance `X`, 
external chemical stressor concentration `C_W` and population size `N`. 

Global state variables can be extended, modified or replaced in the same way as individual-level state variables. 
"""
function initialize_global_statevars(p::ComponentVector)::ComponentVector
    ComponentArray( # initial states
        X = p.glb.dX_in, # initial resource abundance equal to influx rate
        C_W = p.glb.C_W, # external stressor concentrations
        N = p.glb.N0
    )
end

"""
    initialize_statevars(p::ComponentVector)::ComponentVector

For initialization of ODE simulator, initialize the component vector of state variables, `u`, based on common oaraeter collection `p`.
"""
function initialize_statevars(p::ComponentVector)::ComponentVector
    return ComponentVector(     
        glb = initialize_global_statevars(p),
        ind = initialize_individual_statevars(p)
    )
end