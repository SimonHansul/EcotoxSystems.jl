# default_statevars.jl 
# functions to generate default state variables as ComponentArrays
# state variables are differentiated into global and individual-level
# this differentiation is essential for the functionality of the generic IBM schedule, 
# and to allow for extensions with mulitple types of agents

const X_EMB_INT_REL = 0.001 # for the default DEB model, this global constant determines the initial structural mass, relative to the mass of an egg


function initialize_individual_statevars(
    p::ComponentVector; 
    id = 1., 
    cohort = 0.)::ComponentVector
    
    ComponentVector(
        is_embryo = 1.,
        is_juvenile = 0.,
        is_adult = 0.,

        X_emb = p.ind.X_emb_int, # initial mass of vitellus
        S = p.ind.X_emb_int * X_EMB_INT_REL, # initial structure is a small fraction of initial reserve // mass of vitellus
        H = 0., # maturity
        R = 0., # reproduction buffer
        I = 0., # cumulative ingestion
        A = 0., # cumulative assimilation
        M = 0., # cumulative somatic maintenance
        J = 0., # cumulative maturity maintenance 
        f_X = 1., # scaled functional response 
        I_emb = 0., # cumulative ingestion from vitellus
        I_p = 0., # cumulative ingestion from external food resource

        S_max_hist = p.ind.X_emb_int * X_EMB_INT_REL, # initial reference structure
        id = id, 
        cohort = cohort,
        age = 0.,
        cause_of_death = 0.,
        time_since_last_repro = 0.,
        cum_repro = 0.
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