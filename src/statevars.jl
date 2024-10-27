const X_EMB_INT_REL = 0.001 

function initialize_individual_statevars(p::ComponentVector, id = 1, cohort = 0)::ComponentVector
    ComponentVector(
        embryo = 1,
        juvenile = 0,
        adult = 0,

        X_emb = p.ind.X_emb_int, # initial mass of vitellus
        S = p.ind.X_emb_int * X_EMB_INT_REL, # initial structure is a small fraction of initial reserve // mass of vitellus
        S_max_hist = p.ind.X_emb_int * X_EMB_INT_REL, # initial reference structure
        H = 0, # maturity
        R = 0, # reproduction buffer
        f_X = 1, # scaled functional response 
        I_emb = 0, # ingestion from vitellus
        I_p = 0, # ingestion from external food resource
        I = 0, # total ingestion
        A = 0, # assimilation
        M = 0, # somatic maintenance
        J = 0, # maturity maintenance 

        D_G = MVector{length(p.ind.k_D_G), Float64}(zeros(length(p.ind.k_D_G))), # scaled damage | growth efficiency
        D_M = MVector{length(p.ind.k_D_M), Float64}(zeros(length(p.ind.k_D_M))), # scaled damage | maintenance costs 
        D_A = MVector{length(p.ind.k_D_A), Float64}(zeros(length(p.ind.k_D_A))), # scaled damage | assimilation efficiency
        D_R = MVector{length(p.ind.k_D_R), Float64}(zeros(length(p.ind.k_D_R))), # scaled damage | reproduction efficiency
        D_h = MVector{length(p.ind.k_D_h), Float64}(zeros(length(p.ind.k_D_h))), # scaled damage | hazard rate

        y_G = 1., # relative response | growth efficiency
        y_M = 1., # relative response | maintenance costs 
        y_A = 1., # relative response | assimilation efficiency
        y_R = 1., # relative response | reproduction efficiency
        y_T = 1., # relative response | temperature effects
        h_z = 0., # hazard rate | chemical stressors
        S_z = 1., # chemical-related survival probability

        h_fX = 0., # starvation-induced hazard rate

        # these are only needed in the IBM version
        id = id, 
        cohort = cohort,
        age = 0,
        cause_of_death = 0,
        time_since_last_repro = 0,
        cum_offspring = 0,
    )
end

function initialize_global_statevars(p::ComponentVector)
    ComponentArray( # initial states
        X_p = p.glb.Xdot_in, # initial resource abundance equal to influx rate
        C_W = p.glb.C_W # external stressor concentrations
    )
end

"""
    initialize_statevars(p::::ComponentVector, pindt::ComponentVector{Float64})::ComponentArray

For initialization of ODE simulator, initialize the component vector of state variables, `u`, based on common oaraeter collection `p`.
"""
function initialize_statevars(p::ComponentVector)::ComponentVector
    return ComponentVector(
        glb = initialize_global_statevars(p),
        ind = initialize_individual_statevars(p)
    )
end