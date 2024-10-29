
mutable struct IndividualBasedModel <: AbstractDEBIBM
    global_ode!::Function
    individuals::Vector{DEBIndividual}
    du::ComponentVector
    u::ComponentVector
    p::ComponentVector
    t::Real
    dt::Real
    idcount::Int
    saveat::Real
    individual_record::Vector{ComponentVector}

    record_individuals::Bool

    """
        IndividualBasedModel(
            p::ComponentVector; 
            global_ode! = DEBODE_global!,
            individual_ode! = DEBODE_individual!,
            individual_rules! = default_individual_rules!,
            dt = 1/24, 
            saveat = 1, 
            record_individuals::Bool = true
        )::IndividualBasedModel
    
    Initialization of the IndividualBasedModel structure.
    
    args
        - `p`: A compatible parameter collection (`ABMParamCollection` or correspondingly structured named tuple)
    
    kwargs
        - `dt`: Model time step [t]
    """
    function IndividualBasedModel(
        p::ComponentVector; 
        global_ode! = DEBODE_global!,
        individual_ode! = DEBODE_individual!,
        individual_rules! = default_individual_rules!,
        dt = 1/24, 
        saveat = 1, 
        record_individuals::Bool = true
        )::IndividualBasedModel

        m = new()
        m.global_ode! = global_ode!
        m.individuals = Vector{DEBIndividual}(undef, Int(p.glb.N0))
        m.u = ComponentVector(glb = initialize_global_statevars(p))
        
        #m.global_statevar_names = Symbol[keys(m.u)...]
        #m.global_statevar_indices = collect(eachindex(m.global_statevar_names))
        
        m.du = similar(m.u) 
        m.p = p
        m.t = 0
        m.dt = dt
        m.saveat = saveat
        m.idcount = 0

        m.record_individuals = record_individuals
        m.individual_record = ComponentVector[]

        # initialize individuals
        for i in 1:Int(p.glb.N0)
            m.idcount += 1
            m.individuals[i] = DEBIndividual(
                p, 
                m.u.glb;
                individual_ode! = individual_ode!,
                individual_rules! = individual_rules!,
                id = m.idcount
                )
        end

        #m.individual_statevar_names = Symbol[keys(initialize_individual_statevars(m.individuals[1].p))...]
        #
        #if p.glb.N0 > 0
        #    m.individual_statevar_indices = findall(x -> x in m.individual_statevar_names, keys(m.individuals[1].u))
        #else
        #    m.individual_statevar_indices = Int64[]
        #end

        return m
    end
end