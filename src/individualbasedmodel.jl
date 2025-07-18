# individualbasedmodels.jl
# defines a generic composite type for IBMs
# this type should be generic enough for most applications, 
# but if extensions are required, it is best to make use of the associated AbstractIBM type

mutable struct IndividualBasedModel <: AbstractIBM
    global_ode!::Function
    global_rules!::Function
    init_global_statevars::Function
    individuals::Vector{Individual}
    integrator::SciMLBase.DEIntegrator
    dt::Real
    idcount::Int
    saveat::Real
    global_record::Vector{ComponentVector}
    individual_record::Vector{ComponentVector}

    record_individuals::Bool

    """
        IndividualBasedModel(
            p::ComponentVector; 
            global_ode! = DEBODE_global!,
            global_rules! = default_global_rules!,
            init_global_statevars = initialize_global_statevars,
            individual_ode! = DEBkiss_individual!,
            individual_rules! = default_individual_rules!,
            init_individual_statevars = initialize_individual_statevars,
            dt = 1/24, 
            saveat = 1, 
            record_individuals::Bool = true
        )::IndividualBasedModel
    
    Initialization of the IndividualBasedModel structure.
    
    args
        - `p`: A compatible parameter collection (`ABMParamCollection` or correspondingly structured named tuple)
        
    kwargs: 
        - `global_ode!`: Global ODE-based portion of the model
        - `global_rules!`: Glboal rule-based portion of the model
        - `init_global_statevars`: Function that initializes global state variables as `ComponentVector`

        - `individual_ode!`: Individual-level ODE-based portion of the model 
        - `individual_rules!`: Individual-level rule-based portion of the model 
        - `init_individual_statevars`: Function that initializes individual-level state variables as `ComponentVector`
    
        - `dt`: Model time step [t]
        - `saveat`: Time interavals at which to save output
        - `record_individuals`: Whether to record states of each individual 
    """
    function IndividualBasedModel(
        p::ComponentVector; 
        init_global_statevars = initialize_global_statevars,
        global_ode! = DEBODE_global!,
        global_rules! = default_global_rules!,

        individual_ode! = DEBkiss_individual!,
        individual_rules! = default_individual_rules!,
        init_individual_statevars = initialize_individual_statevars,
        gen_ind_params = generate_individual_params,

        dt = 1/24, 
        saveat = 1, 
        record_individuals::Bool = true, 
        ode_alg = Tsit5()
        )::IndividualBasedModel

        m = new()
        m.global_ode! = global_ode!
        m.global_rules! = global_rules!
        m.individuals = Vector{Individual}(undef, Int(p.glb.N0))
        
        # setting up global integrator

        ode_problem = ODEProblem(
            global_ode!,
            ComponentVector(glb = init_global_statevars(p)), 
            (0.0, Inf),
            p
        )
        m.integrator = OrdinaryDiffEq.init(
            ode_problem,
            ode_alg;
            advance_to_tstop = true
        )
        
        m.dt = dt
        m.saveat = saveat
        m.idcount = 0

        m.record_individuals = record_individuals
        m.global_record = ComponentVector[]
        m.individual_record = ComponentVector[]

        # initialize individuals

        for i in 1:Int(p.glb.N0)
            m.idcount += 1
            m.individuals[i] = Individual(
                p, 
                m.integrator.u.glb;
                individual_ode! = individual_ode!,
                individual_rules! = individual_rules!,
                init_individual_statevars = init_individual_statevars,
                gen_ind_params = gen_ind_params,
                id = m.idcount
                )
        end

        return m
    end
end