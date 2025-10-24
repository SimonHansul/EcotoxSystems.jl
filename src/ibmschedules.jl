# ibmschedules.jl
# generic scheduling for individuals

get_recorded_individual_var_indices(m::AbstractIBM) = map(x -> x in m.recorded_individual_vars,  keys(ind)) |> BitVector

"""
    get_global_statevars!(a::AbstractDEBIndividual, m::AbstractIBM)::Nothing

Retrieve global state variables and derivatives for use in the individual-level model step.
"""
function get_global_statevars!(a::AbstractDEBIndividual, m::AbstractIBM)::Nothing
    
    a.du.glb = m.du.glb
    a.u.glb = m.u.glb

    return nothing
end

"""
    set_global_statevars!(m::AbstractIBM, a::AbstractDEBIndividual)::Nothing

Update global state variables and derivatives based on individual-level model step.
"""
function set_global_statevars!(m::AbstractIBM, a::AbstractDEBIndividual)::Nothing

    m.du.glb = a.du.glb 
    m.u.glb = a.u.glb

    return nothing
end

"""
    individual_step!(a::Agent, m::Model)::Nothing

The individual-level model step follows a generic pattern:

First the ODE-portion of the model is executed, and the corresponding state variables are updated using the Euler scheme. 

Then the rule-based portion of the model is executed. These are all the functions which cannot / should not be expressed as part of an ODE.
At the minimum, this will include life stage transitions, reproduction and death of individuals. 

For a spatially explicit model, movement should also most likely be part of the rule-based portion, 
as well as functions which require direct information exchange between individuals.
"""
function individual_step!(a::AbstractDEBIndividual, m::AbstractIBM)
    
    get_global_statevars!(a, m) # update reference to global states
    
    a.individual_ode!(a.du, a.u, a.p, m.t) # calculate derivatives of the ODE-portion of the individual model
    Euler!(a.u, a.du, m.dt) # update states using Euler scheme -> this could even be replaced with a more sophisticated scheme
    a.individual_rules!(a, m) # apply the rule-based portion of the individual model

    set_global_statevars!(m, a) # update global states

    return nothing
end

"""
    Euler!(u::ComponentVector, du::ComponentVector, dt::Real)::Nothing

Apply Euler scheme to state variables.
"""
function Euler!(u::ComponentVector, du::ComponentVector, dt::Real)::Nothing
    u .+= du .* dt
    return nothing
end


"""
    record_individual!(a::AbstractDEBIndividual, m::AbstractIBM)::Nothing

Store individual-level state variables in `m.individual_record`.
"""
function record_individual!(a::AbstractDEBIndividual, m::AbstractIBM)::Nothing

    if m.record_individuals && isapprox(m.t % m.saveat, 0, atol = m.dt)
        push!(
            m.individual_record,
            ComponentVector(a.u.ind; t = m.t)
        )
    end

    return nothing
end


"""
    record_global!(m::AbstractIBM)::Nothing

Store global state variables in `m.global_record`.
"""
function record_global!(m::AbstractIBM)::Nothing

    if isapprox(m.t % m.saveat, 0, atol = m.dt)
        push!(m.global_record,ComponentVector(m.u; t = m.t))
    end

    return nothing
end


"""
    filter_individuals!(m::AbstractIBM)

Remove individuals which have been flagged to die after the current time-step. 

Individuals for which the condition `u.spc.cause_of_death == 0` applies are retained.
"""
filter_individuals!(m::AbstractIBM) = m.individuals = filter(x -> x.u.spc.cause_of_death == 0, m.individuals)

function step_all_individuals!(m::AbstractIBM)::Nothing
    
    shuffle!(m.individuals)

    for a in m.individuals
        # before an individual step is executed, the global derivatives are reset to 0. this
        # this is so that individuals can modify the global states, e.g. when ingesting food
        m.du.glb .= 0. 
        individual_step!(a, m)
        record_individual!(a, m)
    end

    filter_individuals!(m)

    return nothing
end


"""
    default_global_rules!(m)

Global rule-based portion of the default model. 
"""
function default_global_rules!(m)
    m.u.glb.N = length(m.individuals) # tracking population size
    m.u.glb.X = max.(0, m.u.glb.X) # HOTFIX : negative resource abundances can cause chaos
end


"""
    model_step!(m::AbstractIBM)::Nothing

Generic definition of an individual-based model step.
"""
function model_step!(m::AbstractIBM)::Nothing
    # calculate global derivatives
    # change in resource abundance, chemical stressor exposure etc.
    
    m.global_ode!(m.du, m.u, m.p, m.t)
    Euler!(m.u, m.du, m.dt) 
    m.global_rules!(m)
    step_all_individuals!(m)
    
    # global statevars are updated after individual derivatives are calculated
    # this is important because individuals affect global states using mutating operators
    m.t += m.dt

    record_global!(m)

    return nothing
end