#schedules.jl

get_recorded_individual_var_indices(m::AbstractDEBIBM) = map(x -> x in m.recorded_individual_vars,  keys(ind)) |> BitVector


function get_global_statevars!(a::AbstractDEBIndividual, m::AbstractDEBIBM)::Nothing
    
    a.du.glb = m.du.glb
    a.u.glb = m.u.glb

    return nothing
end

function set_global_statevars!(m::AbstractDEBIBM, a::AbstractDEBIndividual)::Nothing

    m.du.glb = a.du.glb 
    m.u.glb = a.u.glb

    return nothing
end

"""
    individual_step!(a::Agent, m::Model)::Nothing

The individual step follows a generic pattern:

First the ODE-portion of the model is executed, and the corresponding state variables are updated using the Euler scheme. 

Then the rule-based portion of the model is executed. These are all the functions which cannot / should not be expressed as part of an ODE.
At the minimum, this will include life stage transitions, reproduction and death of individuals. 

For a spatially explicit model, movement should also most likely be part of the rule-based portion, 
as well as functions which require direct information exchange between individuals.
"""
function individual_step!(a::AbstractDEBIndividual, m::AbstractDEBIBM)
    
    get_global_statevars!(a, m) # update reference to global states

    a.individual_ode!(a.du, a.u, a.p, m.t) # calculate derivatives of the ODE-portion of the individual model
    Euler!(a.u, a.du, m.dt) # update states using Euler scheme -> this could even be replaced with a more sophisticated scheme
    a.individual_rules!(a, m) # apply the rule-based portion of the individual model

    set_global_statevars!(m, a) # update global states

    return nothing
end

function Euler!(u::ComponentVector, du::ComponentVector, dt::Real)::Nothing
    u .+= du .* dt
    return nothing
end

function record_individual!(a::AbstractDEBIndividual, m::AbstractDEBIBM)::Nothing

    if m.record_individuals && isapprox(m.t % m.saveat, 0, atol = m.dt)
        push!(
            m.individual_record,
            ComponentVector(a.u.ind; t = m.t)
        )
    end

    return nothing
end

function record_global!(m::AbstractDEBIBM)::Nothing

    if isapprox(m.t % m.saveat, 0, atol = m.dt)
        push!(
            m.global_record,
            ComponentVector(m.u; t = m.t)
        )

    end

    return nothing
end

filter_individuals!(m::AbstractDEBIBM) = m.individuals = filter(x -> x.u.ind.cause_of_death == 0, m.individuals)

function step_all_individuals!(m::AbstractDEBIBM)::Nothing
    
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

function model_step!(m::AbstractDEBIBM)::Nothing
    # calculate global derivatives
    # change in resource abundance, chemical stressor exposure etc.
    
    m.global_ode!(m.du, m.u, m.p, m.t)
    Euler!(m.u, m.du, m.dt) 
    step_all_individuals!(m)
    
    # global statevars are updated after individual derivatives are calculated
    # this is important because individuals affect global states using mutating operators
    m.u.glb.X = max(0, m.u.glb.X) # HOTFIX : negative resource abundances can cause chaos
    m.t += m.dt

    record_global!(m)

    return nothing
end