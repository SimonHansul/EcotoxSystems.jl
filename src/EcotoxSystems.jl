module EcotoxSystems

using Parameters
using ComponentArrays
using OrdinaryDiffEq
using Distributions
using DataFrames
using PrecompileTools
using StatsBase
using Random
using Base.Threads

using RecipesBase

abstract type AbstractIBM end
abstract type AbstractParams end
const ComponentVector = Union{ComponentVector,AbstractParams} # make it possible to use component vectors or custom structs to store parameters

include("utils.jl")
include("drcfuncts.jl")

# functions to simulate models

include("simulators.jl")
export @replicates, replicates, treplicates, exposure

include("plotrecipes.jl")
export lineplot, groupedlineplot, rugplot

abstract type AbstractIBM end

module IBM

import ..AbstractIBM, ..getcolnames
using ComponentArrays, Parameters
using DataFrames
using Random
using ProgressMeter


function Individual(
    p::ComponentVector, 
    u_glb::ComponentVector, 
    initialize_individual_statevars,
    generate_individual_params; 
    cohort::Real = 1., id::Real = 1.
    )

    p_ind = generate_individual_params(p)
    u = ComponentVector(
        glb = u_glb,
        ind = ComponentVector(
            initialize_individual_statevars(p_ind),
            age = 0.,
            id = id, 
            cohort = cohort)
        )
  
    du = similar(u)
    du .= 0.

    return (
        u = u, 
        du = du,
        p = p_ind
    )
end


function Model(
    p::ComponentVector;
    init_u_glb,
    init_u_ind,
    gen_p_ind,
    dt::Real,
    N0::Real,
    saveat::Real = 1,
    record_individuals::Bool = false
    )

    u_glb = init_u_glb(p)
    u = ComponentVector(
        glb = ComponentVector(u_glb)
        )
    du = similar(u)
    du .= 0   

    global_record = ComponentVector[]
    individual_record = ComponentVector[]
    individuals = instantiate_population(N0, p, u_glb, init_u_ind, gen_p_ind)

    aux = ComponentVector(
        idcount = N0,
        N = N0,
        t = 0., 
        dt = dt, 
        saveat = saveat, 
        record_individuals = record_individuals
    )
   
    m =  (
        du = du, 
        u = u,
        p = p,
        individuals = individuals,
        aux = aux,
        global_record = global_record,
        individual_record = individual_record
        )
    
   

   return m
end

function instantiate_population(
    N0,
    p, 
    u_glb, 
    init_u_ind, 
    gen_p_ind
    )

    a = IBM.Individual(
            p, 
            u_glb, 
            init_u_ind,
            gen_p_ind
    )

    IndividualType = typeof(a)
    individuals = IndividualType[]
    
    for i in 1:N0
        a = Individual(
                p, 
                u_glb, 
                init_u_ind,
                gen_p_ind
        )

        a.u.ind.id = i

        push!(individuals, a)
    end

    return individuals
end

function push_individual!(m::NamedTuple, a::NamedTuple)::Nothing
    
    m.aux.idcount += 1
    a.u.ind.id = m.aux.idcount
    push!(m.individuals, a)
    
    return nothing
end

function _get_global_statevars!(a, m)::Nothing
    
    a.du.glb = m.du.glb
    a.u.glb = m.u.glb

    return nothing
end

function _Euler!(u::ComponentVector, du::ComponentVector, dt::Real)::Nothing

    u .+= du .* dt
    
    return nothing
end

function _set_global_statevars!(m, a)::Nothing

    m.du.glb = a.du.glb 
    m.u.glb = a.u.glb

    return nothing
end

function _filter_individuals!(m)::Nothing
    filter!(x -> x.u.ind.aux.cause_of_death == 0., m.individuals)
    return nothing
end

function _record_individual!(a, m)::Nothing

    if (m.aux.record_individuals > 0) && isapprox(m.aux.t % m.aux.saveat, 0, atol = m.aux.dt)
        push!(
            m.individual_record,
            ComponentVector(a.u.ind; t = m.aux.t)
        )
    end

    return nothing
end

function _record_global!(m)::Nothing

    if isapprox(m.aux.t % m.aux.saveat, 0, atol = m.aux.dt)
        push!(m.global_record,ComponentVector(m.u; t = m.aux.t, N = m.aux.N))
    end

    return nothing
end



"""
Second-order function to capture individual step.
"""
function make_individual_step(ind_ode!, ind_rules!, init_u_ind, gen_p_ind)

    function individual_step!(a, m)::Nothing

        _get_global_statevars!(a, m)
        ind_ode!(a.du, a.u, a.p, m.aux.t)
        _Euler!(a.u, a.du, m.aux.dt)
        ind_rules!(a, m; init_u_ind = init_u_ind, gen_p_ind = gen_p_ind)
        _set_global_statevars!(m, a)
        return nothing
    end

    return individual_step!
end

function _step_all_individuals!(m, individual_step!)::Nothing
    
    shuffle!(m.individuals)

    for a in m.individuals
        # before an individual step is executed, the global derivatives are reset to 0. this
        # this is so that individuals can modify the global states, e.g. when ingesting food
        m.du.glb .= 0. 
        individual_step!(a, m)
        _record_individual!(a, m)
    end

    _filter_individuals!(m)

    return nothing
end


function make_model_step(
    global_ode!, 
    global_rules!, 
    individual_step!
    )

    function model_step!(m)::Nothing
        global_ode!(m.du, m.u, m.p, m.aux.t)
        _Euler!(m.u, m.du, m.aux.dt)
        global_rules!(m)
        _step_all_individuals!(m, individual_step!)
        m.aux.t += m.aux.dt
        _record_global!(m)
        return nothing
    end

    return model_step!

end


global_record_to_df(m)::DataFrame = DataFrame(hcat(m.global_record...)', getcolnames(m))

function individual_record_to_df(
    m; 
    )::DataFrame

    cols = vcat(
        [:t, :id],
        [keys(m.individual_record[1])...]
    )  |> unique

    
    hcat([map(x -> getproperty(x, y), m.individual_record) for y in cols]...) |> 
    x -> DataFrame(x, cols)
end

function simulate(
    p::ComponentVector;
    global_ode!,
    global_rules!, 
    individual_ode!,
    individual_rules!,
    init_u_glb, 
    init_u_ind, 
    gen_p_ind,
    dt = 1/24, 
    N0 = 10,
    record_individuals = false,
    saveat = 1/24
    )

    m = Model(
        deepcopy(p); 
        init_u_glb = init_u_glb,
        init_u_ind = init_u_ind,
        gen_p_ind = gen_p_ind,
        dt = dt,
        N0 = N0, 
        record_individuals = record_individuals, 
        saveat = saveat
    )

    individual_step! = IBM.make_individual_step(
        individual_ode!, 
        individual_rules!, 
        init_u_ind,
        gen_p_ind
    )

    model_step! = IBM.make_model_step(
        global_ode!, 
        global_rules!, 
        individual_step!
        )

    while !(m.aux.t > m.p.glb.t_max)
        model_step!(m)
    end

    df_spc =  record_individuals ?  IBM.individual_record_to_df(m) : DataFrame()

    return (glb = IBM.global_record_to_df(m), spc = df_spc)
end

end


module DEBkiss

using Parameters, ComponentArrays, OrdinaryDiffEq, Distributions
using DataFrames, StatsBase, Random

import ..sol_to_df
import ..AbstractIBM
import ..IBM: Individual

include("models/debkiss/params.jl")
include("models/debkiss/callbacks.jl")
include("models/debkiss/derivatives.jl")
include("models/debkiss/statevars.jl")
include("models/debkiss/traits.jl")
include("models/debkiss/debkiss.jl")
include("models/debkiss/rulebased.jl")


#import ..IBM
#function sim_ibm(p)
#
#    return IBM.simulate(
#        p; 
#        global_ode! = constant_nutrient_influx!,
#        individual_ode! = debkiss!,
#
#        global_rules! = default_global_rules!,
#        individual_rules! = individual_rules!,
#
#        init_u_glb = initialize_global_statevars,
#        init_u_ind = initialize_individual_statevars,
#        gen_p_ind = generate_individual_params,
#    )
#
#end


export FullDEBkiss, instantiate, simulate
end

end # module EcotoxSystems
