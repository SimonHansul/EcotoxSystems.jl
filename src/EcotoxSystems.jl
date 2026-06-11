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

# definition of the default model (full DEBkiss)

module DEBkiss

using Parameters, ComponentArrays, OrdinaryDiffEq, Distributions
using DataFrames, StatsBase, Random

import ..sol_to_df
import ..AbstractIBM

include("models/debkiss/params.jl")
include("models/debkiss/callbacks.jl")
include("models/debkiss/derivatives.jl")
include("models/debkiss/statevars.jl")
include("models/debkiss/traits.jl")
include("models/debkiss/debkiss.jl")
include("models/debkiss/rules.jl")

export FullDEBkiss, instantiate, simulate
end

module IBM

import ..AbstractIBM
using ComponentArrays, Parameters

#mutable struct Individual
#    du::ComponentVector 
#    u::ComponentVector
#    p::ComponentVector
#    t::ComponentVector
#
#    function Individual(
#        p::ComponentVector, 
#        u_glb::ComponentVector, 
#        initialize_individual_statevars,
#        generate_individual_params; 
#        cohort::Real = 1, id::Real = 1
#        )
#    
#        ind = new()
#        ind.p = generate_individual_params(p)
#        ind.u = ComponentVector(
#            glb = u_glb,
#            ind = ComponentVector(
#                initialize_individual_statevars(ind.p),
#                id = id, 
#                cohort = cohort
#            )
#        )
#
#        ind.du = similar(ind.u)
#        ind.du .= 0.
#
#        #@assert typeof(generate_individual_params) != typeof(initialize_individual_statevars)
#
#        return ind
#    end
#end

function Individual(
    p::ComponentVector, 
    u_glb::ComponentVector, 
    initialize_individual_statevars,
    generate_individual_params; 
    cohort::Real = 1, id::Real = 1
    )

    p_ind = generate_individual_params(p)
    u = ComponentVector(
        glb = u_glb,
        ind = ComponentVector(
            initialize_individual_statevars(p_ind),
            id = id, 
            cohort = cohort)
        )
  
    du = similar(u)
    #@assert typeof(generate_individual_params) != typeof(initialize_individual_stateva
    return (
        u = u, 
        du = du,
        p = p_ind
    )
end

mutable struct Model <: AbstractIBM
    individuals::Vector{NamedTuple}
    du::ComponentVector
    u::ComponentVector
    p::ComponentVector
    t::Real
    dt::Real
    idcount::Int
    saveat::Float64
    global_record::Vector{ComponentVector}
    individual_record::Vector{ComponentVector}

    function Model(
        p::ComponentVector;
        initialize_global_statevars,
        initialize_individual_statevars,
        generate_individual_params,
        dt::Real,
        N0::Real,
        saveat::Real = 1
        )

        m = new()
        m.individuals = NamedTuple[]
        u_glb = initialize_global_statevars(p)
        m.u = ComponentVector(
            glb = u_glb
            )
        m.du = similar(m.u)
        m.du .= 0.

        m.p = p
        m.t = 0.
        m.dt = dt
        m.idcount = 0
        m.saveat = saveat
        m.global_record = ComponentVector[]
        m.individual_record = ComponentVector[]

        for _ in 1:N0
            a = Individual(
                p, 
                u_glb, 
                initialize_individual_statevars, 
                generate_individual_params
            )

            push_individual!(m, a)
        end

        return m
    end
end

function push_individual!(m::Model, a::NamedTuple)::Nothing
    m.idcount += 1
    a.u.ind.id = m.idcount
    push!(m.individuals, a)
    return nothing
end

#include("IBM/individualbasedmodel.jl")
#include("IBM/ibmschedules.jl")

end

end # module EcotoxSystems
