module EcotoxSystems

using Parameters
using ComponentArrays
using OrdinaryDiffEq
using Distributions
using DataFrames
using PrecompileTools
using StaticArrays
using StatsBase
using Random
using Base.Threads

using RecipesBase

abstract type AbstractIBM end
abstract type AbstractParams end
const CVOrParamStruct = Union{ComponentVector,AbstractParams} # make it possible to use component vectors or custom structs to store parameters

include("utils.jl")
include("drcfuncts.jl")

# infrastructure to incorporate ODE-based model into an IBM framework

include("individuals.jl")
include("individualbasedmodel.jl")
include("ibmschedules.jl")


# pre-defined model(s)

include("models/debkiss/params.jl")
include("models/debkiss/derivatives.jl")
include("models/debkiss/statevars.jl")
include("models/debkiss/rules.jl")
include("models/debkiss/traits.jl")
include("models/debkiss/debkiss.jl")

export SimplifiedEnergyBudget, instantiate, simulate_ode, simulate_ibm, simulate


# functions to simulate models

include("simulators.jl")
export @replicates, replicates, treplicates, exposure

include("plotrecipes.jl")
export lineplot, groupedlineplot, rugplot

end # module EcotoxSystems
