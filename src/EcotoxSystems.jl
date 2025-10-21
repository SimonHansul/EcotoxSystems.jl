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
abstract type AbstractParamStruct end
const RealOrDist = Union{Real,Distribution}
const CVOrStruct = Union{ComponentVector,AbstractParamStruct} 

include("utils.jl")
include("drcfuncts.jl")

# definition of the default model

include("models/default/params.jl")
include("models/default/derivatives.jl")
include("models/default/statevars.jl")
include("models/default/traits.jl")

# infrastructure to incorporate ODE-based model into an IBM framework

include("individuals.jl")
include("individualbasedmodel.jl")
include("ibmschedules.jl")

# functions to simulate models

include("simulators.jl")
export @replicates, replicates, treplicates, exposure


include("plotrecipes.jl")
export lineplot, groupedlineplot, rugplot

end # module EcotoxSystems