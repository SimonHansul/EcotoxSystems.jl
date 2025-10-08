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
using Interpolations
#using DocStringExtensions

using RecipesBase

abstract type AbstractIBM end

include("utils.jl")

#
## definition of the default model
#

include("models/default/drcfuncts.jl")
include("models/default/params.jl")
include("models/default/derivatives.jl")
include("models/default/statevars.jl")
include("models/default/individuals.jl")
include("models/default/traits.jl")

# infrastructure to incorporate ODE-based model into an IBM framework

include("individualbasedmodel.jl")
include("ibmschedules.jl")

# functions to simulate models

include("simulators.jl")
export @replicates, replicates, treplicates, exposure

# inferring traits from simulations - this is probably misplaced here and should go live somewhere else
# atm we sill need it


include("plotrecipes.jl")
export lineplot, groupedlineplot, rugplot

end # module EcotoxSystems