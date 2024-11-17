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

abstract type AbstractDEBIBM end

include("utils.jl")

#
## definition of the default model
#
include("drcfuncts.jl")
include("default_params.jl")
include("default_derivatives.jl")
include("default_statevars.jl")

# infrastructure to incorporate ODE-based model into an IBM framework

include("individuals.jl")
include("individualbasedmodel.jl")
include("ibmschedules.jl")

# functions to simulate models

include("simulators.jl")
export ODE_simulator, IBM_simulator, @replicates, replicates, treplicates, exposure

# inferring traits from simulations - this is probably misplaced here and should go live somewhere else
# atm we sill need it

include("traits.jl")

include("plotrecipes.jl")
export lineplot, groupedlineplot, rugplot

end # module EcotoxSystems
