module EnergyBudgetDiffEqs

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

include("utils.jl")
export relative_response

include("drcfuncts.jl")
include("paramstructs.jl")
include("derivatives.jl")
include("statevars.jl")

# infrastructure to incorporate ODE-based model into an IBM framework

abstract type AbstractDEBIBM end

include("individuals.jl")
include("individualbasedmodel.jl")
include("ibmschedules.jl")

include("simulators.jl")
export ODE_simulator, IBM_simulator, @replicates, replicates, treplicates, exposure

include("traits.jl")

end # module EnergyBudgetDiffEqs
