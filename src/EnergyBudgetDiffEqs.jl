module EnergyBudgetDiffEqs

	using Parameters
	using ComponentArrays
	using OrdinaryDiffEq
	using Distributions
	using DataFrames
	using PrecompileTools
	using StaticArrays
	using StatsBase
	using Base.Threads

    include("utils.jl")
    export relative_response
	
	include("drcfuncts.jl")

	include("paramstructs.jl")

	include("derivatives.jl")

	include("statevars.jl")

	include("simulators.jl")
	export @replicates, replicates, treplicates, exposure

	include("traits.jl")

end # module EnergyBudgetDiffEqs
