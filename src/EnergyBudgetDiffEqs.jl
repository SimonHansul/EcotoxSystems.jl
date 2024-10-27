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
	
	include("drcfuncts.jl")

	include("paramstructs.jl")
	export defaultparams

	include("derivatives.jl")
	export sig, clipneg

	include("statevars.jl")
	export initialize_statevars

	include("simulators.jl")
	export @replicates, replicates, treplicates, exposure

	include("traits.jl")

end # module EnergyBudgetDiffEqs
