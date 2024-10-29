using Pkg; Pkg.activate("test")
using Test
using Distributions

using OrdinaryDiffEq
using DataFrames, ComponentArrays
using Plots, StatsPlots
using StatsBase


default(leg = false)

using Revise
@time import EnergyBudgetDiffEqs as DEB

include("test01_defparams.jl") # simulates the default parameters
include("test02_debparamsweeps.jl") # performs a parameter sweep
include("test03_singlestressors.jl") # simulates single stressors
include("test04_tempcorr.jl") # simulates temperature effects
include("test05_mixtures.jl") # simulates chemical mixtures