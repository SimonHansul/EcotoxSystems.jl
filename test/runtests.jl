using Pkg; Pkg.activate("test")
using Test
using Distributions

using OrdinaryDiffEq
using DataFrames, ComponentArrays
using Plots, StatsPlots
using StatsBase
default(leg = false)

rankcor

using Revise
@time import EnergyBudgetDiffEqs as DEB


include("test01_defparams.jl")
include("test02_debparamsweeps.jl")
include("test03_TD_singlestressors.jl")