using Pkg; Pkg.activate("test")
using Test
using Distributions

using OrdinaryDiffEq
using DataFrames, ComponentArrays
using Plots, StatsPlots
default(leg = false)

using Revise
@time import EnergyBudgetDiffEqs as DEB



include("test01_defparams.jl")


mat = rand(2,3)
vec = rand(2,1)

mat .* vec

prod(mat, dims = 2)