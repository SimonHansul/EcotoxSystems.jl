using Pkg; Pkg.activate("test")

using Test
using Distributions

using OrdinaryDiffEq
using DataFrames, ComponentArrays
using Plots, StatsPlots, Plots.Measures
default(leg = false)
using StatsBase

using Pkg; Pkg.develop(path = ".")
using Revise
@time using EcotoxSystems