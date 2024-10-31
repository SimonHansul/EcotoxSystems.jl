using Pkg; Pkg.activate("test")
using Test
using Distributions

using OrdinaryDiffEq
using DataFrames, ComponentArrays
using Plots, StatsPlots
using StatsBase


default(leg = false)

using Revise
@time using EcotoxSystems

# creating a sysimage for the tests
#using PackageCompiler
#create_sysimage(
#    ["Plots", "StatsPlots", "Revise", "OrdinaryDiffEq", "StatsBase", "BenchmarkTools", "Chain", "DataFrames", "DataFramesMeta", "Distributions", "Test"],
#    sysimage_path = projectdir("TestSysimage.so")
#    )


#TODO: include randomized inputs in each test

include("test01_defaults.jl") # simulates the default parameters
include("test02_debparamsweeps.jl") # performs a parameter sweep
include("test03_singlestressors.jl") # simulates single stressors
include("test04_tempcorr.jl") # simulates temperature effects
include("test05_mixtures.jl") # simulates chemical mixtures





