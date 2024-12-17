using Pkg; Pkg.activate("test")

using Test
using Distributions

using OrdinaryDiffEq
using DataFrames, ComponentArrays
using Plots, StatsPlots
default(leg = false)
using StatsBase

using Revise
@time using EcotoxSystems
import EcotoxSystems: defaultparams, ODE_simulator

#TODO: include randomized inputs in each test

include("test01_defaults.jl") # simulates the default parameters
include("test02_food.jl") # performs a parameter sweep
include("test03_toxicity.jl") # simulates single stressors
include("test04_temperature.jl") # simulates temperature effects
include("test05_mixtures.jl") # simulates chemical mixtures
include("test06_ibm_defparams.jl")

include("Aqua.jl")
