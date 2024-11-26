using Pkg; Pkg.activate("test")

import EcotoxSystems: defaultparams, ODE_simulator
using 
sim = ODE_simulator(defaultparams)



using Test
using Distributions

using OrdinaryDiffEq
using DataFrames, ComponentArrays
using Plots, StatsPlots
default(leg = false)
using StatsBase

using Revise
@time using EcotoxSystems

#TODO: include randomized inputs in each test

include("test01_defaults.jl") # simulates the default parameters
include("test02_food.jl") # performs a parameter sweep
include("test03_toxicity.jl") # simulates single stressors
include("test04_temperature.jl") # simulates temperature effects
include("test05_mixtures.jl") # simulates chemical mixtures

include("Aqua.jl")
