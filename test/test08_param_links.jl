# ODE_simulator and IBM_simulator accept a param_links function, 
# which allows to provide dependencies between parameters
# this can be used outside of the simulator function, where the link is applied to the glb or spc components 
# within the simulator function we apply the link to the glb or ind component (or whatever component has been added)

using Pkg; Pkg.activate("test")

using Test
using Distributions

using OrdinaryDiffEq
using DataFrames, ComponentArrays
using Plots, StatsPlots, Plots.Measures
default(leg = false)
using StatsBase

using Revise
@time using EcotoxSystems
import EcotoxSystems: defaultparams, ODE_simulator

