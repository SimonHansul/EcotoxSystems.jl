# ======================================== #
# Achieving allocation-free access
# to deeply nested properties
# ======================================== #

#=
Problem: wrapping component vectors in mutable structs creates type instability

the code below suggests that we can simply convert our `Individual` 
to a `NamedTuple`, and by doing so remove the allocations

we could keep the `Individual(...)` function to genereate an individual, 
    in order to provide a familiar API. 

once we are working with multi-species systems, we will have to think about
    how to replicate multiple dispatch behaviour using this approach - 
    but we can keep that for later
=#

abstract type AbstractAgent end

abstract type DaphniaMagna <: AbstractAgent end
abstract type ScapholeberisMucronata <: AbstractAgent end

using Pkg; Pkg.activate("test")
Pkg.rm("EcotoxSystems")
Pkg.develop(path = ".")
Pkg.precompile()
using Distributions
using DataFramesMeta
using DataFrames
using Test
using Plots, StatsPlots

using Revise
using EcotoxSystems, EcotoxSystems.DEBkiss, EcotoxSystems.IBM

# re-writing the IBM simulator to avoid function wrapping altogether...

global_ode! = DEBkiss.constant_nutrient_influx!
global_rules! = DEBkiss.default_global_rules!

DEBkiss.debkiss!
individual_ode! = DEBkiss.debkiss!
individual_rules! = DEBkiss.individual_rules!

using EcotoxSystems.ComponentArrays, EcotoxSystems.Parameters

p = DEBkiss.defaultparams()
u_glb = DEBkiss.initialize_global_statevars(p)

@time a = IBM.Individual(
    p, 
    u_glb, 
    DEBkiss.initialize_individual_statevars,
    DEBkiss.generate_individual_params
);


modify!(du) = du.ind.S += 1e-3

function foo(a::IBM.Individual)
    for _ in 1:1000
        modify!(a.du)
    end
end

using BenchmarkTools
 # with "barrier" : 16000; without: 125392
allocs = @allocated foo(a)

function foo2(du::ComponentVector, u::ComponentVector)
    for _ in 1:1000
        modify!(du)
    end
end

@allocated foo2(a.du, a.u)

function foo3(a::NamedTuple)
    for _ in 1:1000
        modify!(a.du)
    end
end
atup = (du = a.du, u = a.u)
@allocated foo3(atup)


# ========================================================== #
# Attempting allocation-free access in multi-species systems
# ========================================================== #

#=
Option 1: dummy type for species, only used for dispatch. 
Access allocates a little memory in hot-loop (48); 
might be a test artefact rather than worrysome allocation...TBD
=#

abstract type Species end
struct DaphniaMagna <: Species end
struct DaphniaLongispina <: Species end

a1 = (
    species = DaphniaMagna(), 
    u = ComponentVector(a = 1, b = 2., c = 3.)
)

a2 = (
    species = DaphniaLongispina(), 
    u = ComponentVector(a = 1., b = 2., c = 3.)
)

function fun!(::DaphniaMagna, u::ComponentVector)
    u.a += 1e-3 
    return nothing
end

function fun!(::DaphniaLongispina, u::ComponentVector)
    u.b -= 1e-3
    return nothing
end

function fun_allocs(species, u)
    for _ in 1:1000
        fun!(species, u)
    end
end

@allocated fun_allocs(a1.species, a1.u)
@allocated fun_allocs(a2.species, a1.u)

# would be nice to have sonething like a @species macro that defines all the functions we need
# (no priority)
#@species :DaphniaMagna begin
    
    # struct DaphniaMagna <: Species end 

    #initialize_individual_statevars = ...
    #generate_individual_params = ...
    #individual_ode! = ...
    #individual_rules! = ...

    # makestep(DaphniaMagna)...
#end