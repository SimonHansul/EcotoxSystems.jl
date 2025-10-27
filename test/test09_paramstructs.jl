
# test09_paramstructs.jl
# using custom mutable structs for parameters instead of component vectors


include("test00_setup.jl")

#p = deepcopy(EcotoxSystems.defaultparams)

using Distributions

@testset begin 
    defaultparams = DefaultParams()
    @test defaultparams isa AbstractParamEnsemble
end;

@testset begin
    defaultparams = DefaultParams();
    p = EcotoxSystems.generate_individual_params(defaultparams)
    @test p.spc.Z isa Number
    @test defaultparams.spc.individual_variability.Z isa Distribution
end

# benchmarking: mutable structs are 5 times faster here
#defaultparams = DefaultParams()
#@benchmark p = EcotoxSystems.generate_individual_params(defaultparams)
#
#defaultparams = deepcopy(EcotoxSystems.defaultparams)
#@benchmark p = EcotoxSystems.generate_individual_params(defaultparams)

using StatsPlots
using ComponentArrays
using BenchmarkTools

# rewriting  generate_individual_params

p = DefaultParams()

EcotoxSystems.ODE_simulator(p); 

@benchmark EcotoxSystems.ODE_simulator(p)

p.spc.KD[1] = 1.

myview = @view(p.spc.KD[1])
myview isa SubArray


@view p.spc.Z



@benchmark EcotoxSystems.ODE_simulator(p)

VSCodeServer.@profview_allocs [EcotoxSystems.ODE_simulator(p) for _ in 1:10]

p.spc.Z = 1.
p.spc.individual_variability.Z

sim = @replicates EcotoxSystems.ODE_simulator(p) 10 

@df sim plot(:t, :S, group = :replicate)

EcotoxSystems.ODE_simulator(p, alg = Tsit5())

p = DefaultParams(spc = SpeciesParams(Z = Truncated(Normal(1, 0.1), 0, Inf)))
sim = @replicates EcotoxSystems.ODE_simulator(p) 10
@df sim plot(:t, :S, group = :replicate)


@benchmark EcotoxSystems.ODE_simulator(p)

p = DefaultParams()
EcotoxSystems.ODE_simulator(p)
@benchmark EcotoxSystems.ODE_simulator(p, returntype = EcotoxSystems.odesol)


p = DefaultParams(spc = SpeciesParams(Z = Truncated(Normal(1, 0.1), 0, Inf)))
VSCodeServer.@profview begin 
    for _ in 1:100
        EcotoxSystems.ODE_simulator(p)
    end
end