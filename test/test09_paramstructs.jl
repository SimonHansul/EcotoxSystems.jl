
# test09_paramstructs.jl
# using custom mutable structs for parameters instead of component vectors
include("test00_setup.jl")

p = deepcopy(EcotoxSystems.defaultparams)

using Distributions

@testset begin 
    defaultparams = DefaultParams()
    @test defaultparams isa AbstractParamEnsemble
end;

@testset begin
    defaultparams = DefaultParams(spc = SpeciesParams(Z = Truncated(Normal(1, 0.1), 0, Inf)));
    p = EcotoxSystems.generate_individual_params(defaultparams)
    @test p.spc.Z isa Number
    @test defaultparams.spc.Z isa Distribution
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

p = DefaultParams()
p.spc.Z = Truncated(Normal(1., 0.1), 0, Inf)

EcotoxSystems.ODE_simulator(p); 
@benchmark EcotoxSystems.ODE_simulator(p)

# why is setproperty so slow
#   what if we dont allow for RealOrDist


T = typeof(spc)

new_spc = T()


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