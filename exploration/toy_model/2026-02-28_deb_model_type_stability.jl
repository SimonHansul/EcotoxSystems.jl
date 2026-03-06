using Pkg; Pkg.activate("test")
include(joinpath(pwd(), "test", "test00_setup.jl"))
include(joinpath(pwd(), "test", "test01_defaults.jl"))

using BenchmarkTools

debkiss = SimplifiedEnergyBudget() |> instantiate
debkiss.parameters.spc.Z = 1.
sim = simulate_ode(debkiss, alg = Rodas5())

@df sim plot(
    plot(:t, :S), 
    plot(:t, :H), 
    plot(:t, :R), 
    plot(:t, [:is_embryo :is_adult])
)

#VSCodeServer.@profview_allocs simulate_ode(debkiss, alg = Rodas5(), returntype = EcotoxSystems.odesol)
prob = simulate_ode(debkiss, alg = Rodas5(), returntype = EcotoxSystems.odeprob)
integrator = init(prob, Rodas5())

du = similar(integrator.u)
reset_du(du) = du .= 0.;
@code_warntype reset_du(du)

t = 0.

#=
the derivatives are type stable now ✅
we have small allocs (32), but should not hold us back for now
=#

@code_warntype debkiss.complete_derivatives!(du, integrator.u, integrator.p, t)
@allocated debkiss.complete_derivatives!(du, integrator.u, integrator.p, t)

#=
now for the IBM...

- the individual rules need to use the callbacks

=#

# [2026-03-06]
# 8 weeks with peak 2000 individuals = 8.6 seconds
#     not enough ❌
#       the goal should be clearly below 3 seconds (based on comparison with netlogo) for small populations (<10k indivuals), then linear increase in comp time
#     look at profview   
#       reveals that comp time is spent mostly in individual_rules!
#       "age" could go into derivatives
#           minor improvement => 7.6 seconds ☑️
#      more comp time is spent in getindex
#       try to use @unpack instead?
#           applied @unpack for starvation 
#           additional minor improvement => 6 seconds ☑️
#           appled @unpack for remaining statevars relevant for debkiss rules
#           additional improvement => 5 seconds ☑️
#      profiling highligts  `u.ind.time_since_last_repro += m.dt # track reproduction period`
#      I suspect thqt the profiler can't know the type of m.dt
#       how do we achieve that?
#       ⚠️ would it make sense if `Individual` is a struct, not a mutual struct?
#       changing to struct actually decreased performance?
#           we havce 90% performance time       


begin
    debkiss = SimplifiedEnergyBudget() |> instantiate
    p = debkiss.parameters

    p.glb.dX_in = 100_000 #100_000
    p.glb.k_V = 0.1
    p.glb.V_patch = 0.5
    p.glb.N0 = 10.
    p.glb.t_max = 56.

    p.spc.Z = Truncated(Normal(1, 0.1), 0, Inf)
    p.spc.tau_R = truncated(Normal(2., 0.2), 0, Inf)
end

EcotoxSystems.simulate_IBM(debkiss, saveat = 1, showinfo = 14); # precompile

using BenchmarkTools
@benchmark EcotoxSystems.simulate_IBM(debkiss, saveat = 1, showinfo = 14) # measure

@df sim_ibm.glb plot(:t, :N) 

VSCodeServer.@profview EcotoxSystems.simulate_IBM(debkiss, saveat = 1, showinfo = 14)