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
now for the IBM.

- the individual rules need to use the callbacks

=#

# [2026-03-06]
# 8 weeks with peak 2000 individuals = 8.6 seconds
#     - not enough ❌
#       - the goal should be clearly below 3 seconds (based on comparison with netlogo) for small populations (<10k indivuals), then linear increase in comp time
#     - inspect profview   
#       - reveals that comp time is spent mostly in individual_rules!
#       - "age" could go into derivatives
#           - minor improvement => 7.6 seconds ☑️
#      - more comp time is spent in getindex
#      - try to use @unpack instead?
#           - applied @unpack for starvation 
#           - additional minor improvement => 6 seconds ☑️
#           - appled @unpack for remaining statevars relevant for debkiss rules
#           - additional improvement => 5 seconds ☑️
#      - profiling highligts  `u.ind.time_since_last_repro += m.dt # track reproduction period`
#      - I suspect thqt the profiler can't know the type of m.dt
#      - how do we achieve that?
#           - changing the signature of debkiss_individual_rues! to use concrete types
#               - AbstractIBM => IndividualBasedModel: No improvement
#               - AbstractIndividual => Individual: No improvement
#      - it might be an issue that we have functions as fields of structs
#           - check if https://github.com/JuliaLang/FunctionWrappers.jl is useful?
#      - for now: try to write a more specialized function that skips any wrapping of functions.
#           - simulate_debkiss_population below


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
@time EcotoxSystems.simulate_IBM(debkiss, saveat = 1, showinfo = 14); # measure
#@benchmark EcotoxSystems.simulate_IBM(debkiss, saveat = 1, showinfo = 14) # measure

VSCodeServer.@profview EcotoxSystems.simulate_IBM(debkiss, saveat = 1, showinfo = 7)















#=
💡 trying an entirely differently idea

    it should be ok to define a new struct for each IBM
    maybe we can still exploit the abstract type to do some of the heavy lifting
        TODO: if this approach works, replace global_ode! with specific functions
=#

import EcotoxSystems: AbstractIndividual, get_global_statevars!, Euler!, set_global_statevars!, AbstractIBM, IndividualBasedModel

struct Phil <: AbstractIndividual

    du::ComponentVector
    u::ComponentVector
    p::ComponentVector

end

function Phil(
    p, 
    global_statevars; 
    id = 1., 
    cohort = 0.
    )

    p_ind = debkiss:generate_individual_params(p)

    # individual stores a reference to global states + copy of own states
    u = ComponentVector(
        glb = global_statevars, # global states
        ind = ComponentVector( # own states
            init_individual_statevars(p_ind, id = id, cohort = cohort);
        )
    )

    # derivatives have the same shape as statevars and start at 0
    du = similar(u)
    du .= 0.

    return Phil(
        du, 
        u,
        p_ind
    )
end

function individual_step!(a::Phil, m::AbstractIBM)::Nothing

    get_global_statevars!(a, m) # update reference to global states

    debkiss.individual_derivatives!(a.du, a.u, a.p, m.t) # calculate derivatives of the ODE-portion of the individual model
    Euler!(a.u, a.du, m.dt) # update states using Euler scheme -> this could even be replaced with a more sophisticated scheme
    debkiss.individual_rules!(a, m) # apply the rule-based portion of the individual model

    set_global_statevars!(m, a) # update global states

    return nothing
end

mutable struct PhilsModel <: AbstractIBM
    individuals::Vector{Phil}
    du::ComponentVector
    u::ComponentVector
    p::ComponentVector
    t::Float64
    dt::Float64
    idcount::Int
    saveat::Float64
    global_record::Vector{ComponentVector}
    individual_record::Vector{ComponentVector}

    record_individuals::Bool
end

function PhilsModel(
    p::ComponentVector; 

    dt = 1/24, 
    saveat = 1, 
    record_individuals::Bool = true
    )

    return PhilsModel(
        
    )


end

function IndividualBasedModel(
 
    )

    return IndividualBasedModel(
        global_ode!, 
    )

end

