using OrdinaryDiffEq
using Parameters

adult_condition(u, t, integrator) = u.ind.L - integrator.p.ind.Lp
adult_affect!(integrator) = integrator.u.ind.is_adult = 1.
adult_callback = ContinuousCallback(adult_condition, adult_affect!, nothing)

# individual-level ODE
function toy_individual!(du, u, p, t)::Nothing
    
    @unpack N = u.glb
    @unpack L, is_adult = u.ind
    @unpack Lm, rB, N_e, b_e, Rm = p.ind

    # density-dependent of individual maximum length
    Lm_N = Lm * 1/(1 + (N/N_e)^b_e)

    # individual change in body length 
    du.ind.L = rB * L * (Lm_N - L)

    # reproduction rate
    du.ind.R = is_adult * Rm * (L/Lm)^2

    return nothing
end

# global ODE (none needed here)
function toy_global!(du, u, p, t)::Nothing
    return nothing
end

function toy_ODE!(du, u, p, t)::Nothing
    toy_global!(du, u, p, t)
    toy_individual!(du, u, p, t)

    return nothing
end