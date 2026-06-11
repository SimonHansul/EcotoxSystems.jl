# ======================================== #
# Rule-based portion of IBM simulation
# ======================================== #

"""
    default_global_rules!(m)

Global rule-based portion of the default model. 
"""
function default_global_rules!(m)
    m.aux.N = length(m.individuals) # tracking population size
    m.u.glb.X = max.(0, m.u.glb.X) # HOTFIX : negative resource abundances can cause chaos
end

abstract type AbstractIndividual end

CAUSE_OF_DEATH = Dict(
    0 => "none",
    1 => "age"
)

@inline function determine_S_max_hist(
    S::Float64,
    S_max_hist::Float64
    )::Float64

    return max(S, S_max_hist)

end

@inline function death_by_loss_of_structure(
    S::Float64,
    S_max_hist::Float64,
    W_S_rel_crit::Float64,
    h_S::Float64,
    dt::Float64
    )::Bool

    return ((S/S_max_hist) < W_S_rel_crit) && (rand() > exp(-h_S * dt))

end

@inline function death_by_aging(
    age::Float64,
    a_max::Float64
    )::Bool

    return age >= a_max

end

@inline function check_reproduction_period(time_since_last_repro::Real, tau_R::Real)::Bool
    return time_since_last_repro >= tau_R 
end

@inline function calc_num_offspring(R::Float64, X_emb_int::Float64)::Int64
    return trunc(R / X_emb_int)
end

"""
    individual_rules!(a, m)::Nothing

Defines a minimal individual-level rule-based component for the DBEkiss model.
Note that this is just a the minimum needed to get some plausible population dynamics, 
not a fully-fledged population model. 

A crude rule for starvation mortality is implemented, applying a constant hazard rate of a certain relative amount of structural mass is lost.
Reproduction is assumed to occur in fixed time intervals, according to `spc.tau_R`.
"""
function individual_rules!(a, m; init_u_ind, gen_p_ind)::Nothing

    p = a.p
    u = a.u
    a.du.ind.aux_IBM .= 0.

    @unpack is_embryo, is_juvenile, is_adult, X_emb, H = u.ind
    @unpack H_p = p.ind

    # ======================================== #
    # life stage transitions
    # ======================================== #

    if X_emb <= 0
        u.ind.is_embryo = 0.
        if H >= H_p
            u.ind.is_juvenile = 0. 
            u.ind.is_adult = 1.
        else
            u.ind.is_juvenile = 1. 
            u.ind.is_adult = 0.
        end
    else
        u.ind.is_embryo = 1
    end

    # ======================================== #
    # mortality
    # ======================================== #

    @unpack a_max, tau_R = a.p.ind.aux_IBM

    # aging is implemented in a non-mechanistic manner 
    # individuals die when they exceed their maximum age a_max
    # a_max is subject to individual variability
    if death_by_aging(u.ind.age, a_max)
        u.ind.aux_IBM.cause_of_death = 1.
    end
    
    # for starvation mortality, currently only a limit is set on the amount of mass that can be lost
    # this is basically only a sanity check, and the actual starvation rules should be assessed on a species-by-species basis
    u.ind.aux_IBM.S_max_hist = determine_S_max_hist(u.ind.S, u.ind.aux_IBM.S_max_hist)

    if death_by_loss_of_structure(u.ind.S, u.ind.aux_IBM.S_max_hist, p.ind.aux_IBM.S_rel_crit, p.ind.aux_IBM.h_S, m.dt)
        u.ind.aux_IBM.cause_of_death = 2.
    end

    # ======================================== #
    # reproduction
    # ======================================== #
    
    # reproduction only occurs if the reproduction period has been exceeded
    if check_reproduction_period(u.ind.aux_IBM.time_since_last_repro, tau_R) 
        # if that is the case, calculate the number of offspring, 
        # based on the reproduction buffer and the dry mass of an egg

        if isnan(u.ind.R)
            @show u.ind.S u.ind.age u.glb.X u.ind.I
        end

        N = calc_num_offspring(u.ind.R, p.ind.X_emb_int)
        if isnan(N)
            println(u.ind.R, p.ind.X_emb_int)
        end

        for _ in 1:N
            m.aux.idcount += 1 # increment individual counter
            push!(
                m.individuals, 
                Individual(
                    m.p, 
                    a.u.glb, 
                    init_u_ind, 
                    gen_p_ind; 
                    cohort = a.u.ind.aux_IBM.cohort += 1,
                    id = m.aux.idcount
                )
            )
           
            u.ind.R -= p.ind.X_emb_int # decrease reproduction buffer
            u.ind.aux_IBM.cum_repro += 1 # keep track of cumulative reproduction of the mother individual
        end
        u.ind.aux_IBM.time_since_last_repro = 0. # reset reproduction period
    # if reproduction period has not been exceeded,
    else
        u.ind.aux_IBM.time_since_last_repro += m.dt # track reproduction period
    end

    return nothing
end

