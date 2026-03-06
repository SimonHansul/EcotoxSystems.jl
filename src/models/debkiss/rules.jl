
"""
    debkiss_individual_rules!(a::AbstractIndividual, m::AbstractIBM)::Nothing

Defines the default rule-based portion for DEBIndividuals. <br>

The event functions which are used as callbacks during ODE solving are here re-used to apply rules for life stage transitions.

A crude rule for starvation mortality is implemented, applying a constant hazard rate of a certain relative amount of structural mass is lost.

Reproduction is assumed to occur in fixed time intervals, according to `spc.tau_R`.
"""
function debkiss_individual_rules!(a::AbstractIndividual, m::AbstractIBM)::Nothing

    p = a.p
    u = a.u

    @unpack X_emb, H, S, S_max_hist = u.ind
    @unpack H_p, W_S_rel_crit, h_S = p.ind

    if u.ind.X_emb <= 0.
        u.ind.is_embryo = 0.
    end

    if H >= H_p
        u.ind.is_adult = 1.
    end

 
    # aging is implemented in a non-mechanistic manner 
    # individuals die when they exceed their maximum age a_max
    # a_max is subject to individual variability
    if death_by_aging(u.ind[:age], p[:ind][:a_max])
        u.ind.cause_of_death = 1.
    end
    
    # for starvation mortality, currently only a limit is set on the amount of mass that can be lost
    # this is basically only a sanity check, and the actual starvation rules should be assessed on a species-by-species basis
    u.ind.S_max_hist = determine_S_max_hist(S, S_max_hist)

    if death_by_loss_of_structure(S, S_max_hist, W_S_rel_crit, h_S, m.dt)
        u.ind.cause_of_death = 2.
    end

    # reproduction, assuming a constant reproduction period
    
    # reproduction only occurs if the reproduction period has been exceeded
    if check_reproduction_period(u.ind[:time_since_last_repro], p[:ind][:tau_R]) 
        # if that is the case, calculate the number of offspring, 
        # based on the reproduction buffer and the dry mass of an egg
        for _ in 1:calc_num_offspring(u.ind[:R], p[:ind][:X_emb_int])
            m.idcount += 1 # increment individual counter
            push!(m.individuals, Individual( # create new individual and push to individuals vector
                m.p, 
                m.u.glb; 
                id = m.idcount, 
                cohort = Int(u.ind.cohort + 1),
                individual_ode! = a.individual_ode!,
                individual_rules! = a.individual_rules!,
                init_individual_statevars = a.init_individual_statevars,
                gen_ind_params = a.generate_individual_params,
                )
            )
            u.ind.R -= p[:ind][:X_emb_int] # decrease reproduction buffer
            u.ind.cum_repro += 1 # keep track of cumulative reproduction of the mother individual
        end
        u.ind.time_since_last_repro = 0. # reset reproduction period
    # if reproduction period has not been exceeded,
    else
        u.ind.time_since_last_repro += m.dt # track reproduction period
    end

    return nothing
end