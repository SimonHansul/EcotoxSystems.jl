# 2026-03-17_agentsjl_integration.jl
# idea of this script: 
#   skip function wrapping and write a specialized implementaion of the debkiss model
#   however, do this without re-typing the derivatives!
#   optionally, do this via agents.jl?

using Agents
import Agents: NoSpaceAgent
using Distributions
import OrdinaryDiffEq: OrdinaryDiffEqTsit5
#using CairoMakie


@agent struct Individual(NoSpaceAgent)
    du::ComponentVector{Real}
    u::ComponentVector{Real}
    p::ComponentVector{Real}
    i#::Integrator
end

function initialize_model()




end

function model_step!(m)

    du::ComponentVector{Real}
    u::ComponentVector{Real}
    p::ComponentVector{Real}
    i#::Integrator
    m.t += m.dt

    return nothing
end

import EcotoxSystems: debkiss!

function individual_step!(a, m)

    debkiss!(a.du, a.u, a.p, m.t)
    
    p = a.p
    u = a.u

    @unpack X_emb, H, S, S_max_hist, age, time_since_last_repro, R = u.ind
    @unpack H_p, W_S_rel_crit, h_S, a_max, tau_R, X_emb_int = p.ind

    if u.ind.X_emb <= 0.
        u.ind.is_embryo = 0.
    end

    if H >= H_p
        u.ind.is_adult = 1.
    end

    # aging is implemented in a non-mechanistic manner 
    # individuals die when they exceed their maximum age a_max
    # a_max is subject to individual variability
    if death_by_aging(age, a_max)
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
    if check_reproduction_period(time_since_last_repro, tau_R) 
        # if that is the case, calculate the number of offspring, 
        # based on the reproduction buffer and the dry mass of an egg
        for _ in 1:calc_num_offspring(R, X_emb_int)
            m.idcount += 1 # increment individual counter
            push!(m.individuals, Individual( # create new individual and push to individuals vector
                m.p, 
                m.u.glb; 
                id = m.idcount, 
                cohort = Int(u.ind.cohort + 1),
                individual_ode! = a.individual_ode!,
                individual_rules! = a.individual_rules!,
                init_individual_statevars = a.init_individual_statevars,
                generate_individual_params = a.generate_individual_params,
                )
            )
            u.ind.R -= X_emb_int # decrease reproduction buffer
            u.ind.cum_repro += 1 # keep track of cumulative reproduction of the mother individual
        end
        u.ind.time_since_last_repro = 0. # reset reproduction period
    # if reproduction period has not been exceeded,
    else
        u.ind.time_since_last_repro += m.dt # track reproduction period
    end

    return agent.du.I
end


function simulate_debkiss_population(

    )



end



m = EcotoxSystems.IndividualBasedModel(
    debkiss.parameters;
    init_global_statevars = debkiss.initialize_global_statevars,
    global_ode! = debkiss.global_derivatives!
)

