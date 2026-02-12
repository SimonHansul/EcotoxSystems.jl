function toy_individual_rules!(a, m)::Nothing

    @unpack RperN, am, Lp = a.p.ind
    @unpack R, age, L = a.u.ind

    num_offspring = trunc(R/RperN)

    if age >= am
        u.ind.cause_of_death = 1.
    end

    if L >= Lp
        a.u.ind.is_adult = 1. 
    end

    for i in 1:num_offspring
        m.idcount += 1 # increment individual counter
            push!(m.individuals, EcotoxSystems.Individual( # create new individual and push to individuals vector
                m.p, 
                m.u.glb; 
                id = m.idcount, 
                individual_ode! = a.individual_ode!,
                individual_rules! = a.individual_rules!,
                init_individual_statevars = a.init_individual_statevars,
                gen_ind_params = a.generate_individual_params,
                )
            )
            a.u.ind.R -= a.p[:ind][:RperN] # decrease reproduction buffer
    end
    return nothing
end

function toy_global_rules!(m)::Nothing
    m.u.glb.N = length(m.individuals)
    return nothing
end