using ComponentArrays

function toy_individual_statevars(p; cohort = 0., id = 1.)
    return ComponentArray(
        L = p.ind.L0, 
        R = 0., 
        is_adult = 0., 
        age = 0., 
        cohort = cohort, 
        id = id, 
        cause_of_death = 0.
    )

end

function toy_global_statevars(p)
    return ComponentArray(
        N = p.glb.N0
    )
end

function toy_statevars(p)
    return ComponentArray(
        glb = toy_global_statevars(p), 
        ind = toy_individual_statevars(p)
    )
end