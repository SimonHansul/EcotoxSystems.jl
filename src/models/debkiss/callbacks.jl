birth_condition(u, t, integrator) = u.ind.X_emb
function birth_affect!(integrator)
    integrator.u.ind.embryo = 0.
    integrator.u.ind.larva = 1.
    integrator.u.ind.metamorph = 0.
    integrator.u.ind.juvenile = 0.
    integrator.u.ind.adult = 0.
end

function birth_affect_terminal!(integrator)
    birth_affect!(integrator)
    terminate!(integrator)
end

birth = ContinuousCallback(
    birth_condition, 
    birth_affect!, 
    nothing, # we need `neg_affect! = nothing` to tell the solver that the affect should only occur for upcrossings (condition function switches from negative to positive) 
    save_positions = (true,true)
    )

birth_terminal = ContinuousCallback(
    birth_condition,
    birth_affect_terminal!,
    save_positions = (true,true)
)

puberty_condition(u, t, integrator) = u.ind.H - integrator.p.ind.H_p

function puberty_affect!(integrator)
    integrator.u.ind.embryo = 0.
    integrator.u.ind.larva = 0.
    integrator.u.ind.metamorph = 0.
    integrator.u.ind.juvenile = 0.
    integrator.u.ind.adult = 1.
end

function puberty_affect_terminal!(integrator)
    puberty_affect!(integrator)
    terminate!(integrator)
end

puberty = ContinuousCallback(
    puberty_condition, 
    puberty_affect!, 
    nothing,
    save_positions = (true,true)
    )
puberty_terminal = ContinuousCallback(
    puberty_condition, 
    puberty_affect_terminal!, 
    nothing,
    save_positions = (true,true)
    )
