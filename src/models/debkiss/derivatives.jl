# derivatives.jl
# derivatives for the full DEBkiss model

birth_condition(u, t, integrator) = u.ind.X_emb
function birth_affect!(integrator)
    integrator.u.ind.is_embryo = 0.
    integrator.u.ind.is_juvenile = 1.
    integrator.u.ind.is_adult = 0.
end

birth = ContinuousCallback(
    birth_condition, 
    birth_affect!, 
    nothing # we need `neg_affect! = nothing` to tell the solver that the affect should only occur for upcrossings (condition function switches from negative to positive) 
    )

puberty_condition(u, t, integrator) = u.ind.H - integrator.p.ind.H_p
function puberty_affect!(integrator)
    integrator.u.ind.is_embryo = 0.
    integrator.u.ind.is_juvenile = 1.
    integrator.u.ind.is_adult = 1.
end
puberty = ContinuousCallback(puberty_condition, puberty_affect!, nothing)

debkiss_callback_set = CallbackSet(birth, puberty)

function isoutofdomain(u, p, t)
    return u.ind.S < 0
end

function embryo!(du, u, p, t)::Nothing

    # TODO: to make the implementation compatible with gradient-based solvers, it might be better to use NaNMath instead
    u = max.(u, 0) # this is valid for implicit solvers and using isoutofdomain()

    @unpack V_patch = p.glb
    @unpack X = u.glb
    @unpack kappa, H_p, K_X = p.ind
    @unpack X_emb, H, S = u.ind
    
    yT = f_T(p.ind.T_A, p.ind.T_ref, p.glb.T)

    # ingestion (is_embryo & post-birth)

    dI_emb = dI_embryo(S, p.ind.dI_max_emb) *  yT
    du.ind.I = dI_emb

    # resource in patch decreases by post-birth ingestion 
    du.glb.X = 0.

    # vitellus (X_emb) change: original set du.X_emb = -dI_emb
    du.ind.X_emb = -dI_emb

    # assimilation, maintenance and maturity fluxes
    du.ind.A = dA(du.ind.I, p.ind.eta_IA)      # assimilation
    du.ind.M = dM(S, p.ind.k_M) * yT  # somatic maintenance
    du.ind.J = dJ(H, p.ind.k_J) * yT  # maturity maintenance

    # structural growth (dS) and maturation / reproduction (dH, dR)
    du.ind.S = dS(kappa, du.ind.A, du.ind.M, p.ind.eta_SA, p.ind.eta_AS)
    du.ind.H = dH(kappa, du.ind.A, du.ind.J)
    du.ind.R = 0.

    return nothing 
end

function juvenile!(du, u, p, t)::Nothing

    # TODO: to make the implementation compatible with gradient-based solvers, it might be better to use NaNMath instead
    u = max.(u, 0) # this is valid for implicit solvers and using isoutofdomain()

    @unpack V_patch = p.glb
    @unpack X = u.glb
    @unpack eta_AR, kappa, H_p, K_X = p.ind
    @unpack X_emb, H, S, is_embryo, is_adult = u.ind

    yT = f_T(p.ind.T_A, p.ind.T_ref, p.glb.T)
    fX = f_X(X, V_patch, K_X)

    # ingestion (is_embryo & post-birth)

    # TODO: to make the implementation compatible with gradient-based solvers, it might be better to use NaNMath instead
    u = max.(u, 0) # this is valid for implicit solvers and using isoutofdomain()
 
    dI = dI_postbirth(fX, p.ind.dI_max, S) * yT # assuming ingestion rate to scale like a metabolic rate
    du.ind.I = dI

    du.glb.X -= dI # resource in patch decreases by post-birth ingestion
    du.ind.X_emb = 0. # vitellus (X_emb) remains constant after birth

    # remaining fluxes
    du.ind.A = dA(du.ind.I, p.ind.eta_IA)   # assimilation
    du.ind.M = dM(S, p.ind.k_M) * yT  # somatic maintenance
    du.ind.J = dJ(H, p.ind.k_J) * yT  # maturity maintenance
    du.ind.S = dS(kappa, du.ind.A, du.ind.M, p.ind.eta_SA, p.ind.eta_AS)
    du.ind.H = dH(kappa, du.ind.A, du.ind.J)
    du.ind.R = 0.

    return nothing
end

function adult!(du, u, p, t)::Nothing

    # TODO: to make the implementation compatible with gradient-based solvers, it might be better to use NaNMath instead
    u = max.(u, 0) # this is valid for implicit solvers and using isoutofdomain()
    
    @unpack V_patch = p.glb
    @unpack X = u.glb
    @unpack eta_AR, kappa, H_p, K_X = p.ind
    @unpack X_emb, H, S, is_embryo, is_adult = u.ind

    yT = f_T(p.ind.T_A, p.ind.T_ref, p.glb.T)
    fX = f_X(X, V_patch, K_X)

    # ingestion (is_embryo & post-birth)
 
    dI = dI_postbirth(fX, p.ind.dI_max, S) * yT # assuming ingestion rate to scale like a metabolic rate
    du.ind.I = dI

    du.glb.X -= dI # resource in patch decreases by post-birth ingestion
    du.ind.X_emb = 0. # vitellus (X_emb) remains constant after birth

    # remaining fluxes
    du.ind.A = dA(du.ind.I, p.ind.eta_IA)   # assimilation
    du.ind.M = dM(S, p.ind.k_M) * yT  # somatic maintenance
    du.ind.J = dJ(H, p.ind.k_J) * yT  # maturity maintenance
    du.ind.S = dS(kappa, du.ind.A, du.ind.M, p.ind.eta_SA, p.ind.eta_AS)
    du.ind.H = 0.
    du.ind.R = dR(eta_AR, kappa, du.ind.A, du.ind.J)

    return nothing
end

function sys_embryo!(du, u, p, t)
    du.glb .= 0.
    embryo!(du, u, p, t)
end

function sys_juvenile!(du, u, p, t)
    constant_nutrient_influx!(du, u, p, t)
    juvenile!(du, u, p, t)
end

function sys_adult!(du, u, p, t)
    constant_nutrient_influx!(du, u, p, t)
    adult!(du, u, p, t)
end

function sim_embryo(p::ComponentVector; saveat = [], alg = Rodas5P(), return_sol = false, kwargs...)

    p_ind = debkiss_individual_params(p)
    u0 = initialize_statevars(p_ind)
    tspan = (0,p.glb.t_max)
    prob = ODEProblem(sys_embryo!, u0, tspan, p_ind)
    sol = solve(prob, callback = birth_terminal, saveat = saveat, alg = alg, isoutofdomain = isoutofdomain)

    if return_sol
        return sol, sol.u[end], p_ind
    end

    return sol_to_df(sol), sol.u[end], p_ind
end

function sim_juvenile(p_ind::ComponentVector, u0::ComponentVector; saveat = [], alg = Rodas5P(), kwargs...)
    
    tspan = (0,p_ind.glb.t_max)
    prob = ODEProblem(sys_juvenile!, u0, tspan, p_ind, saveat = saveat, alg = alg, isoutofdomain = isoutofdomain)
    sol = solve(prob, callback = puberty_terminal)

    return sol_to_df(sol), sol.u[end], p_ind
end

function sim_adult(p_ind, u0; saveat = [], alg = Rodas5P(), kwargs...)

    tspan = (0,p_ind.glb.t_max)
    prob = ODEProblem(sys_adult!, u0, tspan, p_ind)
    sol = solve(prob, alg = alg, saveat = saveat, isoutofdomain = isoutofdomain)

    return sol_to_df(sol), sol.u[end], p_ind
end

    
"""
Simulate all life stages consecutively as separate ODE systems. 
If a life stage is not reached, the previous ones are returned.
"""
function sim_all(p::ComponentVector; kwargs...)

    sim_emb, u0juv, p_ind = sim_embryo(p; kwargs...)

    if sim_emb.t[end] >= p.glb.t_max
        return sim_emb
    end

    sim_juv, u0ad, p_ind = sim_juvenile(p_ind, u0juv; kwargs...)
    sim_juv[!,:t] = sim_juv.t .+ sim_emb.t[end]
    
    if sim_juv.t[end] >= p.glb.t_max
        return vcat(sim_emb, sim_juv)
    end

    sim_ad, uend, p_ind = sim_adult(p_ind, u0ad; kwargs...)
    sim_ad[!,:t] = sim_ad.t .+ sim_juv.t[end]
    
    return vcat(sim_emb, sim_juv, sim_ad)
end

@inline function constant_nutrient_influx!(du, u, p, t)::Nothing
    du.glb.X = p.glb.dX_in - p.glb.k_V * u.glb.X  
    return nothing
end

@inline function f_T(T_A::Real, T_ref::Real, T::Real)::Real
    return exp((T_A / T_ref) - (T_A / T))
end

@inline function f_X(X::Real, V_patch::Real, K_X::Real)::Real
    return (X / V_patch) / ((X / V_patch) + K_X)
end

@inline function dI_embryo(S::Real, dI_max_emb::Real)::Real
    return dI_max_emb * S^(2/3)
end

@inline function dI_postbirth(f::Real, dI_max::Real, S::Real)::Real
    return f * dI_max * S^(2/3)
end

@inline function dA(dI::Real, eta_IA::Real)::Real
    return dI * eta_IA
end

@inline function dM(S::Real, k_M::Real)::Real
    return S * k_M
end

@inline function dJ(H::Real, k_J::Real)::Real
    return H * k_J
end

@inline function dS(kappa::Real, dA::Real, dM::Real, eta_SA::Real, eta_AS::Real)::Real
    return Base.ifelse(
        kappa * dA >= dM, 
        eta_AS * (kappa * dA - dM),
        -(dM / eta_SA - kappa * dA)
    )
end

@inline function dR(eta_AR::Real, kappa::Real, dA::Real, dJ::Real)::Real
    return eta_AR * ((1 - kappa) * dA - dJ)
end

@inline function dH(kappa::Real, dA::Real, dJ::Real)::Real
    return ((1 - kappa) * dA) - dJ
end

@inline function minimal_TK(
    k_D::Real, 
    C_W::Real,
    D::Real
    )::Real
    return k_D * (C_W - D)
end