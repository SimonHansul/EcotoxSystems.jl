
"""
DEBkiss model without tox component.
"""
function debkiss!(du, u, p, t)::Nothing

    @unpack C_W, V_patch = p.glb
    @unpack X = u.glb
    @unpack kappa, H_p, K_X = p.ind
    @unpack y_j, h_z = p.ind.intermediates
    @unpack X_emb, H, S, D_z, D_h = u.ind
    
    # store smoothed indicator states back into ind (same as original determine_life_stage!)
    is_embryo = switch(X_emb, 0.0; sharpness = 100.0)
    is_adult = switch(H, H_p; sharpness = 100.0)

    Spos = softclamp(u.ind.S, 1e-12) # soft-clamping to avoid odd behaviour at very small values

    # sublethal toxicity: separate damage values + for different PMoAss

    # lethal toxicity based on GUTS-RED-SD

    du.ind.D_h[1] = (1 - is_embryo) * p.ind.KD_h[1] * (C_W[1] - D_h[1]) # calculate change in damage for lethal effects
    h_z += LL2GUTS(D_h[1], p.ind.E_h[1], p.ind.B_h[1])  # update hazard rate
    du.ind.S_z = -h_z * u.ind.S_z # survival probability according to GUTS-RED-SD

    # temperature correction and feeding functional response

    yT = f_T(p.ind.T_A, p.ind.T_ref, p.glb.T)
    fX = f_X(X, V_patch, K_X)

    # ingestion (embryo & post-birth)

    dI_emb = dI_embryo(is_embryo, Spos, p.ind.dI_max_emb) *  yT
    dI_all = dI_postbirth(is_embryo, fX, p.ind.dI_max, Spos) * yT
    du.ind.I = dI_emb + dI_all

    # resource in patch decreases by post-birth ingestion (same as original)
    du.glb.X -= dI_all

    # vitellus (X_emb) change: original set du.X_emb = -dI_emb
    du.ind.X_emb = -dI_emb

    # assimilation, maintenance and maturity fluxes
    du.ind.A = dA(du.ind.I, p.ind.eta_IA) * y_j[3]     # assimilation
    du.ind.M = dM(Spos, p.ind.k_M) * y_j[2] * yT  # somatic maintenance
    du.ind.J = dJ(H, p.ind.k_J) * y_j[2] * yT  # maturity maintenance

    # structural growth (dS) and maturation / reproduction (dH, dR)
    du.ind.S = dS(kappa, du.ind.A, du.ind.M, p.ind.eta_SA, p.ind.eta_AS * y_j[1])
    du.ind.H = dH(is_adult, kappa, du.ind.A, du.ind.J)
    du.ind.R = dR(is_adult, p.ind.eta_AR * y_j[4], kappa, du.ind.A, du.ind.J)

    return nothing
end


"""
DEBkiss model with mixture toxicity, assuming independent action (IA) across all components. 
Supports an arbitrary number of stressors and stressor/PMoA combinations.
"""
function debkiss_mixture_IA!(du, u, p, t)::Nothing

    # unpacking is optional; we can unpack at will in any order using @unpack macro
    @unpack C_W, V_patch = p.glb
    @unpack X = u.glb
    @unpack kappa, H_p, K_X = p.ind
    @unpack y_j, h_z = p.ind.intermediates
    @unpack X_emb, H, S, D_z, D_h = u.ind
    
    is_embryo = switch(X_emb, 0.0; sharpness = 100.0)
    is_adult = switch(H, H_p; sharpness = 100.0)

    Spos = softclamp(u.ind.S, 1e-12) # soft-clamping to avoid odd behaviour at very small values
   
    h_z = 0.
    y_j .= 1.

    for z in eachindex(C_W) # for every chemical
        for j in eachindex(y_j) # for every PMoA
            du.ind.D_z[z,j] = minimal_TK(is_embryo, p.ind.KD[z,j], C_W[z], D_z[z,j]) # calculate change in damage
            # update relative response with respect to PMoA j
            if j != 2 # PMoAs with decreasing response
                y_j[j] *= LL2(D_z[z,j], p.ind.E[z,j], p.ind.B[z,j])
            else # PMoAs with increasing response
                y_j[j] *= LL2pos(D_z[z,j], p.ind.E[z,j], p.ind.B[z,j])
            end
        end

        # calculate change in damage for lethal effects
        du.ind.D_h[z] = (1 - is_embryo) * p.ind.KD_h[z] * (C_W[z] - D_h[z])
        
        # update hazard rate
        h_z += LL2GUTS(D_h[z], p.ind.E_h[z], p.ind.B_h[z])
    end

    du.ind.S_z = -h_z * u.ind.S_z # survival probability according to GUTS-RED-SD
    
    # store smoothed indicator states back into ind (same as original determine_life_stage!)
    embryo = is_embryo
    adult = is_adult

    fT = f_T(p.ind.T_A, p.ind.T_ref, p.glb.T) # temp corr 
    fX = f_X(X, V_patch, K_X) # funct resp

    # ingestion (embryo & post-birth)
    dI_emb = dI_embryo(embryo, Spos, p.ind.dI_max_emb) * fT
    dI_all = dI_postbirth(embryo, fX, p.ind.dI_max, Spos) * fT
    du.ind.I = dI_emb + dI_all

    # resource in patch decreases by post-birth ingestion
    du.glb.X -= dI_all

    # change in vitellus (X_emb)
    du.ind.X_emb = -dI_emb

    # assimilation, maintenance and maturity fluxes
    du.ind.A = dA(du.ind.I, p.ind.eta_IA) * y_j[3] # assimilation
    du.ind.M = dM(Spos, p.ind.k_M * y_j[2] * fT)  # somatic maintenance
    du.ind.J = dJ(H, p.ind.k_J * y_j[2] * fT)  # maturity maintenance

    # structural growth (dS) and maturation / reproduction (dH, dR)
    du.ind.S = dS(kappa, du.ind.A, du.ind.M, p.ind.eta_SA, p.ind.eta_AS * y_j[1])
    du.ind.H = dH(adult, kappa, du.ind.A, du.ind.J)
    du.ind.R = dR(adult, p.ind.eta_AR * y_j[4], kappa, du.ind.A, du.ind.J)

    return nothing
end

@inline function constant_nutrient_influx!(du, u, p, t)::Nothing
    du.glb.X = p.glb.dX_in - p.glb.k_V * u.glb.X  
    return nothing
end

@inline softclamp(x, xmin; δ=1e-12) =
    xmin + 0.5 * ((x - xmin) + sqrt((x - xmin)^2 + δ))

@inline function switch(x, thr; sharpness=100.0)
    0.5 * (1 + tanh(sharpness * (x - thr)))
end

@inline function f_T(T_A::Real, T_ref::Real, T::Real)::Real
    return exp((T_A / T_ref) - (T_A / T))
end

@inline function f_X(X::Real, V_patch::Real, K_X::Real)::Real
    return (X / V_patch) / ((X / V_patch) + K_X)
end

@inline function dI_embryo(embryo::Real, S::Real, dI_max_emb::Real)::Real
    return embryo * S^(2/3) * dI_max_emb
end

@inline function dI_postbirth(embryo::Real, f::Real, dI_max::Real, S::Real)::Real
    return (1 - embryo) * f * dI_max * S^(2/3)
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
    if kappa * dA >= dM
        return eta_AS * (kappa * dA - dM)
    else
        return -(dM / eta_SA - kappa * dA)
    end
end

@inline function dR(adult::Real, eta_AR::Real, kappa::Real, dA::Real, dJ::Real)::Real
    return adult * softclamp(eta_AR * ((1 - kappa) * dA - dJ), 1e-12)
end

@inline function dH(adult::Real, kappa::Real, dA::Real, dJ::Real)::Real
    return (1 - adult) * softclamp(((1 - kappa) * dA) - dJ, 1e-12)
end

@inline function minimal_TK(
    embryo::Real,
    KD::Real, 
    C_W::Real,
    D::Real
    )::Real
    return (1-embryo) * KD * (C_W - D)
end