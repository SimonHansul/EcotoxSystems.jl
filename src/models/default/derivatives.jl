# derivatives.jl
# derivatives of the default model(s)

"""
    sig(
        x::Real, 
        x_thr::Real,
        y_left::Real, 
        y_right::Real; 
        beta::Real = 30
        )::Real

Sigmoid switch function. 
Used to replace simple if-statements with a continuous function in ODE models. 

`y_left` and `y_right` are the function values left and right of the threshold `x_thr`.

"""
@inline function sig(
    x::Real, 
    x_thr::Real,
    y_left::Real, 
    y_right::Real; 
    beta::Real = 30.0
    )::Real

    return 1 / (1 + exp(-beta*(x - x_thr))) * (y_right - y_left) + y_left
end

"""
Clip negative values at 0 as a continuous function, using `sig`.
"""
@inline function clipneg(x::Real)::Real
    return sig(x, 0., 0., x)
end

# definition of the conditions for life stage transitions

#condition_juvenile(u, t, integrator) = u.ind.X_emb # transition to juvenile when X_emb hits 0
#function effect_juvenile!(integrator) 
#    integrator.u.ind.embryo = 0.
#    integrator.u.ind.juvenile = 1.
#    integrator.u.ind.adult = 0.
#end
# 
#condition_adult(u, t, integrator) = integrator.p.ind.H_p - u.ind.H # condition to adult when H reaches H_p
#function effect_adult!(integrator) 
#    integrator.u.ind.embryo = 0.
#    integrator.u.ind.juvenile = 0.
#    integrator.u.ind.adult = 1.
#end

@inline function is_embryo(X_emb::Real)::Real
    return sig(X_emb, 0., 0., 1.; beta = 100.)
end

@inline function is_juvenile(X_emb::Real, H_p::Real, H::Real)::Real
    return sig(X_emb, 0., 1., 0.; beta = 100.) * sig(H, H_p, 1., 0.)
end

@inline function is_adult(H_p::Real, H::Real)::Real
    return sig(H, H_p, 0., 1.)
end

@inline function determine_life_stage!(du::ComponentVector, u::ComponentVector, p::ComponentVector, t::Real)::Nothing 
    
    u.embryo = is_embryo(u.X_emb)
    u.juvenile = is_juvenile(u.X_emb, p.H_p, u.H)
    u.adult = is_adult(p.H_p, u.H)

    return nothing
end

#cb_juvenile = ContinuousCallback(condition_juvenile, effect_juvenile!, save_positions = (false,false))
#cb_adult = ContinuousCallback(condition_adult, effect_adult!, save_positions = (false,false))

DEBODE_callbacks = CallbackSet() # while life stages are defined through sigmoid functions, CB is empty 

@inline function default_global_ODE!(du, u, p, t)::Nothing
    du.glb.X = p.glb.dX_in - p.glb.k_V * u.glb.X  
    return nothing
end

@inline function minimal_TK(
    embryo::Real,
    KD::Real, 
    C_W::Real,
    D::Real
    )::Real
    return (1-embryo) * KD * (C_W - D)
end

"""
    function default_TKTD!(du, u, p, t)::Nothing

Mixture-TKTD for an arbitrary number of stressors, assuming Independent Action.
"""
@inline function default_TKTD!(du, u, p, t)::Nothing

    @unpack glb, ind = u

    # scaled damage dynamics based on the minimal model
    
    ind.y_j .= 1.0 # reset relative responses 
    ind.h_z = 0.0 # reset GUTS-SD hazard rate

    for z in eachindex(glb.C_W) # for every chemical
        for j in eachindex(ind.y_j) # for every PMoA
            # calculate change in damage
            du.ind.D_z[z,j] = minimal_TK(ind.embryo, p.ind.KD[z,j], glb.C_W[z], ind.D_z[z,j]) #(1 - ind.embryo) * p.ind.KD[z, j] * (glb.C_W[z] - ind.D_z[z, j])
            # update relative response with respect to PMoA j
            if j != 2 # PMoAs with decreasing response
                ind.y_j[j] *= LL2(ind.D_z[z,j], p.ind.E[z,j], p.ind.B[z,j])
            else # PMoAs with increasing response
                ind.y_j[j] *= LL2pos(ind.D_z[z,j], p.ind.E[z,j], p.ind.B[z,j])
            end
        end
        # calculate change in damage for lethal effects
        du.ind.D_h[z] = (1 - ind.embryo) * p.ind.KD_h[z] * (glb.C_W[z] - ind.D_h[z])
        # update hazard rate
        ind.h_z += LL2GUTS(ind.D_h[z], p.ind.E_h[z], p.ind.B_h[z])
    end

    du.ind.S_z = -ind.h_z * ind.S_z # survival probability according to GUTS-RED-SD
    
    return nothing
end

@inline function y_T(
    T_A::Real,
    T_ref::Real,
    T::Real
    )::Real

    return exp((T_A / T_ref) - (T_A / T)) 

end

@inline function f_X(
    X::Real,
    V_patch::Real,
    K_X::Real
    )::Real

    return (X / V_patch) / ((X / V_patch) + K_X)

end

"""
    dS(
        kappa::Real,
        dA::Real, 
        dM::Real, 
        eta_SA::Real,
        y_G::Real,
        eta_AS::Real
        )::Real

Structural growth rate. 

## Arguments 

- `kappa`: Allocation fraction to somatic growth and maintenance
- `dA`: Assimilation rate
- `dM`: Somatic maintenance rate 
- `eta_SA`: Growth efficiency (yield of somatic mass on assimilates)
- `y_G`: Relative response of growth efficiency
- `eta_SA`: Shrinking efficiency
"""
@inline function dS(
    kappa::Real,
    dA::Real, 
    dM::Real, 
    eta_SA::Real,
    y_G::Real,
    eta_AS::Real
    )::Real

    if kappa * dA >= dM
        return y_G * eta_AS * (kappa * dA - dM)
    else
        return -(dM / eta_SA - kappa * dA)
    end 

end

"""
    dI_embryo(
        embryo::Real,
        S::Real,
        dI_max_emb::Real,
        y_T::Real
        )::Real

Resource ingestion by embryos, i.e. uptake of yolk/vitellus. 

## Arguments

- `embryo`: State variable indicating whether current state is embryonic
- `S`: Structural mass
- `dI_max_emb`: Maxmimum size-specific ingestion rate of embryos
- `y_T`: Temperature coefficient
"""
@inline function dI_embryo(
    embryo::Real,
    S::Real,
    dI_max_emb::Real,
    y_T::Real
    )::Real

    return embryo * (Complex(S)^(2/3)).re * dI_max_emb * y_T

end

@inline function dI(
    embryo::Real,
    f_X::Real,
    dI_max::Real,
    S::Real,
    y_T::Real
    )::Real

    return (1-embryo) * f_X * dI_max * (Complex(S)^(2/3)).re * y_T

end

@inline function dR(
    adult::Real, 
    eta_AR::Real,
    y_R::Real,
    kappa::Real,
    dA::Real,
    dJ::Real
    )::Real

    return adult * clipneg(eta_AR * y_R * ((1 - kappa) * dA - dJ))  

end

@inline function dH(
    adult::Real,
    kappa::Real,
    dA::Real,
    dJ::Real
    )::Real

    return (1 - adult) * clipneg(((1 - kappa) * dA) - dJ) # maturiation flux

end

@inline function dA(
    dI::Real,
    eta_IA::Real,
    y_A::Real
    )::Real 

    return dI * eta_IA * y_A

end

@inline function dM(
    S::Real, 
    k_M::Real, 
    y_M::Real,
    y_T::Real
    )::Real

    return S * k_M * y_M * y_T

end


@inline function dJ(
    H::Real, 
    k_J::Real, 
    y_M::Real,
    y_T::Real
    )::Real

    return H * k_J * y_M * y_T

end

"""
    default_physiology!(du, u, p, t)::Nothing

Dynamics of DEBkiss model with maturity and explicit simulation of resource dynamics.
The density of structure is ignored, and instead `S^(2/3)` is applied for surface-area scaling. 
This affects the dimension of `dI_max`, but has no effect on the model dynamics.
If model output is to be compared to length data, a statistical weight-lenght relationship 
has to be applied to the model output.
"""
function default_physiology!(du, u, p, t)::Nothing

    @unpack glb, ind = u

    determine_life_stage!(du.ind, ind, p.ind, t)
    ind.y_T = y_T(p.ind.T_A, p.ind.T_ref, p.glb.T) # temperature correction
    ind.f_X = f_X(glb.X, p.glb.V_patch, p.ind.K_X) # ingestion rates and feedback with resource pools    
    dI_emb = dI_embryo(ind.embryo, ind.S, p.ind.dI_max_emb, ind.y_T) # resource uptake for embryos
    dI_all = dI(ind.embryo, ind.f_X, p.ind.dI_max, ind.S, ind.y_T) # resource uptake after birth    
    du.ind.I = dI_emb + dI_all # current resource uptake
    du.glb.X -= dI_all  # change in external resource abundance
    du.ind.X_emb = -dI_emb  # change in vitellus (yolk)
    du.ind.A = dA(du.ind.I, p.ind.eta_IA, ind.y_j[3]) # assimilation rate
    du.ind.M = dM(ind.S, p.ind.k_M, ind.y_j[2], ind.y_T) # somatic maintenance rate
    du.ind.J = dJ(ind.H, p.ind.k_J, ind.y_j[2], ind.y_T) # maturity maintenance rate
    du.ind.S = dS(p.ind.kappa, du.ind.A, du.ind.M, p.ind.eta_SA, ind.y_j[1], p.ind.eta_AS)
    du.ind.H = dH(ind.adult, p.ind.kappa, du.ind.A, du.ind.J) # maturation rate
    du.ind.R = dR(ind.adult, p.ind.eta_AR, ind.y_j[4], p.ind.kappa, du.ind.A, du.ind.J) # reproduction rate

    return nothing
end

"""
Individual-level part of the DEB-ODE model with arbitrary number of stressors, assuming IA to compute combined effects.
"""
@inline function default_individual_ODE!(du, u, p, t)::Nothing
    default_TKTD!(du, u, p, t)
    default_physiology!(du, u, p, t)
    return nothing
end

"""
DEB-ODE model with arbitrary number of stressors, assuming IA to compute combined effects. 
"""
function default_ODE!(du, u, p, t)::Nothing
    default_global_ODE!(du, u, p, t)
    default_individual_ODE!(du, u, p, t)
    return nothing
end