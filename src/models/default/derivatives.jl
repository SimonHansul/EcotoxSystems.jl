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

#condition_juvenile(u, t, integrator) = u.spc.X_emb # transition to juvenile when X_emb hits 0
#function effect_juvenile!(integrator) 
#    integrator.u.spc.embryo = 0.
#    integrator.u.spc.juvenile = 1.
#    integrator.u.spc.adult = 0.
#end
# 
#condition_adult(u, t, integrator) = integrator.p.spc.H_p - u.spc.H # condition to adult when H reaches H_p
#function effect_adult!(integrator) 
#    integrator.u.spc.embryo = 0.
#    integrator.u.spc.juvenile = 0.
#    integrator.u.spc.adult = 1.
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

@inline function determine_life_stage!(du, u, p, t)::Nothing 
    
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

    @unpack glb, spc = u

    # scaled damage dynamics based on the minimal model
    
    spc.y_j .= 1.0 # reset relative responses 
    spc.h_z = 0.0 # reset GUTS-SD hazard rate

    for z in eachindex(glb.C_W) # for every chemical
        for j in eachindex(spc.y_j) # for every PMoA
            # calculate change in damage
            du.spc.D_z[z,j] = minimal_TK(spc.embryo, p.spc.KD[z,j], glb.C_W[z], spc.D_z[z,j]) #(1 - spc.embryo) * p.spc.KD[z, j] * (glb.C_W[z] - spc.D_z[z, j])
            # update relative response with respect to PMoA j
            if j != 2 # PMoAs with decreasing response
                spc.y_j[j] *= LL2(spc.D_z[z,j], p.spc.E[z,j], p.spc.B[z,j])
            else # PMoAs with increasing response
                spc.y_j[j] *= LL2pos(spc.D_z[z,j], p.spc.E[z,j], p.spc.B[z,j])
            end
        end
        # calculate change in damage for lethal effects
        du.spc.D_h[z] = (1 - spc.embryo) * p.spc.KD_h[z] * (glb.C_W[z] - spc.D_h[z])
        # update hazard rate
        spc.h_z += LL2GUTS(spc.D_h[z], p.spc.E_h[z], p.spc.B_h[z])
    end

    du.spc.S_z = -spc.h_z * spc.S_z # survival probability according to GUTS-RED-SD
    
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

    @unpack glb, spc = u

    determine_life_stage!(du.spc, spc, p.spc, t)
    spc.y_T = y_T(p.spc.T_A, p.spc.T_ref, p.glb.T) # temperature correction
    spc.f_X = f_X(glb.X, p.glb.V_patch, p.spc.K_X) # ingestion rates and feedback with resource pools    
    dI_emb = dI_embryo(spc.embryo, spc.S, p.spc.dI_max_emb, spc.y_T) # resource uptake for embryos
    dI_all = dI(spc.embryo, spc.f_X, p.spc.dI_max, spc.S, spc.y_T) # resource uptake after birth    
    du.spc.I = dI_emb + dI_all # current resource uptake
    du.glb.X -= dI_all  # change in external resource abundance
    du.spc.X_emb = -dI_emb  # change in vitellus (yolk)
    du.spc.A = dA(du.spc.I, p.spc.eta_IA, spc.y_j[3]) # assimilation rate
    du.spc.M = dM(spc.S, p.spc.k_M, spc.y_j[2], spc.y_T) # somatic maintenance rate
    du.spc.J = dJ(spc.H, p.spc.k_J, spc.y_j[2], spc.y_T) # maturity maintenance rate
    du.spc.S = dS(p.spc.kappa, du.spc.A, du.spc.M, p.spc.eta_SA, spc.y_j[1], p.spc.eta_AS)
    du.spc.H = dH(spc.adult, p.spc.kappa, du.spc.A, du.spc.J) # maturation rate
    du.spc.R = dR(spc.adult, p.spc.eta_AR, spc.y_j[4], p.spc.kappa, du.spc.A, du.spc.J) # reproduction rate

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