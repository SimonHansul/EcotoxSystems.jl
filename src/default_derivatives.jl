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
    x::Float64, 
    x_thr::Float64,
    y_left::Float64, 
    y_right::Float64; 
    beta::Float64 = 30.0
    )::Float64

    return 1 / (1 + exp(-beta*(x - x_thr))) * (y_right - y_left) + y_left
end

"""
Clip negative values at 0 as a continuous function, using `sig`.
"""
@inline function clipneg(x::Float64)::Float64
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

@inline function is_embryo(X_emb::Float64)::Float64
    return sig(X_emb, 0., 0., 1.; beta = 100.)
end

@inline function is_juvenile(X_emb::Float64, H_p::Float64, H::Float64)::Float64
    return sig(X_emb, 0., 1., 0.; beta = 100.) * sig(H, H_p, 1., 0.)
end

@inline function is_adult(H_p::Float64, H::Float64)::Float64
    return sig(H, H_p, 0., 1.)
end

@inline function determine_life_stage!(du::ComponentVector, u::ComponentVector, p::ComponentVector, t::Float64)::Nothing 
    u.embryo = is_embryo(u.X_emb)
    u.juvenile = is_juvenile(u.X_emb, p.H_p, u.H)
    u.adult = is_adult(p.H_p, u.H)

    return nothing
end

#cb_juvenile = ContinuousCallback(condition_juvenile, effect_juvenile!, save_positions = (false,false))
#cb_adult = ContinuousCallback(condition_adult, effect_adult!, save_positions = (false,false))

DEBODE_callbacks = CallbackSet() # while life stages are defined through sigmoid functions, CB is empty 

@inline function DEBODE_global!(du, u, p, t)::Nothing
    du.glb.X = p.glb.dX_in - p.glb.k_V * u.glb.X  
    return nothing
end

@inline function minimal_TK(
    embryo::Float64,
    KD::Float64, 
    C_W::Float64,
    D::Float64
    )::Float64
    return (1-embryo) * KD * (C_W - D)
end

"""
    function TKTD_mix_IA!(du, u, p, t)::Nothing

Mixture-TKTD for an arbitrary number of stressors, assuming Independent Action.
"""
@inline function TKTD_mix_IA!(du, u, p, t)::Nothing

    @unpack glb, ind = u

    # scaled damage dynamics based on the minimal model
    
    p.ind.y_j .= 1.0 # reset relative responses 
    p.ind.h_z = 0.0 # reset GUTS-SD hazard rate

    for z in eachindex(glb.C_W) # for every chemical
        for j in eachindex(p.ind[:y_j]) # for every PMoA
            # calculate change in damage
            du.ind.D_z[z,j] = minimal_TK(ind.embryo, p.ind.KD[z,j], glb.C_W[z], ind.D_z[z,j]) #(1 - ind.embryo) * p.ind.KD[z, j] * (glb.C_W[z] - ind.D_z[z, j])
            # update relative response with respect to PMoA j
            if j != 2 # PMoAs with decreasing response
                p.ind[:y_j][j] *= LL2(ind.D_z[z,j], p.ind.E[z,j], p.ind.B[z,j])
            else # PMoAs with increasing response
                p.ind[:y_j][j] *= LL2pos(ind.D_z[z,j], p.ind.E[z,j], p.ind.B[z,j])
            end
        end
        # calculate change in damage for lethal effects
        du.ind.D_h[z] = (1 - ind.embryo) * p.ind.KD_h[z] * (glb.C_W[z] - ind.D_h[z])
        # update hazard rate
        p.ind.h_z += LL2GUTS(ind.D_h[z], p.ind.E_h[z], p.ind.B_h[z])
    end

    du.ind.S_z = -p.ind[:h_z] * ind.S_z # survival probability according to GUTS-RED-SD
    
    return nothing
end

@inline function y_T(
    T_A::Float64,
    T_ref::Float64,
    T::Float64
    )::Float64

    return exp((T_A / T_ref) - (T_A / T)) 

end

@inline function f_X(
    X::Float64,
    V_patch::Float64,
    K_X::Float64
    )::Float64

    return (X / V_patch) / ((X / V_patch) + K_X)

end

"""
    dS(
        kappa::Float64,
        dA::Float64, 
        dM::Float64, 
        eta_SA::Float64,
        y_G::Float64,
        eta_AS::Float64
        )::Float64

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
    kappa::Float64,
    dA::Float64, 
    dM::Float64, 
    eta_SA::Float64,
    y_G::Float64,
    eta_AS::Float64
    )::Float64

    if kappa * dA >= dM
        return y_G * eta_AS * (kappa * dA - dM)
    else
        return -(dM / eta_SA - kappa * dA)
    end 

end

"""
    dI_embryo(
        embryo::Float64,
        S::Float64,
        dI_max_emb::Float64,
        y_T::Float64
        )::Float64

Resource ingestion by embryos, i.e. uptake of yolk/vitellus. 

## Arguments

- `embryo`: State variable indicating whether current state is embryonic
- `S`: Structural mass
- `dI_max_emb`: Maxmimum size-specific ingestion rate of embryos
- `y_T`: Temperature coefficient
"""
@inline function dI_embryo(
    embryo::Float64,
    S::Float64,
    dI_max_emb::Float64,
    y_T::Float64
    )::Float64

    return embryo * (Complex(S)^(2/3)).re * dI_max_emb * y_T

end

@inline function dI(
    embryo::Float64,
    f_X::Float64,
    dI_max::Float64,
    S::Float64,
    y_T::Float64
    )::Float64

    return (1-embryo) * f_X * dI_max * (Complex(S)^(2/3)).re * y_T

end

@inline function dR(
    adult::Float64, 
    eta_AR::Float64,
    y_R::Float64,
    kappa::Float64,
    dA::Float64,
    dJ::Float64
    )::Float64

    return adult * clipneg(eta_AR * y_R * ((1 - kappa) * dA - dJ))  

end

@inline function dH(
    adult::Float64,
    kappa::Float64,
    dA::Float64,
    dJ::Float64
    )::Float64

    return (1 - adult) * clipneg(((1 - kappa) * dA) - dJ) # maturiation flux

end

@inline function dA(
    dI::Float64,
    eta_IA::Float64,
    y_A::Float64
    )::Float64 

    return dI * eta_IA * y_A

end

@inline function dM(
    S::Float64, 
    k_M::Float64, 
    y_M::Float64,
    y_T::Float64
    )::Float64

    return S * k_M * y_M * y_T

end


@inline function dJ(
    H::Float64, 
    k_J::Float64, 
    y_M::Float64,
    y_T::Float64
    )::Float64

    return H * k_J * y_M * y_T

end

"""
    DEBkiss!(du, u, p, t)::Nothing

Dynamics of DEBkiss model with maturity and explicit simulation of resource dynamics.

The density of structure is ignored, and instead `S^(2/3)` is applied for surface-area scaling. 
This affects the dimension of `dI_max`, but has no effect on the model dynamics.

If model output is to be compared to length data, a statistical weight-lenght relationship 
has to be applied to the model output.
"""
function DEBkiss_physiology!(du, u, p, t)::Nothing

    @unpack glb, ind = u

    determine_life_stage!(du.ind, ind, p.ind, t)

    # temperature correction
    p.ind.y_T = y_T(p.ind.T_A, p.ind.T_ref, p.glb.T) # temperature correction

    # ingestion rates and feedback with resource pools

    ind.f_X = f_X(glb.X, p.glb.V_patch, p.ind[:K_X])

    # calculation of resource uptake for embryos vs hatched individuals
    
    dI_emb = dI_embryo(ind.embryo, ind.S, p.ind.dI_max_emb, p.ind[:y_T])
    dI_all_but_emb = dI(ind.embryo, ind.f_X, p.ind.dI_max, ind.S, p.ind[:y_T])

    # ingestion rate is the sum of both (dI_emb and dI are mutually exclusive)
    
    du.ind.I = dI_emb + dI_all_but_emb

    du.glb.X -= dI_all_but_emb  # Change in external resource abundance
    du.ind.X_emb = -dI_emb  # Change in vitellus (yolk)

    # remaining derivatives

    du.ind.A = dA(du.ind.I, p.ind.eta_IA, p.ind[:y_j][3]) # Assimilation flux
    du.ind.M = dM(ind.S, p.ind.k_M, p.ind[:y_j][2], p.ind[:y_T]) # Somatic maintenance flux
    du.ind.J = dJ(ind.H, p.ind.k_J, p.ind[:y_j][2], p.ind[:y_T]) # Maturity maintenance flux
    du.ind.S = dS(p.ind.kappa, du.ind.A, du.ind.M, p.ind.eta_SA, p.ind[:y_j][1], p.ind.eta_AS)
 
    du.ind.H = dH(ind.adult, p.ind.kappa, du.ind.A, du.ind.J) # maturiation flux
    du.ind.R = dR(ind.adult, p.ind.eta_AR, p.ind[:y_j][4], p.ind.kappa, du.ind.A, du.ind.J) # reproduction flux

    return nothing
end

"""
Individual-level part of the DEB-ODE model with arbitrary number of stressors, assuming IA to compute combined effects.
"""
@inline function DEBkiss_individual!(du, u, p, t)::Nothing
    
    TKTD_mix_IA!(du, u, p, t)
    DEBkiss_physiology!(du, u, p, t)

    return nothing
end

"""
DEB-ODE model with arbitrary number of stressors, assuming IA to compute combined effects. 
"""
function DEBODE!(du, u, p, t)
    DEBODE_global!(du, u, p, t)
    DEBkiss_individual!(du, u, p, t)
end