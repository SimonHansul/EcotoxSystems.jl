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
    beta::Real = 30
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

condition_juvenile(u, t, integrator) = u.ind.X_emb # transition to juvenile when X_emb hits 0
function effect_juvenile!(integrator) 
    integrator.u.ind.embryo = 0.
    integrator.u.ind.juvenile = 1.
    integrator.u.ind.adult = 0.
end
 
condition_adult(u, t, integrator) = integrator.p.ind.H_p - u.ind.H # condition to adult when H reaches H_p
function effect_adult!(integrator) 
    integrator.u.ind.embryo = 0.
    integrator.u.ind.juvenile = 0.
    integrator.u.ind.adult = 1.
end

# putting the conditions together into a callback set

cb_juvenile = ContinuousCallback(condition_juvenile, effect_juvenile!, save_positions = (false,false))
cb_adult = ContinuousCallback(condition_adult, effect_adult!, save_positions = (false,false))

DEBODE_callbacks = CallbackSet(cb_juvenile, cb_adult)

@inline function DEBODE_global!(du, u, p, t)::Nothing

    du.glb.X = p.glb.dX_in - p.glb.k_V * u.glb.X  

    return nothing
end

"""
    function TKTD_mix_IA!(du, u, p, t)::Nothing

Mixture-TKTD for an arbitrary number of stressors, assuming Independent Action.
"""
@inline function TKTD_mix_IA!(du, u, p, t)::Nothing

    @unpack glb, ind = u

    # scaled damage dynamics based on the minimal model

    for z in eachindex(glb.C_W)
        # for sublethal effects, we broadcost over all PMoAs
        @. du.ind.D_j[z,:] = (1 - ind.embryo) * @view(p.ind.k_D_j[z,:]) * (glb.C_W[z] - @view(ind.D_j[z,:]))
        # for lethal effects, we have only one value per stressor
        du.ind.D_h[z] = (1 - ind.embryo) * p.ind.k_D_h[z] * (glb.C_W[z] - ind.D_h[z])
    end

    @. ind.y_z = softNEC2neg(ind.D_j, p.ind.e_z, p.ind.b_z) # relative responses per stressor and PMoA
    
    ind.y_j .= reduce(*, ind.y_z; dims=1) # relative responses per PMoA are obtained as the product over all chemical stressors
    ind.y_j[2] /= ind.y_j[2]^2 # for pmoas with increasing responses (M), the relative response has to be inverted  (x/x^2 == 1/x) 

    #ind.h_z = sum(@. softNEC2GUTS(ind.D_h, p.ind.e_h, p.ind.b_h)) # hazard rate according to GUTS-RED-SD
    ind.h_z = 0 
    @inbounds for z in eachindex(ind.D_h)
        ind.h_z += softNEC2GUTS(ind.D_h[z], p.ind.e_h[z], p.ind.b_h[z])
    end
    
    du.ind.S_z = -ind.h_z * ind.S_z # survival probability according to GUTS-RED-SD
    
    return nothing
end

"""
    DEBkiss!(du, u, p, t)::Nothing

Dynamics of DEBkiss model with maturity and explicit simulation of resource dynamics.

The density of structure is ignored, and instead `S^(2/3)` is applied for surface-area scaling. 
This affects the dimension of `dI_max`, but has no effect on the model dynamics.

If model output is to be compared to length data, a statistical weight-lenght relationship 
has to be applied to the model output.
"""
function DEBkiss!(du, u, p, t)::Nothing

    @unpack glb, ind = u

    ind.y_T = @fastmath exp((p.ind.T_A / p.ind.T_ref) - (p.ind.T_A / p.glb.T)) # temperature correction

    # ingestion rates and feedback with resource pools

    ind.f_X = @fastmath (glb.X / p.glb.V_patch) / ((glb.X / p.glb.V_patch) + p.ind.K_X)

    # calculation of resource uptake for embryos vs hatched individuals
    
    dI_emb = ind.embryo * (Complex(ind.S)^(2/3)).re * p.ind.dI_max_emb * ind.y_T
    dI = (1-ind.embryo) * ind.f_X * p.ind.dI_max * (Complex(ind.S)^(2/3)).re * ind.y_T

    # ingestion rate is the sum of both (dI_emb and dI are mutually exclusive)
    
    du.ind.I = dI_emb + dI

    du.glb.X -= dI  # Change in external resource abundance
    du.ind.X_emb = -dI_emb  # Change in vitellus (yolk)

    # remaining derivatives

    du.ind.A = du.ind.I * p.ind.eta_IA * ind.y_j[3] # Assimilation flux
    du.ind.M = ind.S * p.ind.k_M * ind.y_j[2] * ind.y_T # Somatic maintenance flux
    du.ind.J = ind.H * p.ind.k_J * ind.y_j[2] * ind.y_T # Maturity maintenance flux
    du.ind.S = sig( # Somatic growth
        p.ind.kappa * du.ind.A, 
        du.ind.M, -(du.ind.M / p.ind.eta_SA - p.ind.kappa * du.ind.A), 
        ind.y_j[1] * p.ind.eta_AS * (p.ind.kappa * du.ind.A - du.ind.M)
        )    
    du.ind.H = (1 - ind.adult) * clipneg(((1 - p.ind.kappa) * du.ind.A) - du.ind.J) # maturiation flux
    du.ind.R = ind.adult * clipneg(p.ind.eta_AR * ind.y_j[4] * ((1 - p.ind.kappa) * du.ind.A - du.ind.J))  # reproduction flux

    return nothing
end

"""
Individual-level part of the DEB-ODE model with arbitrary number of stressors, assuming IA to compute combined effects.
"""
@inline function DEBODE_individual!(du, u, p, t)::Nothing
    
    TKTD_mix_IA!(du, u, p, t)
    DEBkiss!(du, u, p, t)

    return nothing
end

"""
DEB-ODE model with arbitrary number of stressors, assuming IA to compute combined effects. 
"""
function DEBODE!(du, u, p, t)
    DEBODE_global!(du, u, p, t)
    DEBODE_individual!(du, u, p, t)
end