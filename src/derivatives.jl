#derivatives.jl
# derivatives of the default model(s)

"""
function sig(
    x::Real, 
    x_thr::Real,
    y_left::Real, 
    y_right::Real; 
    beta::Real = 30
    )::Real

Sigmoid switch function. 
This can be useful to replace simple if-statements with a continuous function. 

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
Clip negative values as a continuous function.
"""
function clipneg(x::Real)::Real
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


cb_juvenile = ContinuousCallback(condition_juvenile, effect_juvenile!)
cb_adult = ContinuousCallback(condition_adult, effect_adult!)

DEBODE_callbacks = CallbackSet(cb_juvenile, cb_adult)


@inline function DEBODE_global!(du, u, p, t)::Nothing

    du.glb.X_p = p.glb.Xdot_in - p.glb.k_V * u.glb.X_p  

    return nothing
end

"""
Individual-level part of the DEB-ODE model with arbitrary number of stressors, assuming IA to compute combined effects.
"""
@inline function DEBODE_individual!(du, u, p, t)::Nothing

    @unpack glb, ind = u

    # scaled damage dynamics based on the minimal model

    du.ind.D_G = @. (1 - ind.embryo) * p.ind.k_D_G * (glb.C_W - ind.D_G)
    du.ind.D_M = @. (1 - ind.embryo) * p.ind.k_D_M * (glb.C_W - ind.D_M)
    du.ind.D_A = @. (1 - ind.embryo) * p.ind.k_D_A * (glb.C_W - ind.D_A)
    du.ind.D_R = @. (1 - ind.embryo) * p.ind.k_D_R * (glb.C_W - ind.D_R)
    du.ind.D_h = @. (1 - ind.embryo) * p.ind.k_D_h * (glb.C_W - ind.D_h)

    # response to chemical stressors based on independent action

    ind.y_G = prod(@. softNEC2neg(ind.D_G, p.ind.e_G, p.ind.b_G))  
    ind.y_M = prod(@. softNEC2pos(ind.D_M, p.ind.e_M, p.ind.b_M))  
    ind.y_A = prod(@. softNEC2neg(ind.D_A, p.ind.e_A, p.ind.b_A))  
    ind.y_R = prod(@. softNEC2neg(ind.D_R, p.ind.e_R, p.ind.b_R))
    ind.h_z = sum(@. softNEC2GUTS(ind.D_h, p.ind.e_h, p.ind.b_h)) 
    du.ind.S_z = -ind.h_z * ind.S_z

    # calculation of temperature correction factor

    ind.y_T = exp((p.ind.T_A / p.ind.T_ref) - (p.ind.T_A / p.glb.T)) 

    # ingestion rates and feedback with resource pools

    ind.f_X = (glb.X_p / p.glb.V_patch) / ((glb.X_p / p.glb.V_patch) + p.ind.K_X)
    
    dI_emb = ind.embryo * (Complex(ind.S)^(2/3)).re * p.ind.Idot_max_rel_emb * ind.y_T
    dI = (1-ind.embryo) * ind.f_X * p.ind.Idot_max_rel * (Complex(ind.S)^(2/3)).re * ind.y_T
    
    du.ind.I = dI_emb + dI

    du.glb.X_p -= dI  # Change in external resource abundance
    du.ind.X_emb = -dI_emb  # Change in vitellus

    # remaining derivatives

    du.ind.A = du.ind.I * p.ind.eta_IA * ind.y_A # Assimilation flux
    du.ind.M = ind.S * p.ind.k_M * ind.y_M * ind.y_T # Somatic maintenance flux
    du.ind.J = ind.H * p.ind.k_J * ind.y_M * ind.y_T #  
    du.ind.S = sig(p.ind.kappa * du.ind.A, du.ind.M, -(du.ind.M / p.ind.eta_SA - p.ind.kappa * du.ind.A), ind.y_G * p.ind.eta_AS * (p.ind.kappa * du.ind.A - du.ind.M))    
    du.ind.H = (1 - ind.adult) * clipneg(((1 - p.ind.kappa) * du.ind.A) - du.ind.J)
    du.ind.R = ind.adult * clipneg(p.ind.eta_AR * ind.y_R * ((1 - p.ind.kappa) * du.ind.A - du.ind.J))  # reproduction for adults

    return nothing
end

"""
DEB-ODE model with arbitrary number of stressors, assuming IA to compute combined effects. 
"""
function DEBODE!(du, u, p, t)

    DEBODE_global!(du, u, p, t)
    DEBODE_individual!(du, u, p, t)
end

