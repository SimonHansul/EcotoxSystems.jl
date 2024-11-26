using ComponentArrays
using Distributions

# Global parameters with ComponentVector
global_params = ComponentVector(
    N0 = 1,                    # initial number of individuals [#]
    t_max = 21.0,              # maximum simulation time [days]
    dX_in = 1200.0,          # nutrient influx [μg C d-1]
    k_V = 0.0,                 # chemostat dilution rate [d-1]
    V_patch = 0.05,            # volume of a patch [L]
    C_W = [0.0],               # external chemical concentrations [μg L-1]
    T = 293.15                 # ambient temperature [K]
)

# Species-level DEB and TKTD parameters with ComponentVector
species_params = ComponentVector(
    Z = Dirac(1.0), # individual variability through zoom factor
    propagate_zoom = ComponentVector( # lists parameters which are affected by the zoom factor and the corresponding scaling exponent
        dI_max = 1/3, 
        Idot_max_rel_emb = 1/3,
        X_emb_int = 1,
        H_p = 1, 
        K_X = 1
    ),
    T_A = 8000.0,    # Arrhenius temperature [K]
    T_ref = 293.15,  # reference temperature [K]
    X_emb_int = 19.42, # initial vitellus [μgC]
    K_X = 1000.0, # half-saturation constant for food uptake [μgC L-1]
    dI_max = 22.9,      # maximum size-specific ingestion rate [μgC μgC^-(2/3) d-1]
    dI_max_emb = 22.9,  # size-specific embryonic ingestion rate [μgC μgC^-(2/3) d-1]
    kappa = 0.539,  # somatic allocation fraction [-]
    eta_IA = 0.33,  # assimilation efficiency [-]
    eta_AS = 0.8,   # growth efficiency [-]
    eta_SA = 0.8,  # shrinking efficiency [-]
    eta_AR = 0.95,  # reproduction efficiency [-]
    k_M = 0.59,     # somatic maintenance rate constant [d^-1]
    k_J = 0.504,    # maturity maintenance rate constant [d^-1]
    H_p = 100,    # maturity at puberty [μgC]
    k_D_z = [0 0 0 .38;], # k_D - value per PMoA (G,M,A,R) and stressor (1 row = 1 stressor)
    b_z = [0 0 0 0.93;], # slope parameters
    e_z = [0 0 0 167;], # sensitivity parameters (thresholds)
    k_D_h = [0;], # k_D - value for GUTS-Sd module (1 row = 1 stressor)
    e_h = [0;], # sensitivity parameter (threshold) for GUTS-SD module
    b_h = [0;], # slope parameter for GUTS-SD module 
    # these are curently only used in an individual-based context, but could find application in the pure-ODE implementation 
    # for example by triggering emptying of the reproduction buffer through callbacks
    f_Xthr = 0.5,  # functional response threshold for starvation mortality
    s_min = 0.25,  # daily survival mortality at complete food deprivation
    a_max = Truncated(Normal(60, 6), 0, Inf), # maximum life span 
    tau_R = 2.0 # reproduction interval
)


"""
Default parameter object
"""
defaultparams = ComponentVector(
    glb = global_params,   
    spc = species_params     
)

Params() = deepcopy(defaultparams) # for backward compat with old struct-based implementation
params() = deepcopy(defaultparams)

getval(x::Distribution) = rand(x) 
getval(x::Any) = x
#getval(x::Vector{Distribution}) = rand.(x)

propagate_zoom(ind::ComponentVector) = begin
    zprop = ind[keys(ind.propagate_zoom)] .* ind.Z .^ ind.propagate_zoom
    ind = ComponentVector(
        ind;
        zprop...
    )
end

"""
Generate individual-specific parameter set from species-specific parameter set. 
If a parameter entry is a distribution, a random sample is taken. 

This also works for Vectors of distributions.
"""
generate_individual_params(p::ComponentVector; kwargs...) = begin
    ind = getval.(p.spc) |> propagate_zoom
    return ComponentVector(
        glb = p.glb, 
        ind = ind;
        kwargs...
    )
end
