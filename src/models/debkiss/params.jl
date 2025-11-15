# defaultparams.jl 
# default parameter sets for the default model
# the default parameters are a plausible starting point for *Daphnia* with model currency μgC, 
# but their main purpose is to serve as a reference during development and testing (including extensions of the base model)

# Global parameters with ComponentVector
debkiss_global_params = ComponentVector(
    N0 = 1,                    # initial number of individuals [#]
    t_max = 21.0,              # maximum simulation time [days]
    dX_in = 1200.0,          # nutrient influx [μg C d-1]
    k_V = 0.0,                 # chemostat dilution rate [d-1]
    V_patch = 0.05,            # volume of a patch [L]
    C_W = [0.0],               # external chemical concentrations [μg L-1]
    T = 293.15                 # ambient temperature [K]
)

# Species-level DEB and TKTD parameters with ComponentVector
debkiss_species_params = ComponentVector(
    Z = Dirac(1.0), # individual variability through zoom factor
    propagate_zoom = ComponentVector( # lists parameters which are affected by the zoom factor and the corresponding scaling exponent
        dI_max = 1/3, 
        dI_max_emb = 1/3,
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
    KD = Float64[0. 0. 0. 0.;], # KD - value per PMoA (G,M,A,R) and stressor (1 row = 1 stressor)
    B = Float64[2. 2. 2. 2.;], # slope parameters
    E = Float64[1e10 1e10 1e10 167;], # sensitivity parameters (thresholds)
    KD_h = Float64[0.;], # KD - value for GUTS-SD module (1 row = 1 stressor)
    E_h = Float64[1e10;], # sensitivity parameter (threshold) for GUTS-SD module
    B_h = Float64[2.;], # slope parameter for GUTS-SD module 
    # these are curently only used in an individual-based context, but could find application in the pure-ODE implementation 
    # for example by triggering emptying of the reproduction buffer through callbacks
    W_S_rel_crit = 0.66,  # relative amount of structure which can be lost before hazard rate kicks in
    h_S = 0.7, # starvation hazard rate caused be shrinking below W_S_rel_crit
    a_max = Truncated(Normal(60, 6), 0, Inf), # maximum life span 
    tau_R = 2.0, # reproduction interval
    # additional component for intermediate quantities, 
    # these are not technically model parameters (although they could be treated so), 
    # but this is nevertheless the best place to put them
    intermediates = ComponentVector( 
        y_j = ones(4), 
        h_z = 0.0
    )
)

default_debkiss_params = ComponentVector(
    glb = debkiss_global_params,
    spc = debkiss_species_params
)


# the getval function makes it possible that any parameter can also be a distribution
# package users shouldn't have to diretly interact with this
getval(x::Distribution) = rand(x) 
getval(x::Any) = x

# propagate_zoom applies the zoomfactor to the parameters which are indicated in propagate_zoom, with the appropriate scaling factor
propagate_zoom(ind::CVOrParamStruct) = begin
    zprop = ind[keys(ind.propagate_zoom)] .* ind.Z .^ ind.propagate_zoom
    ind = ComponentVector(
        ind;
        zprop...
    )
end

"""
    generate_individual_params(p::CVOrParamStruct; kwargs...)

    Generate individual-specific parameter set from species-specific parameter set. 
If a parameter entry is a distribution, a random sample is taken. 

This also works for Vectors of distributions. 
The kwargs need to be supplemented with additional components if there are more than just a global and an individual-level component.
"""
function generate_individual_params(p::CVOrParamStruct; kwargs...)::CVOrParamStruct
    ind = getval.(p.spc) |> propagate_zoom
    return ComponentVector(
        glb = p.glb, 
        ind = ind;
        kwargs...
    )
end

"""
    link_params!(p::CVOrParamStruct, links::NamedTuple = (spc = linkfun,))::Nothing 

Apply functions to link parameters. <br>

- `p`: ComponentVector containing parameters
- `links`: Component-specific functions to define links between parameters. 

## Examples 

We can apply the link before running a simulation, in which case the link is applied to the `spc` component 
(species-level parameters).-

```Julia 


link_ind_params!(p) = begin
    p.k_J_emb = (1-p.kappa_emb)/p.kappa_emb * p.k_M_emb
end

p = deepcopy(defaultparams)
link_params!(p, (spc = link_ind_params!,) # apply link ahead of simulation 
```

Alternatively, we provide the links as a keyword argument to `ODE_simulator`, 
in which case the link is applied to the `ind` component. <br>
This is useful if one of the parameters involved in the link is subject to individual variability, 
and we need to update the link for each simulated individual. <br>

```Julia

sim = ODE_simulator(p, param_links = (ind = link_ind_params!,))
```
"""
link_params!(p::CVOrParamStruct, links::NamedTuple = (spc = link_ind_params!,))::Nothing = begin

    for (component,linkfun!) in pairs(links)
        linkfun!(p[component])
    end

    return nothing
end

link_params!(p::CVOrParamStruct, links::Nothing)::Nothing = nothing