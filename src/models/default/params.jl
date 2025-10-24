# defaultparams.jl 
# default parameter sets for the default model
# the default parameters are a plausible starting point for *Daphnia* with model currency μgC, 
# but their main purpose is to serve as a reference during development and testing (including extensions of the base model)

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
    tau_R = 2.0 # reproduction interval
)

"""
Default parameter object
"""
defaultparams = ComponentVector(
    glb = global_params,   
    spc = species_params     
)

Params(;kwargs...) = ComponentVector(defaultparams; kwargs...) # for backward compat with old struct-based implementation

"""
    params(;kwargs...) = ComponentVector(defaultparams; kwargs...)

Initialize the default parameter set, modifying and/or adding parameter with kwargs. 

## Examples 

```Julia

# simulate the default parameters
p = params()
sim = ODE_simulator(p)  

# simulate the default parameters, but modify kappa
p = params(kappa = 0.5) 
sim = ODE_simulator

# simulate the default parameters with individual variabilty in kappa
p = params(kappa = truncated(Normal(0.5, 0.1), 0, 1))
sim = @replicates ODE_simulator(p) 10
``` 

"""
function params(;kwargs...)
    defparams = deepcopy(defaultparams)
    return ComponentVector(defparams; kwargs...)
end

# the getval function makes it possible that any parameter can also be a distribution
# package users shouldn't have to diretly interact with this
getval(x::Distribution) = rand(x) 
getval(x::Any) = x

# propagate_zoom applies the zoom factor to the parameters which are indicated in propagate_zoom, with the given value as exponent
propagate_zoom(spc::ComponentVector) = begin
    zprop = spc[keys(spc.propagate_zoom)] .* spc.Z .^ spc.propagate_zoom
    spc = ComponentVector(
        spc;
        zprop...
    )
end

"""
    generate_individual_params(p::ComponentVector; kwargs...)

Generate individual-specific parameter set from species-specific parameter set. 
If a parameter entry is a distribution, a random sample is taken. 

This also works for Vectors of distributions. 
The kwargs need to be supplemented with additional components if there are more than just a global and an individual-level component.
"""
function generate_individual_params(p::ComponentVector; kwargs...)::ComponentVector
    spc = getval.(p.spc) |> propagate_zoom
    return ComponentVector(
        glb = p.glb, 
        spc = spc;
        kwargs...
    )
end

"""
    link_params!(p::ComponentVector, links::NamedTuple = (spc = linkfun,))::Nothing 

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
link_params!(p::ComponentVector, links::NamedTuple = (spc = link_ind_params!,))::Nothing = begin

    for (component,linkfun!) in pairs(links)
        linkfun!(p[component])
    end

    return nothing
end

link_params!(p::ComponentVector, links::Nothing)::Nothing = nothing


# functions for using mutable structs as params


abstract type AbstractParams end
abstract type AbstractSpeciesParams <: AbstractParams end
abstract type AbstractGlobalParams <: AbstractParams end
abstract type AbstractComponentParams <: AbstractParams end
abstract type AbstractAuxParams <: AbstractParams end
abstract type AbstractParamEnsemble <: AbstractParams end

using Distributions

const RealOrDist = Union{Real, Distribution}

Base.@kwdef mutable struct GlobalParams <: AbstractGlobalParams
    N0::RealOrDist = 1                    # initial number of individuals [#]
    t_max::RealOrDist = 21.0              # maximum simulation time [days]
    dX_in::RealOrDist = 1200.0            # nutrient influx [μg C d⁻¹]
    k_V::RealOrDist = 0.0                 # chemostat dilution rate [d⁻¹]
    V_patch::RealOrDist = 0.05            # volume of a patch [L]
    C_W::Vector{RealOrDist} = [0.0]       # external chemical concentrations [μg L⁻¹]
    T::RealOrDist = 293.15                # ambient temperature [K]
end

Base.@kwdef mutable struct PropagateZoom <: AbstractAuxParams
    dI_max::RealOrDist = 1/3
    dI_max_emb::RealOrDist = 1/3
    X_emb_int::RealOrDist = 1
    H_p::RealOrDist = 1
    K_X::RealOrDist = 1
end

Base.@kwdef mutable struct SpeciesParams <: AbstractSpeciesParams
    Z::RealOrDist = Dirac(1.0)                   # individual variability through zoom factor
    propagate_zoom::PropagateZoom = PropagateZoom() # zoom factor propagation exponents  
    T_A::RealOrDist = 8000.0                    # Arrhenius temperature [K]
    T_ref::RealOrDist = 293.15                  # reference temperature [K]
    X_emb_int::RealOrDist = 19.42               # initial vitellus [μgC]
    K_X::RealOrDist = 1000.0                    # half-saturation constant for food uptake [μgC L⁻¹]
    dI_max::RealOrDist = 22.9                   # max size-specific ingestion rate [μgC μgC⁻(2/3) d⁻¹]
    dI_max_emb::RealOrDist = 22.9               # embryonic ingestion rate [μgC μgC⁻(2/3) d⁻¹]
    kappa::RealOrDist = 0.539                   # somatic allocation fraction [-]
    eta_IA::RealOrDist = 0.33                   # assimilation efficiency [-]
    eta_AS::RealOrDist = 0.8                    # growth efficiency [-]
    eta_SA::RealOrDist = 0.8                    # shrinking efficiency [-]
    eta_AR::RealOrDist = 0.95                   # reproduction efficiency [-]
    k_M::RealOrDist = 0.59                      # somatic maintenance rate constant [d⁻¹]
    k_J::RealOrDist = 0.504                     # maturity maintenance rate constant [d⁻¹]
    H_p::RealOrDist = 100                       # maturity at puberty [μgC]

    KD::Matrix{RealOrDist} = [0. 0. 0. 0.;]     # KD per PMoA (G,M,A,R) × stressor
    B::Matrix{RealOrDist} = [2. 2. 2. 2.;]      # slope parameters
    E::Matrix{RealOrDist} = [1e10 1e10 1e10 167.;] # sensitivity thresholds

    KD_h::Vector{RealOrDist} = [0.0]            # KD for GUTS-SD
    E_h::Vector{RealOrDist} = [1e10]            # sensitivity threshold for GUTS-SD
    B_h::Vector{RealOrDist} = [2.0]             # slope parameter for GUTS-SD

    W_S_rel_crit::RealOrDist = 0.66             # relative structure threshold before hazard
    h_S::RealOrDist = 0.7                       # starvation hazard rate
    a_max::RealOrDist = Truncated(Normal(60, 6), 0, Inf) # max life span
    tau_R::RealOrDist = 2.0                     # reproduction interval
end

Base.@kwdef mutable struct DefaultParams <: AbstractParamEnsemble
    glb::GlobalParams = GlobalParams()
    spc::SpeciesParams = SpeciesParams()
end

getval(x::Pair{Symbol,Distribution}) = x # no randomization for pair types
getval(x::PropagateZoom) = x # no randomization for zoom propagation exponents

function propagate_zoom!(spc::AbstractSpeciesParams, propagate_zoom::PropagateZoom)::Nothing
    for field in fieldnames(PropagateZoom)
        setproperty!(spc, field, getproperty(spc, field) * spc.Z .^ getproperty(spc.propagate_zoom, field))
    end
    return nothing
end

function generate_individual_params(spc::AbstractSpeciesParams) 

    T = typeof(spc) # identify type
    new_spc = T() # instantiate default parameter object of the given type

    for param in fieldnames(T) # iterate over all parameters
        val = getfield(spc, param) # get the current value - this might be a scalar or distribution
        indval = EcotoxSystems.getval(val) # infer the individual-level value

        setfield!(new_spc, param, indval) # update the value
    end

    # propagate zoom factor
    propagate_zoom!(new_spc, new_spc.propagate_zoom)

    return new_spc

end

function generate_individual_params(p::AbstractParamEnsemble) 
    p = typeof(p)(
        spc = generate_individual_params(p.spc)
    )
    return p
end

function link_params!(p::AbstractParams, links::NamedTuple = (spc = link_ind_params!,))::Nothing 

    for (component, linkfun!) in pairs(links)
        linkfun!(getproperty(p, component))
    end

    return nothing
end

link_params!(p::AbstractParams, links::Nothing) = nothing