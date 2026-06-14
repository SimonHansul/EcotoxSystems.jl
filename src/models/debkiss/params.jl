# defaultparams.jl 
# default parameter sets for the default model
# the default parameters are a plausible starting point for *Daphnia* with model currency μgC, 
# but their main purpose is to serve as a reference during development and testing (including extensions of the base model)

# Global parameters with ComponentVector
debkiss_global_params() = ComponentVector(
    t_max = 21.0,            # maximum simulation time [days]
    dX_in = 1200.0,          # nutrient influx [μg C d-1]
    k_V = 0.0,               # chemostat dilution rate [d-1]
    V_patch = 0.05,          # volume of a patch [L]
    T = 293.15,              # ambient temperature [K]
    C_W1 = 0.,               # aqueous concentration of chemical 1 [e.g. μg L^-1]
    C_W2 = 0.,               # aqueous concentration of chemical 2 [e.g. μg L^-1]
)

# Species-level DEB and TKTD parameters with ComponentVector
# TODO: harmonize parameter naming with debkiss literature
# TODO: add density of structure, set to 1 by default
debkiss_species_params() = ComponentVector(
    Z = Dirac(1.0), # mass-based zoom factor [-]
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
    # TKTD parameters for two substance and four standard PMoAs
    TKTD = ComponentVector( 
        k_D1_G = 0., # dominant rate constant PMoA growth efficiency [d^-1]
        k_D1_M = 0., # PMoA maintenance costs
        k_D1_A = 0., # PMoA assimilation efficiency
        k_D1_R = 0., # PMoA reproduction efficiency

        e1_G = 1., # sensitivity parameter (ED50); e.g. [μg L^-1]
        e1_M = 1.,
        e1_A = 1.,
        e1_R = 1.,

        b1_G = 1., # slope parameter [-]
        b1_M = 1.,
        b1_A = 1.,
        b1_R = 1.,
        
        k_D2_G = 0.,
        k_D2_M = 0.,
        k_D2_A = 0.,
        k_D2_R = 0.,

        e2_G = 1.,
        e2_M = 1.,
        e2_A = 1.,
        e2_R = 1.,

        b2_G = 1.,
        b2_M = 1.,
        b2_A = 1.,
        b2_R = 1.,
    ),
    aux = ComponentVector( # addition parameters used in the IBM part
        a_max = Truncated(Normal(60, 6), 0, Inf), # maximum life span 
        tau_R = 2.0, # reproduction interval  
        S_rel_crit = 0.5, # critical fraction of structure that can be lost before starvation hazard rate sets in [-]
        h_S = 1e-3, # starvation hazar, showvalues = [("N",N)])d rate [d^-1]
    )
)

defaultparams() = ComponentVector{Union{Real,Distribution}}(
    glb = debkiss_global_params(),
    spc = debkiss_species_params()
)


# the getval function makes it possible that any parameter can also be a distribution
# package users shouldn't have to diretly interact with this
getval(x::Distribution) = rand(x) 
getval(x::Any) = x

# propagate_zoom applies the zoomfactor to the parameters which are indicated in propagate_zoom, with the appropriate scaling factor
propagate_zoom(ind::ComponentVector) = begin
    zprop = ind[keys(ind.propagate_zoom)] .* ind.Z .^ ind.propagate_zoom
    ind = ComponentVector(
        ind;
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
    ind = getval.(p.spc) |> propagate_zoom 
    # type must not bee to specific for autodiff compatibility
    return ComponentVector(
        glb = p.glb, 
        ind = ind;
        kwargs...
    )
end

function generate_individual_params(p::ComponentVector; kwargs...)

    ind = getval.(p.spc) |> 
    x -> Real.(x) |> 
    x -> begin
        x.dI_max *= x.Z^(1/3)
        x.dI_max_emb *= x.Z^(1/3)
        x.H_p *= x.Z
        x.X_emb_int *= x.Z
        x.K_X *= x.Z^(1/3)
        x
    end

    return ComponentVector(
        glb = p.glb, 
        ind = ind; 
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
link_params!(p, (spc = link_ind_params!,)) # apply link ahead of simulation 
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

link_params!(::ComponentVector, ::Nothing)::Nothing = nothing

"""
Convenience function to update a TKTD parameter.

## args

- `p_tktd::ComponentVector`: A parameter Vector that contains TKTD parameters, with parameters named according to a fixed pattern: `[par][chemical_index]_[PMoA_label]`.             
- `par::Symbol` is the name of the parameter (e.g. `:k_D`).
- `chemical_index::Int` is the number of the component in the mixture (1, 2, ..n).
- `PMoA_label` is the abbreviation of the PMoA (e.g. `:G`, `:A`, `:A`, `:R`)

## Examples 

```
debkiss = FullDEBkiss()
p = debkiss.parameters # contains a component called `spc.TKTD` 
set_TKTD_param!(p, )
```
"""
function set_TKTD_param!(
    p_tktd::ComponentVector, 
    par::Union{Symbol,String}, 
    chemical_index::Int, 
    PMoA_label::Union{Symbol,String}, 
    value::Real
    )::Nothing

    parlabel = Symbol("$par$(chemical_index)_$(PMoA_label)")
    setindex!(p_tktd, value, parlabel)

    return nothing
end