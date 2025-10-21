
# test09_paramstructs.jl
# using custom mutable structs for parameters instead of component vectors

p = deepcopy(EcotoxSystems.defaultparams)

abstract type AbstractParams end
using Distributions

begin # define default params as mutable struct

    const RealOrDist = Union{Real, Distribution}

    Base.@kwdef mutable struct GlobalParams <: AbstractParams
        N0::RealOrDist = 1                    # initial number of individuals [#]
        t_max::RealOrDist = 21.0              # maximum simulation time [days]
        dX_in::RealOrDist = 1200.0            # nutrient influx [μg C d⁻¹]
        k_V::RealOrDist = 0.0                 # chemostat dilution rate [d⁻¹]
        V_patch::RealOrDist = 0.05            # volume of a patch [L]
        C_W::Vector{RealOrDist} = [0.0]       # external chemical concentrations [μg L⁻¹]
        T::RealOrDist = 293.15                # ambient temperature [K]
    end

    Base.@kwdef mutable struct PropagateZoom <: AbstractParams
        dI_max::RealOrDist = 1/3
        dI_max_emb::RealOrDist = 1/3
        X_emb_int::RealOrDist = 1
        H_p::RealOrDist = 1
        K_X::RealOrDist = 1
    end

    Base.@kwdef mutable struct SpeciesParams <: AbstractParams
        Z::RealOrDist = Dirac(1.0)                   # individual variability through zoom factor
        propagate_zoom::NamedTuple = ( # zoom propagation mapping
            dI_max = 1/3, 
            dI_max_emb = 1/3, 
            X_emb_int = 1.0, 
            H_p = 1.0, 
            K_X = 1.0
        )  
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

    Base.@kwdef mutable struct DefaultParams <: AbstractParams
        glb::GlobalParams = GlobalParams()
        spc::SpeciesParams = SpeciesParams()
    end

    defaultparams = DefaultParams()

end


function generate_individual_params(spc::AbstractParams) 

    spc = deepcopy(spc)

    # sample from distributions 
    param_names = fieldnames(typeof(spc))
    [setproperty!(spc, param, EcotoxSystems.getval.(getproperty(spc, param))) for param in param_names]

    # propagate zoom factor
    for field in keys(spc.propagate_zoom)
        setproperty!(spc, field, getproperty(spc, field) * spc.Z .^ spc.propagate_zoom[field])
    end

    return spc

end

@benchmark begin 
    p = DefaultParams(spc = SpeciesParams(Z =  Truncated(Normal(1, 0.1), 0, 1), KD = [Truncated(Normal(1, 0.1), 0, Inf); ]))
    p.spc = generate_individual_params = p.spc
end

@benchmark begin
    p = deepcopy(EcotoxSystems.defaultparams)
    p.spc = EcotoxSystems.generate_individual_params(p)
end


Base.@kwdef mutable struct SomeParams
    aparams::AParams
    bparams::BParams
end

x = 1.
@benchmark $x += p.spc[:dI_max]

mutable struct BarPar
    dI_max::Float64
end

mutable struct FooPar
    spc::BarPar
end

p = FooPar(BarPar(1.))

x = 1.
@benchmark $x += p.spc.dI_max