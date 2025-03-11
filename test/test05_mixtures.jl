#=
Simulate mixtures with different PMoAs
=#

using DataFramesMeta
import EcotoxSystems: exposure, relative_response


# changing the number of stressors is curently a bit clunky
# we cannot do it dynamically, because ComponentArrays does not allow us to resize a Vector which is a field of a CA
# so the only solution currently is to manually reconstruct the CV

p = ComponentVector( 
    glb = ComponentVector(
        N0 = 1,                    # initial number of individuals [#]
        t_max = 21.0,              # maximum simulation time [days]
        dX_in = 1200.0,          # nutrient influx [μg C d-1]
        k_V = 0.0,                 # chemostat dilution rate [d-1]
        V_patch = 0.05,            # volume of a patch [L]
        C_W = [0., 0.],               # external chemical concentrations [μg L-1]
        T = 293.15                 # ambient temperature [K]
    ),
    spc = ComponentVector(
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
        H_p = 1 / 3,    # maturity at puberty [μgC]
        k_D_j = zeros((2, 4)), # k_D - value per PMoA (G,M,A,R) and stressor (1 row = 1 stressor)
        e_z = ones((2, 4)), # sensitivity parameters (thresholds)
        b_z = ones((2,4)), # slope parameters
        k_D_h = zeros((2,1)), # k_D - value for GUTS-Sd module (1 row = 1 stressor)
        e_h = ones((2,1)), # sensitivity parameter (threshold) for GUTS-SD module
        b_h = ones((2,1)), # slope parameter for GUTS-SD module 
        # these are curently only used in an individual-based context, but could find application in the pure-ODE implementation 
        # for example by triggering emptying of the reproduction buffer through callbacks
        f_Xthr = 0.5,  # functional response threshold for starvation mortality
        s_min = 0.25,  # daily survival mortality at complete food deprivation
        a_max = Truncated(Normal(60, 6), 0, Inf), # maximum life span 
        tau_R = 2.0 # reproduction interval
    )    
)


# simulating two stressors with different PMoAs

p.spc.k_D_j = [
    1. 0. 0. 0.; # stressor 1 has PMoA G
    0. 0. 0. 1. # stressor 2 has PMoA R
    ]

# setting all dose-response params to the same values
p.spc.e_z = ones(size(p.spc.k_D_j))
p.spc.b_z = fill(1, size(p.spc.k_D_j))

C_Wmat = [ # simulating a ray design
    0.0 0.0; 
    1.5 0.0;
    0.0 1.5;
    1.5 1.5;
]

p.spc.k_D_j 
p.spc.e_z
p.spc.b_z

sim = exposure(p->EcotoxSystems.ODE_simulator(p), p, C_Wmat);

@time ODE_simulator(p);


@df sim plot(
    plot(:t, :S, group = :treatment_id, leg = true),
    plot(:t, :R, group = :treatment_id)
)

# simulating stressors with identical PMoAs
# we chose both stressors to act via PMoA R, 
# so there should be no effects on growth in any of the treatments, 
# but every treatment should affect 

p.glb.C_W

p.spc.k_D_j = [
    0. 0. 0. 1.; 
    0. 0. 0. 1.
    ]

sim = exposure(ODE_simulator, p, C_Wmat);

@df sim plot(
    plot(:t, :S, group = :treatment_id, leg = true),
    plot(:t, :R, group = :treatment_id)
)
    
