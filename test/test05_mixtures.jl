#=
Simulate single stressors with different PMoAs
=#


using DataFramesMeta
import EnergyBudgetDiffEqs: exposure, relative_response

p = DEB.params()
p.spc.k_D_z .= 0

p.glb.t_max = 42.

p.spc.e_z .= 1.
p.spc.b_z .= 0.1

# changing the number of stressors is curently a bit clunky
# we cannot do it dynamically, because ComponentArrays does not allow us to resize a Vector which is a field of a CA
# so the only solution currently is to manually reconstruct the CA

p = ComponentVector( 
    glb = ComponentVector(
        N0 = 1,                    # initial number of individuals [#]
        t_max = 21.0,              # maximum simulation time [days]
        Xdot_in = 1200.0,          # nutrient influx [μg C d-1]
        k_V = 0.0,                 # chemostat dilution rate [d-1]
        V_patch = 0.05,            # volume of a patch [L]
        C_W = [0., 0.],               # external chemical concentrations [μg L-1]
        T = 293.15                 # ambient temperature [K]
    ),
    spc = DEB.species_params
)

# simulating two stressors with different PMoAs

p.spc.k_D_z = [
    1. 0. 0. 0.; # stressor 1 has PMoA G
    0. 0. 0. 1. # stressor 2 has PMoA R
    ]
# setting all dose-response params to the same values
p.spc.e_z = ones(size(p.spc.k_D_z))
p.spc.b_z = fill(1, size(p.spc.k_D_z))

C_Wmat = [ # simulating a ray design
    0.0 0.0; 
    1.5 0.0;
    0.0 1.5;
    1.5 1.5;
]

sim = exposure(DEB.simulator, p, C_Wmat);

@df sim plot(
    plot(:t, :S, group = :treatment_id, leg = true),
    plot(:t, :R, group = :treatment_id)
)


sim.y_z__1


# simulating stressors with identical PMoAs


p.glb.C_W

p.spc.k_D_z = [
    0. 0. 0. 1.; 
    0. 0. 0. 1.
    ]


C_Wmat

sim = exposure(DEB.simulator, p, C_Wmat);

@df sim plot(
    plot(:t, :S, group = :treatment_id, leg = true),
    plot(:t, :R, group = :treatment_id)
)
    

p.glb.C_W


@df sim plot(
    plot(:t, :y_j_4, group = :treatment_id)
)


p.glb.C_W
