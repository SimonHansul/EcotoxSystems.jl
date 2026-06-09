#=
Simulate mixtures with different PMoAs
=#

using DataFramesMeta
import ComponentArrays: ComponentVector; const CV = ComponentVector

debkiss = FullDEBkiss()

# simulating two stressors with different PMoAs

p = DEBkiss.debkiss_defaultparams()
DEBkiss.set_TKTD_param!(p.spc.TKTD, :k_D, 1, "G", 1.)
DEBkiss.set_TKTD_param!(p.spc.TKTD, :k_D, 2, "R", 1.)


# simulating a ray design
Cvec_1 = [0, 0.5, 0.0, 0.5]
Cvec_2 = [0, 0.0, 0.5, 0.0]

p.spc.TKTD.b1_G

sim = DEBkiss.simulate_constant_exposure(
    p, 
    (:C_W1, :C_W2) => [
        0   0; 
        0.5 0; 
        0   0.5; 
        0.5 0.5
    ]
)

@df sim plot(
    plot(:t, :S, group = :treatment_id, leg = true),
    plot(:t, :R, group = :treatment_id)
)

# simulating stressors with identical PMoAs
# we chose both stressors to act via PMoA R, 
# so there should be no effects on growth in any of the treatments, 
# but every treatment should affect repro

p = DEBkiss.debkiss_defaultparams()
DEBkiss.set_TKTD_param!(p.spc.TKTD, :k_D, 1, "R", 1.)
DEBkiss.set_TKTD_param!(p.spc.TKTD, :k_D, 2, "R", 1.)


sim = DEBkiss.simulate_constant_exposure(
    p, 
    (:C_W1, :C_W2) => [
        0   0; 
        0.5 0; 
        0   0.5; 
        0.5 0.5
    ]
)

@df sim plot(
    plot(:t, :S, group = :treatment_id, leg = true),
    plot(:t, :R, group = :treatment_id)
)

# TODO: add proper test