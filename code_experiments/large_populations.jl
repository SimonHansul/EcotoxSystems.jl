# sismulating a large population for several years for benchmarking and profiling

using Pkg; Pkg.activate("test")
using Revise
@time import EcotoxSystems: defaultparams, IBM_simulator
using Plots, StatsPlots
using Distributions


import EcotoxSystems: DEBkiss!, default_global_rules!

function modified_global_rules!(m)::Nothing
    
    default_global_rules!(m)
    
    if isapprox(0, m.t%14, atol = 1/24)
        println((t = m.t, m.u.glb.N))
    end
    
    return nothing
end

# modification of the individual-ODE without TK
function individual_noTK!(du, u, p, t)::Nothing
    
    @. u.ind.y_j = 1.
    u.ind.h_z = 0.
    
    DEBkiss!(du, u, p, t)
    
    return nothing
end

# rough comparison with netlogo
# netlogo took 1:30 for 56d, population leveling off around 900 individuals, dt=1/229
begin
    p = defaultparams

    p.glb.dX_in = 1e7
    p.glb.k_V = 0.1
    p.glb.V_patch = 1
    p.glb.N0 = 10
    p.glb.t_max = 3650

    p.spc.Z = Truncated(Normal(1, 0.05), 0, Inf)
    p.spc.tau_R = Truncated(Normal(2., 0.2), 0, Inf)
    p.spc.eta_IA = 0.33
    p.spc.f_Xthr = 0.9
    p.spc.s_min = 0.25
    p.spc.H_p = 100.
    p.spc.a_max = 60.

    @time sim_ibm = IBM_simulator(
        p,
        global_rules! = modified_global_rules!,
        individual_ode! = individual_noTK!,  
        dt = 1/24, 
        showinfo = 7, 
        saveat = 60
        )

    @df sim_ibm.glb plot(:t, :N)
end


