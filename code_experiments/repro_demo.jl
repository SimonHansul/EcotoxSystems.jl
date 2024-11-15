# demo of simulating repro explicitly in a life-table experiment 

using Pkg; Pkg.activate("test")

using Plots, StatsPlots
using Distributions
using DataFramesMeta
using Plots.Measures

using Revise
@time using EcotoxSystems
import EcotoxSystems: defaultparams, IBM_simulator
import EcotoxSystems: DEBkiss!, default_global_rules!

default(leg = false, thickness_scaling = 1.2)


#=
The rules for individuals are extended 
=#

import EcotoxSystems: default_individual_rules!
function juvenile_removal!(a, m)
    default_individual_rules!(a, m)

    if (a.u.ind.id>1)&&(isapprox(0, m.t%5, atol = m.dt))
        a.u.ind.cause_of_death = 3.
    end
end

#=

=#
begin
    p = copy(defaultparams)

    p.glb.dX_in = 1e7
    p.glb.k_V = 0.1
    p.glb.V_patch = 1
    p.glb.N0 = 1
    p.glb.t_max = 21
    
    p.spc.Z = Truncated(Normal(1, 0.05), 0, Inf)
    p.spc.tau_R = Truncated(Normal(2., 0.2), 0, Inf)
    p.spc.eta_IA = 0.33
    p.spc.f_Xthr = 0.9
    p.spc.s_min = 0.25
    p.spc.H_p = 100.
    p.spc.a_max = 60.
    p.spc.kappa = 0.8
    

    @time sim_ibm = IBM_simulator(
        p,
        dt = 1/24, 
        saveat = 1/24,
        showinfo = 7,
        individual_rules! = juvenile_removal!
        )

    plt = @df sim_ibm.spc plot(
        plot(:t, :S ./ maximum(:S), group = :id, leg = false, color = 1, xlabel = "Time (d)", ylabel = "Scaled dry weight"), 
        plot(:t, :R, group = :id, leg = false, color = 1, xlabel = "Time (d)", ylabel = "dR (μgC)"), 
        size = (600,600), bottommargin = 5mm, leftmargin = 5mm, layout = (2,1)
    )

    annotate!(subplots = 1, 7.5, 0.75, Plots.text("Mother individual", :blue))
    annotate!(subplots = 1, 15, 0.3, Plots.text("Offspring \n (removed in discrete intervals)", :green, 12))

    display(plt)

    savefig(plot(plt, dpi = 300), "repro_demo1.png")
end



#=

Simulating different food levels
=#

using ComponentArrays


import EcotoxSystems: combine_outputs, add_idcol

function param_sweep(
    simulator::Function,
    p::ComponentVector,
    component::Symbol,
    param::Symbol,
    vals::Vector{R}
    ) where R <: Real

    let val_int = p[component][param] # we will modify this value and then reset to the initial value
        sim = []

        for (i,val) in enumerate(vals)
            p[component][param] = val
            sim_i = simulator(p)
            sim_i = add_idcol(sim_i,param, val)
            sim_i = add_idcol(sim_i, :parval_id, i)
            push!(sim, sim_i)
        end
        
        setproperty!(p[component], param, val_int) 

        return combine_outputs(Vector{typeof(sim[1])}(sim); idcol = :treatment_id)
    end
end

import EcotoxSystems: groupedlineplot


#=
Simulating different food levels with daily food addition
=#


global_ODE_none!(du, u, p, t) = return nothing

# simulating daily medium renewal + addition of food
function food_renewal!(m)
    default_global_rules!(m)
    
    if isapprox(0, m.t%1, atol = m.dt)
        m.u.glb.X = m.p.glb.dX_in

        println(m.t)
    end
end



simulate_experiment(p) = IBM_simulator(
    p,
    dt = 1/24, 
    saveat = 1/24,
    showinfo = 7,
    individual_rules! = juvenile_removal!, 
    global_ode! = global_ODE_none!,
    global_rules! = food_renewal! 
    )

using Test
begin 
    p = copy(defaultparams)

    p.glb.dX_in = 1200

    p.spc.Z = 1. #Truncated(Normal(1, 0.1), 0, Inf)
    p.spc.H_p = 20
   
    p.spc.s_min = 1
    p.spc.tau_R = 2
    sim1 = simulate_experiment(p)

    plt = plot(
        (@df sim1.glb plot(:t, :X, xlabel = "t", ylabel = "Dry mass", title = "Food abundance over time")),
        (@df sim1.spc plot(:t, :S, group = :id, color = 1, xlabel = "t", ylabel = "Organism dry mass", label = "Simulation of offspring", title = "Organism growth"))
    )
    
    p.spc.tau_R = Inf
    
    sim2 = simulate_experiment(p)
    @df sim2.glb plot!(plt, subplot = 1, :t, :X, xlabel = "t", ylabel = "Food mass", title = "Daily food addition \n + renewals")
    @df sim2.spc plot!(plt, subplot = 2, :t, :S, group = :id, xlabel = "t", color = 2, label = "No simulation of offspring")
    
    annotate!(plt, 15, 100, Plots.text("Offspring \n simulated", palette(:default)[1], 10), subplot = 2)
    annotate!(plt, 5, 150, Plots.text("Offspring \n not simulated", palette(:default)[2], 10), subplot = 2)

    plt = plot(plt, size = (700,400), leg = false, xlabel = "Time")
    display(plt)

    savefig(plot(plt, dpi = 300), "repro_demo2.png")
end


begin
    p = copy(defaultparams)

    p.glb.dX_in = 1e7
    p.glb.k_V = 0.1
    p.glb.V_patch = 1
    p.glb.N0 = 1
    p.glb.t_max = 21
    
    p.spc.Z = Truncated(Normal(1, 0.05), 0, Inf)
    p.spc.tau_R = Inf #Truncated(Normal(2., 0.2), 0, Inf)
    p.spc.eta_IA = 0.33
    p.spc.f_Xthr = 0.9
    p.spc.s_min = 0.25
    p.spc.H_p = 100.
    p.spc.a_max = 60.
    p.spc.kappa = 0.8


    @time sim_ibm = IBM_simulator(
        p,
        dt = 1/24, 
        saveat = 1/24,
        showinfo = 7,
        individual_rules! = juvenile_removal!
        )

    plt = @df sim_ibm.spc plot(
        plot(:t, :S ./ maximum(:S), group = :id, leg = false, color = 1, xlabel = "Time (d)", ylabel = "Scaled dry weight"), 
        plot(:t, trunc.(:R ./ p.spc.X_emb_int), group = :id, leg = false, color = 1, xlabel = "Time (d)", ylabel = "cR (# offspring)"), 
        size = (600,600), bottommargin = 5mm, leftmargin = 5mm, layout = (2,1)
    )

    #annotate!(subplots = 1, 7.5, 0.75, Plots.text("Mother individual", :blue))
    #annotate!(subplots = 1, 15, 0.3, Plots.text("Offspring \n (removed in discrete intervals)", :green, 12))

    display(plt)

    savefig(plot(plt, dpi = 300), "repro_demo3.png")
end



begin
    p = copy(defaultparams)

    p.glb.dX_in = 1e7
    p.glb.k_V = 0.1
    p.glb.V_patch = 1
    p.glb.N0 = 1
    p.glb.t_max = 21
    
    p.spc.Z = Truncated(Normal(1, 0.05), 0, Inf)
    p.spc.tau_R = Inf #Truncated(Normal(2., 0.2), 0, Inf)
    p.spc.eta_IA = 0.33
    p.spc.f_Xthr = 0.9
    p.spc.s_min = 0.25
    p.spc.H_p = 100.
    p.spc.a_max = 60.
    p.spc.kappa = 0.8


    @time sim_ibm = IBM_simulator(
        p,
        dt = 1/24, 
        saveat = 1/24,
        showinfo = 7,
        individual_rules! = juvenile_removal!
        )

    plt = @df sim_ibm.spc plot(
        plot(:t, :S ./ maximum(:S), group = :id, leg = false, color = 1, xlabel = "Time (d)", ylabel = "Scaled dry weight"), 
        plot(:t, :R, group = :id, leg = false, color = 1, xlabel = "Time (d)", ylabel = "R (μg C)"), 
        size = (600,600), bottommargin = 5mm, leftmargin = 5mm, layout = (2,1)
    )

    #annotate!(subplots = 1, 7.5, 0.75, Plots.text("Mother individual", :blue))
    #annotate!(subplots = 1, 15, 0.3, Plots.text("Offspring \n (removed in discrete intervals)", :green, 12))

    display(plt)

    savefig(plot(plt, dpi = 300), "repro_demo4.png")
end