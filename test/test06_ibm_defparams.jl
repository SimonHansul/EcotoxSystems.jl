using Pkg; Pkg.activate("test")
Pkg.develop(path = ".")

#Pkg.rm("EcotoxSystems")
#Pkg.precompile()

using Distributions
using DataFramesMeta
using DataFrames
using Test
using Plots, StatsPlots

using Revise
using EcotoxSystems, EcotoxSystems.DEBkiss, EcotoxSystems.IBM

# re-writing the IBM simulator to avoid function wrapping altogether...

global_ode! = DEBkiss.constant_nutrient_influx!
global_rules! = DEBkiss.default_global_rules!

individual_ode! = DEBkiss.debkiss!
individual_rules! = DEBkiss.individual_rules!

using EcotoxSystems.ComponentArrays, EcotoxSystems.Parameters

# step 1: instantiate an individual ✅

p = DEBkiss.defaultparams()
u_glb = DEBkiss.initialize_global_statevars(p)

@time a = IBM.Individual(
    p, 
    u_glb, 
    DEBkiss.initialize_individual_statevars,
    DEBkiss.generate_individual_params
);

@testset "initialization of Individual expands the state vector" begin
    @test a.u.ind.cohort == 1
    @test a.u.ind.id == 1
end

@testset "accessing individual properties is allocation-free" begin

    # this caused ≈ 1200 allocations when `a` was a mutable struct
    # expecting 0 when a is `a` NamedTuple
    function alloc_test(a)
        for _ in 1:1000
            a.u.ind.id
        end
    end

    [alloc_test(a) for _ in 1:10];
    allocs = @allocated alloc_test(a) # 
    @test allocs == 0
end

# step 2: instantiate a model ✅
# step 3: put individuals into model ✅
# TODO: same for model as for individuals...=> NamedTuple instead of struct!

@testset "Instantiating a model and population" begin

    global m = IBM.Model(
        p; 
        init_u_glb = DEBkiss.initialize_global_statevars,
        init_u_ind = DEBkiss.initialize_individual_statevars,
        gen_p_ind = DEBkiss.generate_individual_params,
        dt = 1/24, 
        N0 = 10
    )

    @test length(m.individuals) == 10
    @test m.aux.idcount == 10
end

function test_ind_access_allocs(m)
    for _ in 1:1000
        m.individuals[1] # accesing individual should be allocation-free
        m.individuals[1].u.ind.id # accessing individual state should be allocation-free
    end
end

@allocated test_ind_access_allocs(m) # precompile
@test (@allocated test_ind_access_allocs(m)) == 0

# step 4: run an individual step ✅
# step 5: run a  model step ✅

individual_step! = IBM.make_individual_step(
    individual_ode!, 
    individual_rules!, 
    DEBkiss.initialize_individual_statevars,
    DEBkiss.generate_individual_params
    )

m = IBM.Model(
        p; 
        init_u_glb = DEBkiss.initialize_global_statevars,
        init_u_ind = DEBkiss.initialize_individual_statevars,
        gen_p_ind = DEBkiss.generate_individual_params,
        dt = 1/24, 
        N0 = 10_000
    );

model_step! = IBM.make_model_step(global_ode!, global_rules!, individual_step!);

theme(:juno)


begin # simulating a small simulation to check plausibility
    p.glb.t_max = 120.
    p.glb.dX_in = 50_000. 
    p.glb.k_V = 0.1
    p.glb.V_patch = 0.5

    p.spc.Z = Truncated(Normal(1, 0.1), 0, Inf)
    p.spc.aux.tau_R = 2
    p.spc.aux.h_S = 1.
    p.spc.aux.a_max = Truncated(Normal(30, 12), 0, Inf)
    p.spc.K_X = 500.

    @time sim = IBM.simulate(
        p; 
        global_ode! = DEBkiss.constant_nutrient_influx!,
        individual_ode! = DEBkiss.debkiss!,

        global_rules! = DEBkiss.default_global_rules!,
        individual_rules! = DEBkiss.individual_rules!,

        init_u_glb = DEBkiss.initialize_global_statevars,
        init_u_ind = DEBkiss.initialize_individual_statevars,
        gen_p_ind = DEBkiss.generate_individual_params,

        N0 = 10, 
        record_individuals = true, 
        saveat = 1/24, 
        dt = 1/24
    )

    @df sim.glb plot(
        plot(:t, :N), # why does population collapes entirely?
        plot(:t, :X) # why does X suddenly jump up, then increase linearly?
    )

    # proper test needed here
end



# the output for a single individual from the IBM is compared with the pure-ODE solution 
# this is done for a case with food limitation, since that is where there are most likely to be discrepancies
# (it is curently not possible to handle feedback with external food abundance in exactly the same way in both variants, and I am not sure it ever will be)
@testset "IBM vs ODE" begin

    debkiss = DEBkiss.FullDEBkiss() 
    p = debkiss.parameters
    
    p.glb.dX_in = 1e2
    p.glb.k_V = 0.1
    p.glb.V_patch = 0.05
    p.glb.t_max = 56

    p.spc.Z = Dirac(1) # don't simulate individual variability
    
    p.spc.H_p = 100.
    p.spc.K_X = 10
    
    p.spc.aux.tau_R = Inf # don't simulate reproduction for now - just tracking repro buffer
    p.spc.aux.a_max = Inf

    @time global sim_ode = DEBkiss.simulate(debkiss, saveat = 1);
    @time global sim_ibm = simulate_ibm(
        debkiss.parameters, 
        dt = 1/2400, N0 = 1, 
        record_individuals = true
        );

    plt = @df sim_ode plot(
        plot(:t, :X, xlabel = "t", ylabel = "X", leg = true, label = "ODE simulator"),
        plot(:t, :S, xlabel = "t", ylabel = "S", leg = true, label = "ODE simulator"),
        plot(:t, :H, xlabel = "t", ylabel = "H", leg = true, label = "ODE simulator"),
        plot(:t, :R, xlabel = "t", ylabel = "R", leg = true, label = "ODE simulator"),
        )
    
    @df sim_ibm.glb plot!(:t, :X, label = "IBM simulator", subplot = 1)
    @df sim_ibm.spc plot!(:t, :S, label = "IBM simulator", subplot = 2)
    @df sim_ibm.spc plot!(:t, :H, label = "IBM simulator", subplot = 3)
    @df sim_ibm.spc plot!(:t, :R, label = "IBM simulator", subplot = 4)
    display(plt)
    
    sim_ibm.spc.t = round.(sim_ibm.spc.t)
    
    eval_df = leftjoin(sim_ibm.spc, sim_ode, on = :t, makeunique = true) |> 
    x -> select(x, [:S, :S_1, :H, :H_1, :R, :R_1]) |> 
    x -> EcotoxSystems.dropmissing(x)
    
    dS = eval_df.S .- eval_df.S_1
    dH = eval_df.H .- eval_df.H_1
    dR = eval_df.R .- eval_df.R_1
   
    # proper test needed here
end

