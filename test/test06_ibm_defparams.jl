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

DEBkiss.debkiss!
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
#   the example simulation does not work yet, but we are now in model-specification territory...

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


# step 6: put it all together


# FIXME: life stage indicators occasionally go through the roof...
# FIXME: individuals are not dying...?
#       - [x] food levels are not updated...
#           - fixed: forgot to include _set_global_states!()
#       - [x] individual not instantiated correctly
#       - [x] still no death.
#           - [x] do inds lose struct? yes
#           - [x] is cause_of_death set correctly?  no
#           - [x] aging submodule: age does not change correctly
#           - [x] starvation submodule...
#       - [ ] why does the population crash?
#             - added plausible value for V_patch, K_X
#       - [ ] X behaves weirdly...goes to 0 as expected but then jumps up, with now effect on individuals 
#       - [ ] ...


theme(:juno)


begin
    p.glb.t_max = 120.
    p.glb.dX_in = 2000. # 4,000 μgC ≈ 2 mg DW 
    p.glb.k_V = 0.1
    p.glb.V_patch = 0.5

    p.spc.Z = Truncated(Normal(1, 0.1), 0, Inf)
    p.spc.aux.tau_R = 2
    p.spc.aux.h_S = 1.
    p.spc.aux.a_max = Truncated(Normal(30, 12), 0, Inf)
    p.spc.K_X = 250.

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
        plot(:t, :X, marker = true) # why does X suddenly jump up, then increase linearly?
    )
end
using CSV
CSV.write("exploration/example_output_popcrash_glb.csv", sim.glb)
CSV.write("exploration/example_output_popcrash_spc.csv", sim.spc)

@transform! sim.spc :cause_of_death = map(x -> x.cause_of_death, sim.spc.aux)
@transform! sim.spc :Xi  = map(x -> x.Xi, sim.spc.aux)
@test sort(unique(sim.spc.cause_of_death)) == [0., 1., 2.]

using Chain

# FIXME: why does the population not recover at low N...?
#   - no growth whatsoever towards the end
#   - recording experienced food density Xi
#       - in fact, individuals have 0 food...?
#   - recording at every timestep

@chain sim.spc begin
    @subset :t .>= 70 :t .<= 75
    @with _ begin
        plot(
            plot(:t, :S, group = :id),
            plot(:t, :H, group = :id),
            plot(:t, :Xi, group = :id),
            leg = false
        )
    end
end

using ProgressMeter

dt = 1/24; num_steps = p.glb.t_max / dt
N = 1e3



prog = ProgressThresh(num_steps; desc=":")

for step in 1:100
    update!(prog, step, showvalues = ["N", N])
    sleep(0.1)
end

#Progress(10; dt=1) 


#@df sim.spc plot(
#    plot(:t, :S, group = :id),
#    plot(:t, [:is_embryo :is_juvenile :is_adult]), 
#    plot(:t, :X_emb), 
#    leg = false
#)




# the output for a single individual from the IBM is compared with the pure-ODE solution 
# this is done for a case with food limitation, since that is where there are most likely to be discrepancies
# (it is curently not possible to handle feedback with external food abundance in exactly the same way in both variants, and I am not sure it ever will be)
@testset "IBM vs ODE" begin

    debkiss = SimplifiedEnergyBudget() |> instantiate
    p = debkiss.parameters
    
    p.glb.dX_in = 1e2
    p.glb.k_V = 0.1
    p.glb.V_patch = 0.05
    p.glb.N0 = 1.
    p.glb.t_max = 56

    p.spc.Z = Dirac(1) # don't simulate individual variability
    p.spc.tau_R = Inf # don't simulate reproduction for now - just tracking repro buffer

    p.spc.H_p = 100.
    p.spc.a_max = Inf
    p.spc.K_X = 10

    @time global sim_ode = ETS.simulate_ODE(debkiss, saveat = 1);
    @time global sim_ibm = ETS.simulate_IBM(debkiss, dt = 1/240, showinfo = 14);

    plt = @df sim_ode plot(
        plot(:t, :X, xlabel = "t", ylabel = "X", leg = true, label = "ODE simulator"),
        plot(:t, :S, xlabel = "t", ylabel = "S", leg = true, label = "ODE simulator"),
        plot(:t, :H, xlabel = "t", ylabel = "H", leg = true, label = "ODE simulator"),
        plot(:t, :R, xlabel = "t", ylabel = "R", leg = true, label = "ODE simulator"),
        plot(:t, :f_X, xlabel = "t", ylabel = "f(X)", leg = true, label = "ODE simulator")
        )
    
    @df sim_ibm.glb plot!(:t, :X, label = "IBM simulator", subplot = 1)
    @df sim_ibm.spc plot!(:t, :S, label = "IBM simulator", subplot = 2)
    @df sim_ibm.spc plot!(:t, :H, label = "IBM simulator", subplot = 3)
    @df sim_ibm.spc plot!(:t, :R, label = "IBM simulator", subplot = 4)
    @df sim_ibm.spc plot!(:t, :f_X, label = "IBM simulator", subplot = 5)
    display(plt)
    
    sim_ibm.spc.t = round.(sim_ibm.spc.t)
    
    eval_df = leftjoin(sim_ibm.spc, sim_ode, on = :t, makeunique = true) |> 
    x -> select(x, [:S, :S_1, :H, :H_1, :R, :R_1]) |> 
    x -> EcotoxSystems.dropmissing(x)
    
    dS = eval_df.S .- eval_df.S_1
    dH = eval_df.H .- eval_df.H_1
    dR = eval_df.R .- eval_df.R_1
   
    @test sum(dS .> 2) == 0
    @test sum(dH .> 2) == 0
    @test sum(dR .> 2) == 0
end


@testset "Running IBM with default parameters" begin

    debkiss = SimplifiedEnergyBudget() |> instantiate
    p = debkiss.parameters

    p.glb.dX_in = 500_000 #100_000
    p.glb.k_V = 0.1
    p.glb.V_patch = 0.5
    p.glb.N0 = 10
    p.glb.t_max = 56

    p.spc.Z = Truncated(Normal(1, 0.1), 0, Inf)
    p.spc.tau_R = truncated(Normal(2., 0.2), 0, Inf)

    @time global sim_ibm = @replicates ETS.simulate_IBM(debkiss, saveat = 1, showinfo = 14) 1
    #VSCodeServer.@profview_allocs global sim_ibm = treplicates(x -> IBM_simulator(p, saveat = 1, showinfo = 14), p, 1)
    
    plt = plot(
        (@df sim_ibm.glb plot(:t, :N, group = :replicate, ylabel = "N")),
        (@df sim_ibm.spc groupedlineplot(:t, :S, :cohort, ylabel = "S")), 
        (@df sim_ibm.spc groupedlineplot(:t, :H, :cohort, ylabel = "H")), 
        leg = false
    )

    display(plt)

    # we know from experience that these values should be approximately reached for the given parameters
    
    @test 500 < maximum(sim_ibm.glb.N) < 1500
    @test 250 < maximum(sim_ibm.spc.S) < 800
    @test 50 < maximum(sim_ibm.spc.H) < 200

    # checking initial growth rate
    
    t0, N0 = sim_ibm.glb[1,[:t, :N]]
    t_peak, N_peak = sim_ibm.glb[argmax(sim_ibm.glb.N),[:t, :N]]

    r_int = (log(N_peak) - log(N0))/(t_peak-t0)

    @test 0.1 < r_int < 1.0
end

