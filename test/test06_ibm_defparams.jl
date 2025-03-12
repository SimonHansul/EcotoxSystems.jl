using Pkg; Pkg.activate("test")
#using Plots, StatsPlots

using Plots
default(leg = false)

using Distributions
using DataFramesMeta
using DataFrames
using Test

using Revise
@time using EcotoxSystems
@time import EcotoxSystems: params, ODE_simulator, IBM_simulator 
@time import EcotoxSystems: @replicates, DEBIndividual, treplicates

using Plots, StatsPlots

# the output for a single individual from the IBM is compared with the pure-ODE solution 
# this is done for a case with food limitation, since that is where there are most likely to be discrepancies
# (it is curently not possible to handle feedback with external food abundance in exactly the same way in both variants, and I am not sure it ever will be)
@testset "IBM vs ODE" begin

    p = params()

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

    @time global sim_ode = ODE_simulator(p, saveat = 1);
    @time global sim_ibm = IBM_simulator(p, dt = 1/240, showinfo = 14);

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


@testset "Running IBM" begin
    # FIXME: lots of memory allocs in the ODE part
    # maybe a problem in how I used broadcasting in the TKTD model, although removing TKTD part does not change so much
    # (cf https://discourse.julialang.org/t/can-i-avoid-allocations-when-broadcasting-over-slices/102501/2)
    # keyword "broadcast fusion" https://bkamins.github.io/julialang/2023/03/31/broadcast.html

    # we can play a little with the parameters here, but we are not looking for a "correct" prediction here!
    # the IBM_simulator just provides the infrastructure to combine the necessary components

    p = params()

    p.glb.dX_in = 50_000 #100_000
    p.glb.k_V = 0.1
    p.glb.V_patch = 0.5
    p.glb.N0 = 10
    p.glb.t_max = 56

    p.spc.Z = Truncated(Normal(1, 0.1), 0, Inf)
    p.spc.tau_R = truncated(Normal(2., 0.2), 0, Inf)

    @time global sim_ibm = @replicates IBM_simulator(p, saveat = 1, showinfo = 14) 1
    #VSCodeServer.@profview_allocs global sim_ibm = treplicates(x -> IBM_simulator(p, saveat = 1, showinfo = 14), p, 1)
    
    plt = plot(
        (@df sim_ibm.glb plot(:t, :N, group = :replicate, ylabel = "N")),
        (@df sim_ibm.spc groupedlineplot(:t, :S, :cohort, ylabel = "S")), 
        (@df sim_ibm.spc groupedlineplot(:t, :H, :cohort, ylabel = "H")), 
        leg = false
    )

    display(plt)

    # we know from experience that these values should be approximately reached for the given parameters
    
    #@test 500 < maximum(sim_ibm.glb.N) < 1500
    #@test 250 < maximum(sim_ibm.spc.S) < 800
    #@test 50 < maximum(sim_ibm.spc.H) < 200
end


#= 
Starting point: 561M allocations (15GB) for 365 days, 30 seconds comp time
For 56 days: 4-5 s, 81M allocs, 2.6GB

implemented "dS" function
        not much change
implemented "dI" function and inlined
        not much change; 4.3s and 2.36GB allocs; 11% GC 
implemented determine_S_max_hist and death_by_loss_of_structure
        allocations went down to 75M, comp time only slightly (2.3 GB)
        for 365 d? 
            same story, comp time is still at 28s, allocs only slightly down to 12.8 GB
changed life stage determination from CB to sig
        this made things slightly worse for the ODE benchmark
        small improvement for the IBM through
            3.3 - 4.2s, 2GV allocs
implemented dR function, inlined
        maybe minor improvement, allocs still at 1.9GB, comp time 3.3s

for 3655d, we are now down to 25s and 11.6 GB. making small improvements, but the bottleneck is still to be found

implemented dH function, inlined
    no measurable difference

implemented death_by_aging, inlined
    at most minor improvement

implemented dA, dM, dJ; all inlined
    comp time down to 3.2s, 1.6 GB allocs
    for 365d, we are still at ca 25s, 9.6 GB 

temporarily switching to profview. do we get a different result?
    no

trying so-called "direct indexing" (not sure if this is actually a thing). i.e. ind[:S] instead of ind.S etc.
for now, just for death_by_loss_of_structure
    this seems to actually do someting??
    comp time below 3s, just by changing the call to death_by_loss_of_structure

changed death_by_aging to direct indexing
    still 2.39s

same for calculating S_max_hist

refactoring reproduction rule, inlined
    no significant improvement
    for 365d, we are now at ca 22s, 9GB allocs

added a bit more direct indexing in repro 
        didn't help

replaced ind.y_j *= with ind.y_j = ind.y_j 
        no difference 

replaced p.ind with p[:ind]
        at most minor improvement, but also not many occurrences

changed m.dt to dt, assigning local variable
        no no no

=#

let p = params()

    p.glb.dX_in = 50_000 #100_000
    p.glb.k_V = 0.1
    p.glb.V_patch = 0.5
    p.glb.N0 = 10
    p.glb.t_max = 56

    p.spc.Z = Truncated(Normal(1, 0.1), 0, Inf)
    p.spc.tau_R = truncated(Normal(2., 0.2), 0, Inf)

    VSCodeServer.@profview_allocs IBM_simulator(p, saveat = 1, showinfo = 14)
    #@time IBM_simulator(p, saveat = 1, showinfo = 28)
end;


# comparison with netlogo
# running generic DEB-IBM with dt = 1/24 for 56d, population density up to ca 900: 33s
# EcotoxSystems.jl: roughly the same

let p = params()

    p.glb.dX_in = 60_000 #100_000
    p.glb.k_V = 0.1
    p.glb.V_patch = 0.5
    p.glb.N0 = 10
    p.glb.t_max = 56

    p.spc.Z = Truncated(Normal(1, 0.1), 0, Inf)
    p.spc.tau_R = truncated(Normal(2., 0.2), 0, Inf)

    #VSCodeServer.@profview_allocs IBM_simulator(p, saveat = 1, showinfo = 14)
    @time global sim = IBM_simulator(p, saveat = 1, showinfo = 28, dt = 1/240)
end;


# comparison with netlogo --> increasing population size by an order of magnitude
# netlogo crashes after day 19
# EcotoxSystems.jl: 
#  up to 10k individuals over 8 weeks with 1 time step = 36 s ==> 6 minutes
#  same with 1 time step = 1h ==> 82 seconds
 

let p = params()

    p.glb.dX_in = 600_000 #100_000
    p.glb.k_V = 0.1
    p.glb.V_patch = 0.5
    p.glb.N0 = 10
    p.glb.t_max = 56

    p.spc.Z = Truncated(Normal(1, 0.1), 0, Inf)
    p.spc.tau_R = truncated(Normal(2., 0.2), 0, Inf)

    #VSCodeServer.@profview_allocs IBM_simulator(p, saveat = 1, showinfo = 14)
    @time global sim = IBM_simulator(p, saveat = 1, showinfo = 1, dt = 1/24)
end;

@df sim.glb plot(:t, :N)


# simulating 10 years with approx 500 individuals peak, 300 individuals at later timepoints
#   simulation gets slow after 1 year...maybe because individual state vars are recorded?
#   trying again with record_individuals = false
#       simulator keeps going at constant speed now
#       comp time =  few minutes
#   increasing showinfo from 1 to 14
#       2 minutes to run 10 years. 90 GB allocations
#   increasing popsize to peak 1k, 10 years
#       5 minutes to run 10 years, 224 GB allocations

# comparison wit python, peak popsize approx 400, 10 years:  5 minutes
#       --> julia is currently ca. 2 times faster 


# at popsizes of 40-70k, netlogo takes 1 minute for 1 day 
# jula?

let p = params()

    p.glb.dX_in = 1_000_000 
    p.glb.k_V = 0.1
    p.glb.V_patch = 0.5
    p.glb.N0 = 10
    p.glb.t_max = 365 * 1

    p.spc.Z = Truncated(Normal(1, 0.1), 0, Inf)
    p.spc.tau_R = truncated(Normal(2., 0.2), 0, Inf)

    #VSCodeServer.@profview_allocs IBM_simulator(p, saveat = 1, showinfo = 14)
    @time global sim = IBM_simulator(p, saveat = 1, showinfo = 1, dt = 1/24, record_individuals = false)
end;
 