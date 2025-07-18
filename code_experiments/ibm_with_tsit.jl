#=
# Incorporation more efficient solvers in the IBM schedule?
=# 

using Pkg; Pkg.activate("test")
Pkg.instantiate()

using Plots
default(leg = false)

using Distributions
using DataFramesMeta
using DataFrames
using Test

using Revise
@time using EcotoxSystems    
@time import EcotoxSystems: params, ODE_simulator, IBM_simulator 
@time import EcotoxSystems: @replicates, Individual, treplicates

using Plots, StatsPlots

using OrdinaryDiffEq

#=
the generic IBM step

- step global integrator 
- step global rules
=#

#=
Starting point with custom euler function: 
- 5 seconds for 56d with dX_in = 50k, N0 = 10
- 40 seconds for 365d same settings
=#

#=
trying to simulate just one
=#

p = deepcopy(EcotoxSystems.defaultparams)

p.glb.dX_in = 50_000 #100_000
p.glb.k_V = 0.1
p.glb.V_patch = 0.5
p.glb.N0 = 1
p.glb.t_max = 21

p.spc.Z = Truncated(Normal(1, 0.1), 0, Inf)
p.spc.tau_R = truncated(Normal(2., 0.2), 0, Inf)

# the code kind of works but tsit tends to get stuck on the second or third cohort...

@time global sim_ibm = @replicates IBM_simulator(p, saveat = 1, showinfo = 14, dt = 7.) 1;

@df sim_ibm.spc plot(:t, :S, marker = true, group = :id)
@df sim_ibm.glb plot(:t, :N)


@testset "Consistent time stepping" begin
    m = IndividualBasedModel(p)
    OrdinaryDiffEq.u_modified!(m.integrator, true)
    OrdinaryDiffEq.step!(m.integrator, m.dt, true)

    @test m.integrator.t == m.dt
end

m = IndividualBasedModel(p)
OrdinaryDiffEq.u_modified!(m.integrator, true)
OrdinaryDiffEq.step!(m.integrator, m.dt, true)
m.integrator.t 
m.integrator.u
m.integrator.du

m.integrator.p.glb.dX_in*(1/24)

a = m.individuals[1]

less(EcotoxSystems.model_step!)

m.integrator

m.integrator.du.glb .= 0.

EcotoxSystems.step_all_individuals!(m)

#=
## Defining a global ODE problem 

Here I am trying to adopt the approach explained in the Agents.jl docs. 
https://juliadynamics.github.io/Agents.jl/stable/examples/diffeq/#Coupling-DifferentialEquations.jl-to-Agents.jl

Although we are not using Agents.jl, it should be possible to use the same approach.
=#

m = EcotoxSystems.IndividualBasedModel(p)
u = m.integrator.u

u[:glb][:X]

global_ODE_problem = ODEProblem(
    EcotoxSystems.DEBODE_global!, 
    u, 
    (0.0, Inf), 
    p
    )
global_integrator = OrdinaryDiffEq.init(
    global_ODE_problem, 
    Tsit5();
    advance_to_tstop = true
)
OrdinaryDiffEq.u_modified!(global_integrator, true)
OrdinaryDiffEq.step!(global_integrator, m.dt, true)



global_integrator.t
global_integrator.u 


