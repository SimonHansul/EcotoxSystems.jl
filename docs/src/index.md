# EcotoxSystems.jl

## Quickstart

To get acquainted with using `EcotoxSystems.jl`, you can use the default model and default parameters. <br>

Let's start with simulating individual life-history based on the default parameters:

```Julia
p = deepcopy(EcotoxSystems.debkiss_defaultparams)
p.glb.t_max = 21.
sim = EcotoxSystems.ODE_simulator(p)
first(sim, 5)
```

Here, we first create a copy of the default parameters. 
These contain two components:

- `glb` for *global* parameters. This includes simulation settings like the maximum simulated time-span, and forcings such as the food input rate `dX_in`.
- `spc` for *species-specific* parameters. For the default model, these are DEB-TKTD and some auxiliary parameters. 

The second line, 

```Julia 
p.glb.t_max = 21.
```

adjusts the simulated time-span. <br>
With 

```Julia 
sim = EcotoxSystems.ODE_simulator(p)
```
, we then simulate the model. <br>
`ODE_simulator` is the central function to execute models which are purely ODE-based. <br>
The output ``sim` is a `DataFrame` containing all state variables over time. <br>
Storing the solution in a `DataFrame` can be convenient, but, for simple models, consumes a considerable proportion of computation time. 
It is possible to skip the conversion and retrieve the `ODEsolution` instead by setting `returntype=EcoxSystems.odesol`. <br>

Instead of simulating the life-history of a single individual, we can also simulate population dynamics. 
To do so for the default model, we need to adjust the global parameters, 
set plausible value for the starvation parameters (which were irrelevant before), 
and call `IBM_simulator`:

```Julia
using Distributions, ProgressMeter

p = deepcopy(EcotoxSystems.debkiss_defaultparams)

# adjusting global parameters

p.glb.t_max = 365. # simulated timespan [d]
p.glb.V_patch = 0.5 # simulated volume [L]
p.glb.dX_in = 10_000 # provide more food 

# adjusting species-level parameters

p.spc.Z = Truncated(Normal(1, 0.1), 0, Inf) # induce individual variability 
p.spc.W_S_rel_crit = 0.66
p.spc.h_S = -log(0.5)

sims = @replicates EcotoxSystems.IBM_simulator(p, showinfo = Inf) 10 # run replicated simulations

```


To deviate from the default model, there are a number of components to tweak in order to fully specify the system:

- `init_global_statevars`: A function that initializes global state variables.
- `global_ode!`: ODE-portion of the dynamics of global states.
- `global_rules!`: Rule-based portion of the dynamics of global states. 
    
- `init_individual_statevars`: A function that initializes individual-level state variables. 
- `individual_ode!`: ODE-portion of the dynamics of individual-level states. 
- `individual_rules!`: Rule-based portion of the dynamics of individual-level states.
- `gen_ind_params`:  A function that translates species-specicif parameters `spc` into individual-level parameters `ind`. The default is `generate_individual_params` and is generic, apart from assuming that `spc` contains an entry for `Z` and `propagate_zoom`.


## API

```@autodocs
Modules = [EcotoxSystems]
```
