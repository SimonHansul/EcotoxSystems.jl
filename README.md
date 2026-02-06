# EcotoxSystems.jl : Simulating component-structured systems for ecotoxicology and chemical risk assessment


[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://simonhansul.github.io/EcotoxSystems.jl/dev/)
[![CI](https://github.com/SimonHansul/EcotoxSystems.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/SimonHansul/EcotoxSystems.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/SimonHansul/EcotoxSystems.jl/graph/badge.svg?token=BL1CFR86M6)](https://codecov.io/gh/SimonHansul/EcotoxSystems.jl)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

## Scope 

This package aims to 

- Provide a systematic approach for modular modelling of mechanistic models in ecotoxicology
- Facilitate the integration of organism-level ecotoxicological models into individual-based population models
- Provide a default model and parameters to serve as an entry point for teaching, research and development

Explicitly not within the scope of this package:

- Routines for model fitting. These are being developed in [EcotoxModelFitting.jl](https://github.com/simonhansul/ecotoxmodelfitting.jl.git) (work in progress)
- Providing pre-calibrated/validated models and parameter sets. Such things are better out-sourced into independent smaller modules.
- Providing an extensive set of pre-implemented models. Also better done in smaller separate modules, as done for the [AmphiDEB](https://github.com/simonhansul/amphideb) model.

## Implementation logic

The logic of `EcotoxSystems.jl` is based around the notion that a model can be divided into an equation-based portion and a rule-based portion. <br><br>
Furthermore, in the context of models developed for ecotoxicology and ecological risk assessment, 
we need to consider the organism-level as well as the environment (*individual* and *global* in the diagram below). <br>

The equation-based part is a system of ODEs. We make heavy use of [ComponentArrays.jl](https://github.com/jonniedie/ComponentArrays.jl) 
to define parameters and state variables for each component. 
<img src="implementation_logic.png" alt="Implementation logic">


## Installation 

While `EcotoxSystems.jl` is not registered, install directly from github:

```Julia
using Pkg; Pkg.add("https://github.com/simonhansul/ecotoxsystems.jl")
```


## Quickstart

### Using the pre-defined model

The following code executes the *default* model and parameters. <br>
This is a Dynamic Energy Budget Toxicokinetic-Toxicodynamic (DEB-TKTD) model based on the DEBkiss. <br> 
The derivatives are defined in `src/default_derivatives.jl`. 

```Julia
using EcotoxSystems
debkiss = SimplifiedEnergyBudget() |> instantiate
sim = ODE_simulator(debkiss)
```

The definition of `SimplifiedEnergyBudget` shows all the components that make up the fully functional simulation model:

```julia
Base.@kwdef mutable struct SimplifiedEnergyBudget <: AbstractEnergyBudget
    parameters::ComponentVector = debkiss_defaultparams

    # global component

    initialize_global_statevars::Union{Function,Nothing} = debkiss_global_statevars
    global_derivatives!::Union{Function,Nothing} = constant_nutrient_influx!
    global_rules!::Function = default_global_rules!
    
    # individual-level component
    
    initialize_individual_statevars::Function = debkiss_individual_statevars
    individual_derivatives!::Function = debkiss!
    individual_rules!::Function = default_individual_rules!
    generate_individual_params::Function = generate_individual_params

    # composed model

    initialize_all_statevars::Union{Function,Nothing} = nothing
    complete_derivatives!::Union{Function,Nothing} = nothing
end
```

The parameters `debkiss_defaultparams` are a component vector with global (`glb`) and species-level (`spc`) parameters. <br>

Species-level parameters are internally converted to individual-level parameters (`ind`) when passing them onto a simulator.
The function that converts `spc` to `ind` is also a keyword argument of `ODE_simulator`. <br>

### Individual variability

There are some convenience functions, e.g. to run replicated simulations:

```Julia
p = debkiss.parameters
p.spc.Z = Truncated(Normal(1,0.1), 0, Inf) # introduce individual variability through the mass-based zoom factor
sim = @replicates simulate(debkiss) 10 # simulate 10 times, each time sampling from Z
```
This is interesting if one of the parameters is subject to individual variability, 
as done for the mass-based zoom factor `Z` in the example above. <br>
Any species-level parameter can be defined as a probability distribution instead of a scalar value, 
and calling `simulate` on such a parameter set will cause it to take a random sample 
from the distribution to execute the simulation. <br>

### Running the IBM simulator

The following code simulates the default model as individual-based model (IBM):

```
import EcotoxSystems: simulate_IBM
sim = simulate_IBM(debkiss)
``` 

We can thus relatively easily switch between both modes. This can be useful for several reasons:

- We can develop and unit-test test the ODE part independently of the rule-based part
- During calibration, we often only need the ODE part.
- We can more easily re-combine different rules and ODE-systems, respectively, to study model behaviour, perform model selections, Bayesian model averaging, etc.


## Changelog 

### v0.1.0

Initial version.

### v0.1.1 

- Added the link_params! function, and link_params as argument of ODE_simulator
- Fixed a typo which caused tests to fail

### v0.1.2 

- Changed the default individual rules to use `S_rel_hist` and `W_S_rel_crit`: If structural mass decreases by a critical relative amount, constant hazard rate kicks in.
- Default callbacks use `save_positions = (false,false)` => time-points of callbacks will not be additionaly saved

### v0.1.3

- `init_global_statevars` and `init_individual_statevars` are now a keyword argument to `IBM_simulator` - necessary to implement modified models with additional states

### v0.1.4

- small bugfixes in IBM_simulator
- small improvements on utils

## v0.2.1

- optimized default model implementation by improving type stability

## v0.2.2

- fixed issue in default TKTD function

## v0.2.3

- `exposure()` can now take a Vector as argument

## v0.3.0

- updated depdendencies
- updated debkiss model implementation to struct-based interface
- updated tests to accomodate new struct-based itnerface 

