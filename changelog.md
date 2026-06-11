
# Changelog 

## v0.1.0

Initial version.

## v0.1.1 

- Added the link_params! function, and link_params as argument of ODE_simulator
- Fixed a typo which caused tests to fail

## v0.1.2 

- Changed the default individual rules to use `S_rel_hist` and `W_S_rel_crit`: If structural mass decreases by a critical relative amount, constant hazard rate kicks in.
- Default callbacks use `save_positions = (false,false)` => time-points of callbacks will not be additionaly saved

## v0.1.3

- `init_global_statevars` and `init_individual_statevars` are now a keyword argument to `IBM_simulator` - necessary to implement modified models with additional states

## v0.1.4

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

## v0.3.1

- change in dependency versions

## v0.3.2

- added missing method to fix `sol_to_df`

## v0.4.0

- Major update to internals and API! 
- updated DEBkiss implementation to omit function wrapping in struct (matching `EnergyBudgetModelZoo` approach)
- added type `FullDEBkiss <: SimplifiedEnergyBudget` for future integration into `EnergyBudgetModelZoo`
    - **⚠️ as of v0.5.0, model implementations should not be part of EcotoxSystems.jl anymore**
    - the long-term strategy is to separate models via SubModules, rather than structs
- for the IBM part, the use of structs is omitted entirely:
    - Individuals and models are just `NamedTuples`. Some exploratory analysis of the resulting performance showed that this is way more efficient because we can remove most of the type instability issues we previously had. 
    - The same strategy is applied to the Model representation.
    - The IBM implementation also omits wrapping of functions entirely, and instead uses closures to derive model-specific schedules based on a generic step and the component-wise functions.