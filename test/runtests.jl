include("test00_setup.jl")

include("test01_defaults.jl") # simulates the default parameters
include("test02_food.jl") # performs a parameter sweep
include("test03_toxicity.jl") # simulates single stressors
include("test04_temperature.jl") # simulates temperature effects
include("test05_mixtures.jl") # simulates chemical mixtures
include("test06_ibm_defparams.jl") # simulates default individual-based model with defaultparams
include("test07_alternative_ind_params.jl") # adding individual parameters to the defaultparams

include("Aqua.jl")
