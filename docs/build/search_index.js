var documenterSearchIndex = {"docs":
[{"location":"#EcotoxSystems.jl","page":"EcotoxSystems.jl","title":"EcotoxSystems.jl","text":"","category":"section"},{"location":"","page":"EcotoxSystems.jl","title":"EcotoxSystems.jl","text":"Modules = [EcotoxSystems]","category":"page"},{"location":"#EcotoxSystems.defaultparams","page":"EcotoxSystems.jl","title":"EcotoxSystems.defaultparams","text":"Default parameter object\n\n\n\n\n\n","category":"constant"},{"location":"#EcotoxSystems.DEBODE!-NTuple{4, Any}","page":"EcotoxSystems.jl","title":"EcotoxSystems.DEBODE!","text":"DEB-ODE model with arbitrary number of stressors, assuming IA to compute combined effects. \n\n\n\n\n\n","category":"method"},{"location":"#EcotoxSystems.DEBODE_individual!-NTuple{4, Any}","page":"EcotoxSystems.jl","title":"EcotoxSystems.DEBODE_individual!","text":"Individual-level part of the DEB-ODE model with arbitrary number of stressors, assuming IA to compute combined effects.\n\n\n\n\n\n","category":"method"},{"location":"#EcotoxSystems.DEBkiss!-NTuple{4, Any}","page":"EcotoxSystems.jl","title":"EcotoxSystems.DEBkiss!","text":"DEBkiss!(du, u, p, t)::Nothing\n\nDynamics of DEBkiss model with maturity and explicit simulation of resource dynamics.\n\nThe density of structure is ignored, and instead S^(2/3) is applied for surface-area scaling.  This affects the dimension of dI_max, but has no effect on the model dynamics.\n\nIf model output is to be compared to length data, a statistical weight-lenght relationship  has to be applied to the model output.\n\n\n\n\n\n","category":"method"},{"location":"#EcotoxSystems.IBM_simulator-Tuple{ComponentArrays.ComponentVector}","page":"EcotoxSystems.jl","title":"EcotoxSystems.IBM_simulator","text":"IBM_simulator(\n    p::ComponentVector; \n    dt = 1/24, \n    saveat = 1,\n    showinfo::Number = Inf\n    )::DataFrame\n\nSimulate the individual-based version of the default model. \n\nimport EcotoxSystems as DEB\np = DEB.params()\nsim = DEB.IBMsimulator(p)\n\nargs\n\np: The parameter collection with defined global and species parameters.\n\nkwargs\n\ndt: Length of a timestep in the model (unit according to chosen unit of rate parameters)\nsaveat: Time interval at which to record output\nshowinfo: Time interval at which to print an update. Nothing will be printed if showinfo == Inf (the default).\n\n\n\n\n\n","category":"method"},{"location":"#EcotoxSystems.ODE_simulator-Tuple{ComponentArrays.ComponentVector}","page":"EcotoxSystems.jl","title":"EcotoxSystems.ODE_simulator","text":"simulator(\n    p::ComponentVector; \n    alg = Tsit5(),\n    saveat = 1,\n    reltol = 1e-6,\n    model = DEBODE!,\n    callbacks = DEBODE_callbacks,\n    returntype::ReturnType = dataframe,\n    kwargs...\n)::DataFrame\n\nRun the model as ODE system. \n\nExample: \n\nimport EcotoxSystems as DEB\np = DEB.params() # loads the default parameters\nsim = DEB.ODE_simulator(p)\n\n\n\n\n\n","category":"method"},{"location":"#EcotoxSystems.TKTD_mix_IA!-NTuple{4, Any}","page":"EcotoxSystems.jl","title":"EcotoxSystems.TKTD_mix_IA!","text":"function TKTD_mix_IA!(du, u, p, t)::Nothing\n\nMixture-TKTD for an arbitrary number of stressors, assuming Independent Action.\n\n\n\n\n\n","category":"method"},{"location":"#EcotoxSystems.calc_SL_max-Tuple{ComponentArrays.ComponentVector}","page":"EcotoxSystems.jl","title":"EcotoxSystems.calc_SL_max","text":"calc_SL_max(spc::ComponentVector)::Float64\n\nCalculate maximum structural length slmax [m^(1/3)]\n\n\n\n\n\n","category":"method"},{"location":"#EcotoxSystems.calc_S_max-Tuple{ComponentArrays.ComponentVector}","page":"EcotoxSystems.jl","title":"EcotoxSystems.calc_S_max","text":"calc_S_max(spc::ComponentVector)::Float64\n\nCalculate maximum structural mass smax [m]\n\n\n\n\n\n","category":"method"},{"location":"#EcotoxSystems.clipneg-Tuple{Real}","page":"EcotoxSystems.jl","title":"EcotoxSystems.clipneg","text":"Clip negative values at 0 as a continuous function, using sig.\n\n\n\n\n\n","category":"method"},{"location":"#EcotoxSystems.default_individual_rules!-Tuple{EcotoxSystems.AbstractDEBIndividual, EcotoxSystems.AbstractDEBIBM}","page":"EcotoxSystems.jl","title":"EcotoxSystems.default_individual_rules!","text":"default_individual_rules(a::AbstractDEBIndividual, m::AbstractDEBIBM)::Nothing\n\nDefines the default rule-based portion for DEBIndividuals. <br>\n\nThe event functions which are used as callbacks during ODE solving are here re-used to apply rules for life stage transitions. A crude rule for starvation mortality is implemented, applying a dependency between      starvation mortality and the scaled functional response.\n\nReproduction is assumed to occur in fixed time intervals, according to spc.tau_R.\n\n\n\n\n\n","category":"method"},{"location":"#EcotoxSystems.exposure-Union{Tuple{R}, Tuple{Function, ComponentArrays.ComponentVector, Matrix{R}}} where R<:Real","page":"EcotoxSystems.jl","title":"EcotoxSystems.exposure","text":"exposure(\n    simulator::Function, \n    p::ComponentVector, \n    C_Wmat::Matrix{R}\n) where R <: Real\n\nSimulate exposure to an arbitrary number of stressors over a Matrix of constant exposure concentrations C_Wmat. \n\nThe exposure matrix columns are stressors, the rows are treatments.  That means, to simulate a single-stressor experiment, do \n\nC_Wmat = [0.; 1.; 2;]\n\n, creating a n x 1 matrix with exposure concentrations 0, 1 and 2. \n\nIn contrast, a single treatment with multiple stressors would be defined as \n\nC_Wmat = [0 1 2;]\n\n, creating a 1 x n matrix. Here, we would have three stressors with the simultaneous exposure concentrations 0,1,2. <br>\n\nThus, defining four treatments for two stressors looks like this:\n\nC_Wmat = [0.0 0.0; 0.0 0.5; 1.0 0.0; 0.5 1.0]\n\nThis exposure matrix corresponds to a ray design with constant exposure ratios.\n\n\n\n\n\n","category":"method"},{"location":"#EcotoxSystems.generate_individual_params-Tuple{ComponentArrays.ComponentVector}","page":"EcotoxSystems.jl","title":"EcotoxSystems.generate_individual_params","text":"Generate individual-specific parameter set from species-specific parameter set.  If a parameter entry is a distribution, a random sample is taken. \n\nThis also works for Vectors of distributions.\n\n\n\n\n\n","category":"method"},{"location":"#EcotoxSystems.individual_step!-Tuple{EcotoxSystems.AbstractDEBIndividual, EcotoxSystems.AbstractDEBIBM}","page":"EcotoxSystems.jl","title":"EcotoxSystems.individual_step!","text":"individual_step!(a::Agent, m::Model)::Nothing\n\nThe individual step follows a generic pattern:\n\nFirst the ODE-portion of the model is executed, and the corresponding state variables are updated using the Euler scheme. \n\nThen the rule-based portion of the model is executed. These are all the functions which cannot / should not be expressed as part of an ODE. At the minimum, this will include life stage transitions, reproduction and death of individuals. \n\nFor a spatially explicit model, movement should also most likely be part of the rule-based portion,  as well as functions which require direct information exchange between individuals.\n\n\n\n\n\n","category":"method"},{"location":"#EcotoxSystems.initialize_statevars-Tuple{ComponentArrays.ComponentVector}","page":"EcotoxSystems.jl","title":"EcotoxSystems.initialize_statevars","text":"initialize_statevars(p::::ComponentVector, pindt::ComponentVector{Float64})::ComponentArray\n\nFor initialization of ODE simulator, initialize the component vector of state variables, u, based on common oaraeter collection p.\n\n\n\n\n\n","category":"method"},{"location":"#EcotoxSystems.k_J!-Tuple{ComponentArrays.ComponentVector}","page":"EcotoxSystems.jl","title":"EcotoxSystems.k_J!","text":"k_J!(spc::ComponentVector)::Nothing\n\nSet the maturity maintenance rate constant,  assuming that the cumulative investment into maturity maintenance  equals the cumulative investment into somatic maintenance (cf. DEBkiss book by Tjalling Jager).\n\n\n\n\n\n","category":"method"},{"location":"#EcotoxSystems.relative_response-Union{Tuple{D}, Tuple{D, Vector{Symbol}, Symbol}} where D<:DataFrames.AbstractDataFrame","page":"EcotoxSystems.jl","title":"EcotoxSystems.relative_response","text":"relative_response(\n    sim::D, \n    response_vars::Vector{Symbol},\n    treatment_var::Symbol; \n    groupby_vars::Vector{Symbol} = Symbol[],\n    identify_control = minimum\n    ) where D <: AbstractDataFrame\n\nCalculate relative responses. \n\nargs:\n\nsim::AbstractDataFrame: results\nresponse_vars::Vector{Symbol}: response variables for which to calculate the relative responses\ntreatment_var::Symbol: Column indicating the treatment. Column values can be numerical or categorical, but identify_control kwarg has to be specified in the latter case\n\nkwargs:\n\ngroupby_vars::Vector{Symbol}: relative response will be conditioned on these variables (e.g. time, separate experiments...). Empty by default.\nidentify_control: function used to identify reference values from treatment_var. By default, this is minimum() (assuming that column values in treatment_var are numerical).\n\n\n\nExample\n\nusing MechanisticEffectModels.EcotoxSystems, MechanisticEffectModels.Utils\n\nsimfunct(x) = @replicates EcotoxSystems.simulator(x) 10\n\nsim = exposure(simfunct, Params(), [0., 100., 200.]) |>\nx -> relative_response(x, [:S, :H, :R]) # -> data frame will contain columns y_S, y_H, y_R for control-normalized values\n\n\n\n\n\n","category":"method"},{"location":"#EcotoxSystems.replicates-Tuple{Function, ComponentArrays.ComponentVector, Int64}","page":"EcotoxSystems.jl","title":"EcotoxSystems.replicates","text":"replicates(simulator::Function, defaultparams::ComponentVector, nreps::Int64; kwargs...)\n\nPerform replicated runs of simulator with parameters defaultparams (simulator(defaultparams) has to be a valid function call).  Analogous to @replicates, but a bit more flexible.\n\n\n\n\n\n","category":"method"},{"location":"#EcotoxSystems.sig-NTuple{4, Real}","page":"EcotoxSystems.jl","title":"EcotoxSystems.sig","text":"sig(\n    x::Real, \n    x_thr::Real,\n    y_left::Real, \n    y_right::Real; \n    beta::Real = 30\n    )::Real\n\nSigmoid switch function.  Used to replace simple if-statements with a continuous function in ODE models. \n\ny_left and y_right are the function values left and right of the threshold x_thr.\n\n\n\n\n\n","category":"method"},{"location":"#EcotoxSystems.softNEC2GUTS-Tuple{Real, Tuple{Real, Real}}","page":"EcotoxSystems.jl","title":"EcotoxSystems.softNEC2GUTS","text":"Transformation of softNECpos for increasing relationships with domain (0,Inf).\n\n\n\n\n\n","category":"method"},{"location":"#EcotoxSystems.softNEC2neg-Tuple{Real, Tuple{Real, Real}}","page":"EcotoxSystems.jl","title":"EcotoxSystems.softNEC2neg","text":"Transformation of softNEC2pos for decreasing relationships with domain (0,1).\n\n\n\n\n\n","category":"method"},{"location":"#EcotoxSystems.softNEC2pos-Tuple{Real, Tuple{Real, Real}}","page":"EcotoxSystems.jl","title":"EcotoxSystems.softNEC2pos","text":"Linear positive relationship with soft threshold (NEC) and domain (1,Inf).\n\n\n\n\n\n","category":"method"},{"location":"#EcotoxSystems.texposure-Tuple{Function, ComponentArrays.ComponentVector, Vector{Float64}}","page":"EcotoxSystems.jl","title":"EcotoxSystems.texposure","text":"Threaded version of exposure().\n\n\n\n\n\n","category":"method"},{"location":"#EcotoxSystems.treplicates-Tuple{Function, ComponentArrays.ComponentVector, Int64}","page":"EcotoxSystems.jl","title":"EcotoxSystems.treplicates","text":"treplicates(\n    simulator::Function, \n    defaultparams::ComponentVector, \n    nreps::Int64; \n    kwargs...)\n\nMulti-threaded version of replicates. \n\nOnly useful if Julia has been started with multiple threads. \n\nTo check the number of threads, run  using Base.Threads; Threads.nthreads().\n\nIn VSCode, you can use the entry \"julia.NumThreads\" in settings.json to set the default number of threads  (searching for \"julia threads\" in the preferences will lead you there). \n\nCheck the Multi-threading documentation  for more information.\n\n\n\n\n\n","category":"method"},{"location":"#EcotoxSystems.@replicates-Tuple{Expr, Int64}","page":"EcotoxSystems.jl","title":"EcotoxSystems.@replicates","text":"@replicates(simcall::Expr, nreps::Int64)\n\nPerform replicated runs of simcall, where simcall is a call to a simulator function. \n\nExample:\n\n    spc = SpeciesParams(Z = Truncated(Normal(1, 0.1), 0, Inf)) # initialize default parameters with variable zoom factor\n    sim = @replicates MechanisticEffectModels.simulator(Params(spc = spc))) 10 # execute replicated runs to simulator\n\nIn this case, sim will contain the output of 10 replicated simulations. For each replicate, the zoom factor is sampled from a truncated Normal distribution.  sim contains an additional column replicate.\n\n\n\n\n\n","category":"macro"}]
}