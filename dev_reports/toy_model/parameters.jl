using ComponentArrays; const CV = ComponentVector

toy_parameters = CV(
    glb = CV(N0 = 1., t_max = 21.), 
    spc = CV(
        L0 = 0.05,
        Lp = 0.5,
        Lm = 1., 
        rB = 0.2, 
        N_e = 10., 
        b_e = 2., 
        Rm = 5., 
        am = Truncated(Normal(60, 6), 0, Inf), 
        RperN = 1.
    )
)