# if parameters contain additional components, 
# the function to generate individual parameters has to take this into account
# we test the functionality of providing an alternative generate_individual_params definition for both the ODE and IBM simulator


# default parameters with an additional component
@testset "Adding a component to parameters: ODE simulator" begin 
    p = ComponentVector(
    EcotoxSystems.params(),
    foo = ComponentVector(x = 0.5)
    )

    # constructing an ODE system that uses the additional component
    # not important for this test what it actually does, 
    # only that we can access p.foo.x from within the individual step

    function foo_individual!(du, u, p, t)
        EcotoxSystems.default_TKTD!(du, u, p, t)
        EcotoxSystems.default_physiology!(du, u, p, t)
        du.ind.S *= p.foo.x
    end

    function foo!(du,u,p,t)
        EcotoxSystems.default_global_ODE!(du, u, p, t)
        foo_individual!(du, u, p, t)
    end

    gen_ind_params(p) = EcotoxSystems.generate_individual_params(p;foo = p.foo)

    # FIXME: why is this much slower than running the default model?
    @time sim = EcotoxSystems.ODE_simulator(
        p;
        model = foo!,
        gen_ind_params = gen_ind_params
    )

    # if this code runs without throwing an error, 
    # it means we can access p.foo from within foo!, which is all that this test needs to show
    @test true
end


@testset "Adding a component to parameters: IBM simulator" begin 
    
    p = ComponentVector(
        EcotoxSystems.params(),
        foo = ComponentVector(x = 0.5)
    )

    p.spc.tau_R = 2.

    # constructing an ODE system that uses the additional component
    # not important for this test what it actually does, 
    # only that we can access p.foo.x from within the individual step

    function foo_individual!(du, u, p, t)
        EcotoxSystems.default_TKTD!(du, u, p, t)
        EcotoxSystems.default_physiology!(du, u, p, t)
        du.ind.S *= p.foo.x
    end

    gen_ind_params(p) = EcotoxSystems.generate_individual_params(p; foo = p.foo)

    @time sim = EcotoxSystems.IBM_simulator(
        p;
        
        individual_ode! = foo_individual!,
        gen_ind_params = gen_ind_params
    )

    # if this code runs without throwing an error, 
    # it means we can access p.foo from within foo!, which is all that this test needs to show
    @test true
end
