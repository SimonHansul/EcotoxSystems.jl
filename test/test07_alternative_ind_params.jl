# if parameters contain additional components, 
# the function to generate individual parameters has to take this into account
# we test the functionality of providing an alternative generate_individual_params definition for both the ODE and IBM simulator


# default parameters with an additional component
p = ComponentVector(
    params(),
    foo = ComponentVector(x = 3.14)
)

