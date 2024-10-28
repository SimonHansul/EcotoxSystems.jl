
"""
Two-parameter log-logistic function. 
The expression `y = 1 / (1 + (x/p[1])^p[2])` is extended to 
`y = (1 / (1 + Complex(x/p[1])^p[2])).re`. 
This way we deal with domain errors which might occur if `x` or `p[1]` temporarily takes a negative value. 
Negative values should theoretically be impossible, but very small values of `x` (`<= 1e-20`) might occur during ODE solving. 
In this case, the returned real part of the expression evaluates to 1, which is in turn the expected behaviour.  
"""
@inline function LL2(x::Real, p::NTuple{2,Real})
    return (1 / (1 + Complex(x / p[1]) ^ p[2])).re
end

LL2(x, p1, p2) = LL2(x, (p1, p2))

"""
Cumulative hazard function of the log-logistic distribution. Mainly used for application in GUTS.
"""
@inline function LL2h(x::Real, p::NTuple{2,Real})
    -log((1 / (1 + Complex(x / p[1]) ^ p[2])).re)
end


"""
Two-parameter log-logistic function transformed to increasing function 
for application to PMoA maintenance costs.
"""
@inline function LL2M(x::Real, p::NTuple{2,Real})
    1 + (x/p[1])^p[2]
end


"""
Inverse of the two-parameter log-logistic function.
"""
function LL2inv(y::Real, p::NTuple{2,Real})
    return p[1] * (((1 / y) - 1)^(1 / p[2]))
end

"""
Inverse of the cumulative hazard two-parameter log-logistic function.
"""
function LL2hinv(y::Real, p::NTuple{2,Real})
    return p[1] * (((1 / exp(-y)) - 1)^(1 / p[2]))
end

"""
Two-parameter Weibull function.
"""
@inline function WB2(x::Real, p::NTuple{2,Real})
    return exp(-exp(p[2]*(log(x)-log(p[1]))))
end

@inline function WB2(x::Vector{Real}, p::NTuple{2,Real})
    return [WB2(xi, p) for xi in x]
end

"""
Bi-phasic log-logistic function, where each phase is described by a two-parameter log-logistic function. 
Requires additional "breakpoint" parameter.
Arguments:
    - `p1` = EC50_1
    - `p2` = beta_1
    - `p3` = EC50_2
    - `p4` = beta_2
    - `p5` = breakpoint = relative response at which second phase starts
"""
@inline function LLBP5(x::Real, p::NTuple{5,Real})
    y1 = p[5] / ( 1 +(x / p[1])^p[2]) # first-phase response
    y2 = p[5] / (1 + (x / p[3])^p[4]) # second-phase response
    y = y1 + y2 # total response
    return y
end

@inline function LLBP5(x::Vector{Real}, p::NTuple{5,Real})
    return [LLBP5(xi, p) for xi in x]
end


"""
    LLAS3(x::Real, p::NTuple{3,Real})

Asymmetric log-logistic function with additional slope parameter.
Arguments:
    - `p[1]` = EC50
    - `p[2]` = beta
    - `p[3]` = beta_2
"""
@inline function LLAS3(x::Real, p::NTuple{3,Real})
    return 1 / ((1 + ((x / p[1])^p[2]))^p[3])
end

function LLAS3(x::Vector{Real}, p::NTuple{3,Real})
    return [LLAS3(xi, p) for xi in x]
end

"""
    LL3(x::Real, p::NTuple{3,Real})

Three-parameter log-logistic function, where `p[3]` is the upper limit. 
"""
@inline function LL3(x::Real, p::NTuple{3,Real})
    return p[3] / (1 + (x / p[1])^p[2])
end


"""
Cedergreend-Ritz-Streibig model. \\
Parameters are \\

- `p[1]` = α = rate of hormetic increase
- `p[2]` = b = quasi-slope
- `p[3]` = c = lower limit
- `p[4]` = d = response in the control
- `p[5]` = e = inflection point
- `p[6]` = f = hormesis parameter
"""
CRS6(x::Real, p::NTuple{6,Real}) = @.(p[3] + ((p[4] - p[3] + p[6]*exp(-1 / (x^p[1]))) / (1 + exp(p[2] * (log(x) - log(p[5]))))))

"""
    CRS4(x::Real, p::NTuple{4,Real}) 

4-parameter CRS model. Lower limit and response in the control are fixed to 0 and 1, respectively.
"""
CRS4(x::Real, p::NTuple{4,Real}) = @.(((1 + p[4]*exp(-1 / (x^p[1]))) / (1 + exp(p[2] * (log(x) - log(p[3]))))))

"""
4-parameter CRS model transformed to u-shape, where the response in the control is 1, the lower limit is 0 and the maximum is 2.0."""
function CRS4U(x::Real, p::NTuple{4,Real})
    y = @. 2.0 .- (((1 + p[4] * exp(-1 / (x^p[1])) ) /(1 + exp(p[2] * (log(x) - log(p[3]))))))
    return max.(0, y)
end

"""
6-parameter U-shaped CRS model."""
function CRS6U(x::Real, p::NTuple{6,Real})
    y = @. p[4]-(p[3]+((p[4]-p[3]+p[6]*exp(-1/(x^p[1])))/(1+exp(p[2]*(log(x)-log(p[5]))))))
end

"""
Re-scaled 6-parameter U-shaped CRS model."""
function CRS6US(x::Real, p::NTuple{6,Real})
    y = @. 1 + (p[4]-(p[3]+((p[4]-p[3]+p[6]*exp(-1/(x^p[1])))/(1+exp(p[2]*(log(x)-log(p[5])))))))
    return max.(0, y)
end

"""
Re-scaled U-shaped CRS model with parameter C fixed to 0. \\
Parameters are \\
- alpha: rate of hormetic increase
- b: slope of the inclining part of the curve
- d: maximum stress 
- e: inflection point of the inclining part of the curve
- f: hormesis parameter"""
function CRS5US(x::Real, p::NTuple{5,Real})
    y = 1 + (p[3]-(((p[3]+p[5]*exp(-1/(x^p[1])))/(1+exp(p[2]*(log(x)-log(p[4])))))))
    return max.(0, y)
end

"""
"""
function CRS5US(x::Vector{Real}, p::NTuple{5,Real})
    y = @. 1 + (p[3]-(((p[3]+p[5]*exp(-1/(x^p[1])))/(1+exp(p[2]*(log(x)-log(p[4])))))))
    return max.(0, y)
end


"""
Continuous maximum function, using a sigmoid function to determine the maximum.
"""
contmax(a, b; beta = 100.) = 1. / (1. + exp(-beta*((a-b) - 0.))) * (y_right - y_left) + y_left

"""
Increasing NEC model with domain (1, Inf).
"""
@inline NEC2pos(x::Real, p::NTuple{2,Real}) = 1 + (p[2] * max(0, x - p[1]))

"""
Increasing NEC model with domain (0, Inf), mainly used in GUTS models. 
"""
@inline NEC2GUTS(x::Real, p::NTuple{2,Real}) = p[2] * max(0, x - p[1])



"""
Decreasing NEC model with domain (0, 1).
"""
@inline NEC2neg(x::Real, p::NTuple{2,Real}) = 1 / (1 + (p[2] * max(0, x - p[1])))

"""
Differentiable version of the NEC2pos model.
"""
@inline softNEC2pos(x::Real, p::NTuple{2,Real}) = 1 + (p[2] *(x - p[1]) * (0.5 + (1 / pi) * atan(1e6 * (x - p[1]))))
softNEC2pos(x, p1, p2) = softNEC2pos(x, (p1, p2))



"""
Differentiable version of the NEC2neg model.
"""
@inline softNEC2neg(x::Real, p::NTuple{2,Real}) = 1 / (1 + p[2] *(x - p[1]) * (0.5 + (1 / pi) * atan(1e6 * (x - p[1]))))
softNEC2neg(x, p1, p2) = softNEC2neg(x, (p1, p2))

"""
Differentiable version of the NEC2GUTS model.
"""
@inline softNEC2GUTS(x::Real, p::NTuple{2,Real}) = p[2] * (x - p[1]) * (0.5 + (1 / pi) * atan(1e6 * (x - p[1])))
softNEC2GUTS(x, p1, p2) = softNEC2GUTS(x, (p1, p2))
