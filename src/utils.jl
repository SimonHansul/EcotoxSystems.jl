function extract_colnames(c::R, k::Symbol) where R <: Real
    return k
end

function extract_colnames(c::AbstractVector, k::Symbol)
    return [Symbol("$(k)_$(i)") for i in 1:length(c)]
end

function extract_colnames(u::ComponentVector)
    colnames = []
    for k in keys(u)
        push!(colnames, extract_colnames(u[k], k))
    end
    return vcat(colnames...)
end

extract_colnames(sol::ODESolution)::Vector{Symbol} = vcat([:t], extract_colnames(sol.u[1]))

function sol_to_mat(sol::O)::Matrix{Float64} where O <: ODESolution
    return hcat(sol.t, hcat(sol.u...)')
end

function sol_to_df(sol::O)::DataFrame where O <: ODESolution
    simout = DataFrame(
        sol_to_mat(sol), # ODE output converted to matrix
        extract_colnames(sol) # column names inferred from component array keys
    )
    return simout
end

robustmin(x) = length(x)>0 ? minimum(x) : Inf

ismin(x::AbstractVector) = x .== minimum(x)
