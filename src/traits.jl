"""
    calc_SL_max(spc::ComponentVector)::Float64
Calculate maximum structural length slmax [m^(1/3)]
"""
function calc_SL_max(spc::ComponentVector)::Float64
    return ((spc.kappa_0 * spc.Idot_max_rel_0 * spc.eta_IA_0) / spc.k_M_0)
end

"""
    calc_S_max(spc::ComponentVector)::Float64
Calculate maximum structural mass smax [m]
"""
function calc_S_max(spc::ComponentVector)::Float64
    return calc_SL_max(spc)^3
end

"""
    k_J!(spc::ComponentVector)::Nothing

Set the maturity maintenance rate constant, 
assuming that the cumulative investment into maturity maintenance 
equals the cumulative investment into somatic maintenance (cf. DEBkiss book by Tjalling Jager).
"""
function k_J!(spc::ComponentVector)::Nothing
    spc.k_J_0 = ((1 - spc.kappa_0) / spc.kappa_0) * spc.k_M_0

    return nothing
end

age_at_birth(sim::AbstractDataFrame) = sim[ismin(abs.(sim.X_emb)),:t] |> robustmin