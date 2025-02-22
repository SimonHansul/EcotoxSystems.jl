# plotrecipes.jl
# Plots.jl recipess, useful for visualizing output dynamic ecotox models
# the lineplot and groupedlineplot recipes provide a functionality which is similar to seaborn.lineplot, but with differnt defaults and arguments.

@userplot RugPlot
@recipe function f(h::RugPlot)
    if length(h.args)>1
        error("Rugplot does not support multiple arguments.")
    end
    @series begin
        seriestype := :scatter
        markershape := :vline
        x, y = h.args[1], repeat([0], length(h.args[1]))
    end
end

"""
    lineplot(x, y, q; estimator)

Plot aggregated values of `y` over `x`. Similar to `seaborn.lineplot`.
By default, the arithmetic mean of `y` for every value of `x` is plotted, 
with 5th to 95th percentiles as ribbon. 
The plotted percentile range can be controlled via the positional argument `q`, which is a two-element Tuple.
The arithmetic mean can be replaced with another function via the `estimator` argument. 
"""
@userplot LinePlot
@recipe function f(h::LinePlot; estimator = mean)
    let x, y, q, mask
        
        if length(h.args)>=3
            x, y,  q = h.args
        else # if no quantile argument given at all => assume q=(0.05, 0.95)
            x, y, = h.args
            q = (.05, .95)
        end

        mask = isfinite.(x) .& isfinite.(y) 
        x = x[mask]
        y = y[mask]
        
        x_ag = []
        y_mn = []
        y_lo = []
        y_hi = []

        for (j,xi) in enumerate(unique(x))
            mask = x.==xi
            push!(x_ag, xi)
            push!(y_mn, estimator(y[mask]))
            push!(y_lo, quantile(y[mask], q[1]))
            push!(y_hi, quantile(y[mask], q[2]))
        end

        @series begin
            seriestype --> :path
            ribbon := @. (y_mn - y_lo, y_hi - y_mn)
            x_ag, y_mn
        end
    end
end

"""
    groupedlineplot(x, y, g, q; estimator)

Version of lineplot with additional grouping variable `g`.
"""
@userplot GroupedLinePlot
@recipe function f(h::GroupedLinePlot; estimator = mean)
    let x, y, group, q, mask
        # if no quantile argument given at all => assume q=(0.05, 0.95)
        if length(h.args)>=4
            x, y, group, q = h.args
        else
            x, y, group = h.args
            q = (.05, .95)
        end
        if typeof(group[1])<:Number
            mask = isfinite.(x) .& isfinite.(y) .& isfinite.(group) .& (ismissing.(x).==false) .& (ismissing.(y).==false) .& (ismissing.(group).==false)
        else
            mask = isfinite.(x) .& isfinite.(x) .& (ismissing.(x).==false) .& (ismissing.(y).==false)
        end
        x = x[mask]
        y = y[mask]
        group = group[mask]
        # aggregate values
        x_ag = []
        y_mn = []
        y_lo = []
        y_hi = []
        g_ag = []
        for (i,g) in enumerate(unique(group))
            for (j,xi) in enumerate(unique(x[group.==g]))
                mask = (group.==g) .& (x.==xi)
                push!(x_ag, xi)
                push!(y_mn, estimator(y[mask]))
                push!(y_lo, quantile(y[mask], q[1]))
                push!(y_hi, quantile(y[mask], q[2]))
                push!(g_ag, g)
            end
        end
        for (i,g) in enumerate(unique(g_ag))
            @series begin
                mask = g_ag .== g
                seriestype --> :path
                ribbon := (y_mn[mask] .- y_lo[mask], y_hi[mask] .- y_mn[mask])
                x_ag[mask], y_mn[mask]
            end
        end
    end
end