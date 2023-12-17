module TITVisualizations

using TransitionsInTimeseries, Makie

function TransitionsInTimeseries.plot_indicator_changes(res::SlidingWindowResults,
        colors =  ["#7143E0", "#0A9A84", "#191E44", "#AF9327", "#701B80", "#2E6137",],
        skip = Int[],
    )
    config = res.config
    fig = Figure()
    axts = Axis(fig[1,1]; ylabel = "input")
    if isnothing(config.indicators)
        axcha = Axis(fig[2,1]; ylabel = "changes")
        linkxaxes!(axts, axcha)
    else
        axind = Axis(fig[2,1]; ylabel = "indicators")
        axcha = Axis(fig[3,1]; ylabel = "changes", xlabel = "time")
        linkxaxes!(axts, axind, axcha)
        hidexdecorations!(axind; grid = false)
    end
    hidexdecorations!(axts; grid = false)

    lines!(axts, res.t, res.x; color = "black", linewidth = 3)
    for (i, cha) in enumerate(config.change_metrics)
        i ∈ skip && continue
        if !isnothing(config.indicators)
            lines!(axind, res.t_indicator, res.x_indicator[:, i];
                color = colors[i], linewidth = 3, label = string(nameof(config.indicators[i]))
            )
        end
        lines!(axcha, res.t_change, res.x_change[:, i];
            color = colors[i], linewidth = 3, label = string(nameof(cha))
        )
    end

    !isnothing(config.indicators) && axislegend(axind)
    axislegend(axcha)
    return fig
end

function TransitionsInTimeseries.plot_significance!(fig::Figure, res::SlidingWindowResults, signif::SurrogatesSignificance;
        colors = ["#7143E0", "#0A9A84", "#191E44", "#AF9327", "#701B80", "#2E6137",],
        skip = Int[],
        nsurro = 20,
    )
    config = res.config
    for (i, cha) in enumerate(config.change_metrics)
        i ∈ skip && continue

        for _ in 1:nsurro
            s = TimeseriesSurrogates.surrogate(res.x, signif.surrogate)
            if isnothing(config.indicators)
                p = s
                chaj = 2
            else
                p = windowmap(config.indicators[i], s; width = config.width_ind, stride = config.stride_ind)
                lines!(fig[2, 1], res.t_indicator, q; color = (colors[i], 2/nsurro), linewidth = 1)
                chaj = 3
            end
            q = windowmap(cha, p; width = config.width_cha, stride = config.stride_cha)
            lines!(fig[chaj, 1], res.t_change, q; color = (colors[i], 2/nsurro), linewidth = 1)
        end
    end

    return fig
end


end