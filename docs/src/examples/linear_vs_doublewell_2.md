```julia
function plot_sol_doublewell(ttrend, Xtrend, Xres, sol, nx)

    for j in 1:nx
        lines!(axs[1][j], ttrend, Xtrend[j, :], label = L"trend $\,$")
        lines!(axs[2][j], ttrend, Xres[j, :], label = L"residual $\,$")

        lines!(
            axs[3][j],
            sol.tindctr,
            sol.reference_ti["ar1_whitenoise"][j, :],
            label = L"AR1 $\,$",
        )

        lines!(
            axs[4][j],
            sol.tindctr,
            sol.reference_ti["var"][j, :] ./ maximum(sol.reference_ti["var"][j, :]),
            label = L"normalized variance $\,$",
        )

        lines!(axs[5][j],
            sol.tindctr,
            sol.reference_ti["lfps"][j, :],
            label = L"LFPS $\,$",
        )

        [ylims!(axs[i][j], (-.1, 1.1)) for i in 3:nrows]

        lines!(
            axs[3][j],
            sol.tidtrend,
            Array(sol.significance["ar1_whitenoise"][j, :]), 
            label = L"AR1 significance $\,$",
        )

        lines!(
            axs[4][j],
            sol.tidtrend,
            Array(sol.significance["var"][j, :]), 
            label = L"variance significance $\,$",
        )

        lines!(
            axs[5][j],
            sol.tidtrend,
            Array(sol.significance["lfps"][j, :]), 
            label = L"LFPS significance $\,$",
        )

        lines!(
            axs[3][j],
            sol.tidtrend, Array(sol.positive_indicators[j, :]),
            label = L"positive indicators $\,$",
        )

        vlines!(
            axs[1][j],
            sol.predictor_time[j],
            label = L"warnings $\,$",
            color = :red,
            linestyle = :dash,
            linewidth = 1.0,
        )

        [axislegend(axs[i][2], position = :lt) for i in 1:nrows]
    end
end

plot_sol_doublewell(ttrend, Xtrend, Xres, sol)
```

Here the default settings were used to lump the percentile significance of the indicator trends (several continuous values) into a warning (a single binary value). Hence, a warning is only output when **all** the listed transition indicators display a percentile significance above 95%. If one wants a warning to be output when **at least two** indicators display a trend with a percentile significance above 99%, this can be done by running:

```julia
sol = indicate_transition(
    ttrend,
    Xres,
    p_window_indctr,
    p_window_idtrend,
    ns,
    [ar1_whitenoise, var, lfps],
    min_num_indicators = 3,
)
```

