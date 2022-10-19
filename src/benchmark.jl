function benchmark_functions_over_size(
    T::Type,
    nt::Int,
    nl_list::Vector{Int},
    function_list::Vector{Function},
)
    nf = length(function_list)
    nd = length(nl_list)

    speedup_ratio = zeros(T, nf, nd)
    for j in eachindex(nl_list)
        Xcpu = rand(T, nl_list[j], nt)
        Xgpu = CuArray(Xcpu)
        for i in eachindex(function_list)
            bmtime = benchmark_cpu_vs_gpu(Xcpu, Xgpu, function_list[i])
            speedup_ratio[i, j] = bmtime[1] / bmtime[2]
        end
    end
    return speedup_ratio
end

function benchmark_cpu_vs_gpu(Xcpu::Matrix{T}, Xgpu::CuArray{T,2}, func) where {T}
    bmcpu = @benchmark $func($Xcpu)
    bmgpu = @benchmark $func($Xgpu)
    return [StatsBase.mean(bmcpu.times), StatsBase.mean(bmgpu.times)]
end

function lineplot_benchmark(
    nl_list::Vector{Int},
    speedup_ratio,
    function_list::Vector{Function},
)
    labels = string.(function_list)
    fig = Figure()
    ax = Axis(
        fig[1, 1],
        title = "Ratio of CPU vs. GPU run-time",
        xscale = log10,
        yscale = log10,
    )
    for i in axes(speedup_ratio, 1)
        scatterlines!(ax, nl_list, speedup_ratio[i, :], label = labels[i])
    end
    hlines!(ax, [1.0f0], label = "Factor 1", color = :red)
    axislegend(ax, position = :lt)
    return fig, ax
end
