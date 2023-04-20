using TransitionIndicators

abstract type TimeSeries end

struct EvenspacedTS <: TimeSeries
    nt::Int64
    t::AbstractVector
    spacedim::Int64
    x::AbstractArray
    dt::Real
end

struct UnevenspacedTS <: TimeSeries
    nt::Int64
    t::AbstractVector
    spacedim::Int64
    x::AbstractArray
end

function TimeSeries(t, x)
    nt = length(t)
    spacedim = ndims(x) - 1
    if isequispaced(t)
        dt = mean(diff(t))
        return EvenspacedTS(nt, t, spacedim, x, dt)
    else
        println(
            "\n" *
            " Unevenly spaced time series are work in progress. " * 
            "If you spot bugs, please open an issue on GitHub!\n " * 
            "Bare in mind that unevenly spaced time series come with limited functionalities! \n"
        )
        return UnevenspacedTS(nt, t, spacedim, x)
    end
end

struct WindowView
    width::Int
    stride::Int
    windows::Vector{UnitRange{Int}}
end

function window_view(x, wv::WindowView)
    return [view(x, window) for window in wv.windows]
end

function WindowView(ts::EvenspacedTS; width::Int = default_window_width(ts), stride::Int = default_window_stride(ts))
    si = width:stride:ts.nt
    windows = [k-width+1:k for k in si]
    return WindowView(width, stride, windows)
end
default_window_width(ts::EvenspacedTS) = ts.nt ÷ 100
default_window_stride(ts::EvenspacedTS) = 1

function WindowView(ts::UnevenspacedTS; twidth::Real = default_window_width(ts), stride::Int = default_window_stride(ts))
    ibegin = find_nearest(twidth, ts.t)
    if ibegin < 1
        ibegin = 1
    end
    windows = [findspan(t, twidth, ts.t) for t in ts.t[ibegin:stride:end]]
    return WindowView(0, stride, windows)
end
default_window_width(ts::UnevenspacedTS) = last(ts.t) / 100
default_window_stride(ts::UnevenspacedTS) = 1

function findspan(t, twidth, t_uneven)
    i1 = find_nearest(t-twidth, t_uneven)
    i2 = find_nearest(t, t_uneven)
    if i1 == i2
        error("Time series not dense enough before t=$t for the chosen window width")
    else
        return i1:i2
    end
end

function find_nearest(t, tvec::AbstractVector)
    return argmin((tvec .- t).^2)
end

t = 0:0.1:100
x = sin.(t)
ts = TimeSeries(t, x)
wv = WindowView(ts)
y = zeros(eltype(x), length(wv.windows))
@btime map!(var, y, window_view(x, wv))

tu = cumsum(0.2 .* rand(length(t)))
xu = sin.(t)
tsu = TimeSeries(tu, xu)
wvu = WindowView(tsu)
yu = zeros(eltype(x), length(wvu.windows))
@btime map!(var, yu, window_view(xu, wvu))