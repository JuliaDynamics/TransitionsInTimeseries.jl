struct WindowViewer{T<:Real, V<:AbstractVector{<:T}}
    timeseries::V
    halfwidth::Int
    stride::Int
    strided_indices::AbstractVector{Int}
end

"""

    WindowViewer(
        timeseries::V,
        halfwidth::Int,
        stride::Int,
    ) where {T<:Real, V<:AbstractVector{<:T}}

Initialize an iterator that generates views over the given timeseries based on a window with a given `halfwidth`, incrementing the window views with the given `stride`. You can use this directly with `map`, such as `map(std, WindowViewer(x, ...))` would give you the moving-window-timeseries of the `std` of `x`.
"""
function WindowViewer(
    timeseries::V,
    halfwidth::Int,
    stride::Int,
) where {T<:Real, V<:AbstractVector{<:T}}

    n = length(timeseries)
    strided_indices = get_stride_indices(n, halfwidth, stride)
    return WindowViewer(timeseries, halfwidth, stride, strided_indices)
end

"""

    get_stride_indices(l::Int, halfwidth::Int, stride::Int)

Return a vector with strided indices based on windowing parameters.
"""
get_stride_indices(l::Int, halfwidth::Int, stride::Int) = 2*halfwidth+1:stride:l


# Define iterator for WindowViewer.
function Base.iterate(wv::WindowViewer, state::Int = 1)
    if state > length(wv.strided_indices)               # Stop condition: end of vector containing strided indices.
        return nothing
    else
        k = wv.strided_indices[state]
        i1, i2 = (k-2*wv.halfwidth, k)
        return (view(wv.timeseries, i1:i2), state + 1)  # Else: return a view of time series.
    end
end
Base.eltype(::Type{WindowViewer}) = AbstractVector{<:Real}
Base.length(wv::WindowViewer) = length(wv.strided_indices)