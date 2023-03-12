struct WindowViewer{T, V<:AbstractVector{T}, I<:AbstractVector{Int}}
    timeseries::V
    width::Int
    stride::Int
    strided_indices::I
end

"""

    WindowViewer(
        timeseries::V,
        width::Int,
        stride::Int,
    ) where {T<:Real, V<:AbstractVector{<:T}}

Initialize an iterator that generates views over the given timeseries based on a
window with a given `width`, incrementing the window views with the given
`stride`. You can use this directly with `map`, such as `map(std, WindowViewer(x, ...))`
would give you the moving-window-timeseries of the `std` of `x`.
"""
function WindowViewer(
        timeseries::V, width::Int, stride::Int,
    ) where {T<:Real, V<:AbstractVector{<:T}}
    n = length(timeseries)
    strided_indices = get_stride_indices(n, width, stride)
    I = typeof(strided_indices)
    return WindowViewer{T,V,I}(timeseries, width, stride, strided_indices)
end

"""

    get_stride_indices(l::Int, width::Int, stride::Int)

Return a vector with strided indices based on windowing parameters.
"""
get_stride_indices(l::Int, width::Int, stride::Int) = width+1:stride:l


# Define iterator for WindowViewer.
function Base.iterate(wv::WindowViewer, state::Int = 1)
    if state > length(wv.strided_indices)               # Stop condition: end of vector containing strided indices.
        return nothing
    else
        k = wv.strided_indices[state]
        i1, i2 = (k-wv.width, k)
        return (view(wv.timeseries, i1:i2), state + 1)  # Else: return a view of time series.
    end
end

function Base.eltype(::Type{<:WindowViewer{T,V}}) where {T,V}
    return SubArray{T, 1, V, Tuple{UnitRange{Int64}}, true}
end
Base.length(wv::WindowViewer) = length(wv.strided_indices)
Base.size(wv::WindowViewer) = (length(wv),)



"""
    windowmap(f::Function, x::AbstractVector, wv_width, wv_stride) → mapped_f

A shortcut for first generating a `wv = WindowViewer(x, wv_width, wv_stride)` and then
applying `mapped_f = map(f, wv)`. If `x` is accompanied by a time vector `t`,
you probably also want to call this function with `t` instead of `x` and with one of
`mean, first, last` as `f` to obtain a time vector for the `mapped_f` output.
"""
function windowmap(
    f::Function,
    x::AbstractVector,
    wv_width::Int,
    wv_stride::Int,
)
    wv = WindowViewer(x, wv_width, wv_stride)
    return map(f, wv)
end


# TODO: This function shouldn't exist. We really have to stop providing
# convenience functions to replace 2-loc calls.
"""

    get_windowmapview_length(x, wv_width, wv_stride)

Compute the length of the `WindowViewer` induced by `x`, `wv_width` and `wv_stride`.
"""
function get_windowmapview_length(x, wv_width, wv_stride)
    wv = WindowViewer(x, wv_width, wv_stride)
    return length(wv)
end
