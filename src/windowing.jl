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

Returns struct of WindowViewer with fields:
- `WindowViewer.timeseries`
- `WindowViewer.halfwidth`
- `WindowViewer.stride`
- `WindowViewer.strided_indices`: a vector of indices at which the window viewing should be applied.
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

    bracket(idx::Int, wv::WindowViewer)

Returns index span based on target index and parameters of the window viewer.
"""
bracket(idx::Int, wv::WindowViewer) = (idx-2*wv.halfwidth, idx)

"""

    get_stride_indices(l::Int, halfwidth::Int, stride::Int)

Returns a vector with strided indices based on windowing parameters.
"""
get_stride_indices(l::Int, halfwidth::Int, stride::Int) = 2*halfwidth+1:stride:l


# Define iterator for WindowViewer.
function Base.iterate(wv::WindowViewer, state::Int = 1)
    if state > length(wv.strided_indices)               # Stop condition: end of vector containing strided indices.
        return nothing
    else
        k = wv.strided_indices[state]
        i1, i2 = bracket(k, wv)
        return (view(wv.timeseries, i1:i2), state+1)    # Else: return a view of time series.
    end
end
Base.eltype(::Type{WindowViewer}) = AbstractVector{<:Real}
Base.length(wv::WindowViewer) = length(wv.strided_indices)