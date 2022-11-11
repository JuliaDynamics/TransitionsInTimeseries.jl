struct WindowViewer{T<:Real, V<:AbstractVector{<:T}}
    timeseries::V
    halfwidth::Int
    stride::Int
    bracketing::Function
    strided_indices::AbstractVector{Int}
end

"""

    WindowViewer(
        timeseries::V,
        halfwidth::Int,
        stride::Int,
        bracketing::Function,
    ) where {T<:Real, V<:AbstractVector{<:T}}

Returns struct of WindowViewer.
"""
function WindowViewer(
    timeseries::V,
    halfwidth::Int,
    stride::Int,
    bracketing::Function,
) where {T<:Real, V<:AbstractVector{<:T}}

    n = length(timeseries)
    strided_indices = bracketing(n, halfwidth, stride)
    return WindowViewer(timeseries, halfwidth, stride, bracketing, strided_indices)
end

"""

    left_bracketing(idx::Int, wv::WindowViewer)
    center_bracketing(idx::Int, wv::WindowViewer)
    right_bracketing(idx::Int, wv::WindowViewer)

Returns index span based on target index and parameters of the window viewer.
"""
left_bracketing(idx::Int, wv::WindowViewer) = (idx-2*wv.halfwidth, idx)
center_bracketing(idx::Int, wv::WindowViewer) = (idx - wv.halfwidth, idx + wv.halfwidth)
right_bracketing(idx::Int, wv::WindowViewer) = (idx, idx+2*wv.halfwidth)

"""

    left_bracketing(l::Int, halfwidth::Int, stride::Int)
    center_bracketing(l::Int, halfwidth::Int, stride::Int)
    right_bracketing(l::Int, halfwidth::Int, stride::Int)

Returns a vector with strided indices based on windowing parameters.
"""
left_bracketing(l::Int, halfwidth::Int, stride::Int) = 2*halfwidth+1:stride:l
center_bracketing(l::Int, halfwidth::Int, stride::Int) = halfwidth+1:stride:l-halfwidth
right_bracketing(l::Int, halfwidth::Int, stride::Int) = 1:stride:l-2*halfwidth

# Define iterator for WindowViewer.
function Base.iterate(wv::WindowViewer, state::Int = 1)
    if state > length(wv.strided_indices)               # Stop condition: end of vector containing strided indices.
        return nothing
    else
        k = wv.strided_indices[state]
        i1, i2 = wv.bracketing(k, wv)
        return (view(wv.timeseries, i1:i2), state+1)    # Else: return a view based on the chosen bracketing.
    end
end
Base.eltype(::Type{WindowViewer}) = AbstractVector{<:Real}
Base.length(wv::WindowViewer) = length(wv.strided_indices)