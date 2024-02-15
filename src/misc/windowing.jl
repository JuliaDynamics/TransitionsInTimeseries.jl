struct WindowViewer{T, V<:AbstractVector{T}}
    timeseries::V
    width::Int
    stride::Int
    strided_indices::StepRange{Int, Int}
end

"""
    WindowViewer(x; width, stride)

Initialize an iterator that generates views over the given timeseries `x` based on a
window with a given `width`, incrementing the window views with the given
`stride`. You can use this directly with `map`, such as `map(std, WindowViewer(x, ...))`
would give you the moving-window-timeseries of the `std` of `x`.

If not given, the keywords `width, stride` are respectively taken as
`default_window_width(x)` and `1`.
"""
function WindowViewer(
        x::AbstractVector;
        width::Int = default_window_width(x), stride::Int = DEFAULT_WINDOW_STRIDE,
    )
    n = length(x)
    si = width:stride:n
    return WindowViewer{eltype(x),typeof(x)}(x, width, stride, si)
end

# Define iterator for WindowViewer.
function Base.iterate(wv::WindowViewer, state::Int = 1)
    if state > length(wv.strided_indices) # Stop condition: end of vector containing strided indices.
        return nothing
    else # return a view of time series.
        k = wv.strided_indices[state]
        i1, i2 = (k-wv.width+1, k)
        return (view(wv.timeseries, i1:i2), state + 1)
    end
end

function Base.eltype(::Type{<:WindowViewer{T,V}}) where {T,V}
    return SubArray{T, 1, V, Tuple{UnitRange{Int64}}, true}
end
Base.length(wv::WindowViewer) = length(wv.strided_indices)
Base.size(wv::WindowViewer) = (length(wv),)

"""
    windowmap(f::Function, x::AbstractVector; kwargs...) → mapped_f

A shortcut for first generating a `wv = WindowViewer(x; kwargs...)` and then
applying `mapped_f = map(f, wv)`. If `x` is accompanied by a time vector `t`,
you probably also want to call this function with `t` instead of `x` and with one of
`mean, midpoint, midvalue` as `f` to obtain a time vector for the `mapped_f` output.
"""
function windowmap(f::Function, x::AbstractVector; kwargs...)
    wv = WindowViewer(x; kwargs...)
    return map(f, wv)
end

"""
    windowmap!(f::Function, out, x::AbstractVector; kwargs...)

Same as [`windowmap`](@ref), but writes the output in-place in `out`.
"""
function windowmap!(f::F, out::AbstractVector, x::AbstractVector; kwargs...) where {F}
    wv = WindowViewer(x; kwargs...)
    X = eltype(out)
    if length(out) != length(wv)
        throw(ArgumentError(
            "Allocated output doesn't match size of window viewer, "*
            "got $(length(out)) instead of $(length(wv))."
        ))
    end
    @inbounds for (i, v) in enumerate(wv)
        val::X = f(v)
        out[i] = val
    end
    return out
end

"""
    midpoint(x)

Return `x[midindex]` with `midindex = (firstindex(x) + lastindex(x))÷2`.
Typically useful in [`windowmap`](@ref) with a time vector.
"""
midpoint(x) = x[(firstindex(x) + lastindex(x))÷2]
midpoint(x::Vector) = x[length(x)÷2] # faster

"""
    midvalue(x)

Return `(first(x) + last(x))/2`.
Typically useful in [`windowmap`](@ref) with a time vector.
"""
midvalue(x) = (first(x) + last(x))/2
