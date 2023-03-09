"""

    is_equispaced(t; reltol = 1e-6)

Return `true` if time vector `t` is evenly spaced w.r.t. a relative tolerance `reltol`
and `false` otherwise. Return mean time step as second output.
"""
function is_equispaced(t::AbstractVector{T}; reltol = 1e-6) where {T<:Real}
    mean_dt, stddev_dt = mean_and_std(diff(t))
    return stddev_dt < reltol * maximum(t), mean_dt
end