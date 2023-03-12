"""
    isequispaced(t; kwargs...) → true/false

Return `true` if time vector `t` is evenly spaced.
For `AbstractRange` this is always true, while for `AbstractVector` successive
differences in `t` are compared with `isapprox(; kwargs...)` and if all are approximately
the same then `true` is returned.
"""
function isequispaced(t::AbstractVector; kwargs...)
    # TODO: assert 1 based indexing
    s = t[2] - t[1]
    for i in 3:length(t)
        sn = t[i] - t[i-1]
        isapprox(sn, s; kwargs...) || return false
    end
    return true
end
isequispaced(t::AbstractRange; reltol = 1e-6) = true

equispaced_step(t::AbstractRange) = step(t)
function equispaced_step(t::AbstractVector)
    s = zero(eltype(t))
    for i in 2:length(t)
        s += t[i] - t[i-1]
    end
    return s/(length(t)-1)
end
