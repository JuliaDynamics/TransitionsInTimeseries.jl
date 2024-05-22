"""
    isequispaced(t; kwargs...) → true/false

Return `true` if time vector `t` is evenly spaced.
For `AbstractRange` this is always true, while for `AbstractVector` successive
differences in `t` are compared with `isapprox(; kwargs...)` and if all are approximately
the same then `true` is returned.
"""
function isequispaced(t::AbstractVector; kwargs...)
    fi = firstindex(t)
    s = t[fi+1] - t[fi]
    for i in 2:length(t)-1
        sn = t[fi+i] - t[fi+i-1]
        isapprox(sn, s; kwargs...) || return false
    end
    return true
end
isequispaced(t::AbstractRange; reltol = 1e-6) = true

equispaced_step(t::AbstractRange) = step(t)
function equispaced_step(t::AbstractVector)
    fi = firstindex(t)
    s = zero(eltype(t))
    for i in 1:length(t)-1
        s += t[fi+i] - t[fi+i-1]
    end
    return s/(length(t)-1)
end
