function roundint(x::Real)
    return Int( round( x ) )
end

function get_step(time::Real, dt::Real)
    return roundint(time / dt)
end

function check_std_endpoint(x::Vector{T}, xwin::Vector{T}, std_tol::Real) where {T<:Real}
    return max(xwin) > std_tol * std(x)
end