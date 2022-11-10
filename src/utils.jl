function roundint(x::Real)
    return Int( round( x ) )
end

function get_step(time::Real, dt::Real)
    return roundint(time / dt)
end

function check_std_endpoint(x::Vector{T}, xwin::Vector{T}, std_tol::Real) where {T<:Real}
    return max(xwin) > std_tol * std(x)
end

###############################
# Handle dimensions
###############################

function catch_dim_error(X::A, N::Int, operation::String) where {A<:Union{Array{T}, CuArray{T}}} where {T<:Real}
    n = length(size(X))
    if n > N
        error("$operation computation not defined for Arrays of dimension n > $N.")
    end
end
struct Residual{T<:AbstractFloat, A<:Union{Array{T, 4}, CuArray{T, 4}}}
    X::A
    nx::Int
    ny::Int
    nz::Int
    nt::Int
end

"""

    structured_residual(X::A; verbose=true)

Build a struct containing important dimension information about residual.
"""
function structured_residual(X::A; verbose=true) where {A<:Union{Array{T,4}, CuArray{T,4}}} where {T<:AbstractFloat}
    if verbose
        println("-------------------------------------------------------------- \n")
        println("You are about to create a structured residual with dimensions: \n")
        println("nx: $(size(X, 1))")
        println("ny: $(size(X, 2))")
        println("nz: $(size(X, 3))")
        println("nt: $(size(X, 4)) \n")
        println("...and type: $A")
        println("-------------------------------------------------------------- \n")
    end
    return Residual(X, size(X)...)
end

struct FlatResidual{T<:AbstractFloat, A<:Union{Array{T, 2}, CuArray{T, 2}}}
    X::A
    nx::Int
    ny::Int
    nz::Int
    nt::Int
end

"""

    flatten_residual()

Flattens the residual along the time dimension.
"""
function flatten_residual(R::Residual)
    return FlatResidual( reshape(R.X, :, R.nt), R.nx, R.ny, R.nz, R.nt )
end

"""

    reshape_residual()

Reshapes flattened residual into original dimensions.
"""
function reshape_residual(R::FlatResidual)
    return reshape( R.X, R.nx, R.ny, R.nz, : )
end

"""

    get_placeholder_spacedims(X::A) where {A<:Union{Matrix{T}, CuArray{T, 2}}} where {T<:Real}

Extract placeholder for space dimensions of an Array with dims: nx x ny x nz x nt.
"""
function get_placeholder_spacedims(X::A) where {A<:Union{Array{T}, CuArray{T}}} where {T<:Real}
    n = length(size(X))
    space_dims = hcat([[:] for i in 1:n-1]...)
    return space_dims
end