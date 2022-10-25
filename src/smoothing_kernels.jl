#=
Smoothing can be obtained by convolution of kernel with signal.
Here we define the kernels on a normed window [-1, 1].
We then scale in x to match the actual window.
We finally scale in y to have a kernel summing to 1 over the window.
=#

"""
    scaled_kernel(T::Type, p::WindowingParams, kernel::Function)

Computes a smoothing discrete kernel of type `T` and windowing defined by `p`.
Functions available for `kernel`:
    - uniform_kernel
    - triangular_kernel
    - parabolic_kernel
    - biweight_kernel
    - triweight_kernel
    - tricube_kernel
    - gaussian_kernel
    - cosine_kernel
    - logistic_kernel
    - sigmoid_kernel

For now only centered window implemented.
For more details about the kernels, refer to https://en.wikipedia.org/wiki/Kernel_(statistics).
"""
# TODO implement for different window types.
function scaled_kernel(
    T::Type,
    p::WindowingParams,
    kernel::Function,
)
    return scale_y_kernel( kernel( scale_x_kernel(T, p) ) )
end

# Scale such that x is mapped on [-1, 1].
scale_x_kernel(T::Type, p::WindowingParams) = T.(-p.Nwndw:p.Nwndw) ./ p.Nwndw

# Scale such that sum = 1.
scale_y_kernel(u::Vector{T}) where {T} = u ./ sum(u)

# Common smoothing kernels as found in https://en.wikipedia.org/wiki/Kernel_(statistics)
uniform_kernel(u::Vector{T}) where {T} = fill( T(1), length(u) )
triangular_kernel(u::Vector{T}) where {T} = 1 .- abs.( u )
parabolic_kernel(u::Vector{T}) where {T} = T(3/4) .* ( 1 .- u.^2 )
biweight_kernel(u::Vector{T}) where {T} = T(15/16) .* ( 1 .- u.^2 ).^2
triweight_kernel(u::Vector{T}) where {T} = T(35/32) .* ( 1 .- u.^2 ).^3
tricube_kernel(u::Vector{T}) where {T} = T(15/16) .* ( 1 .- abs.(u).^3 ).^3
gaussian_kernel(u::Vector{T}) where {T} = T(1 / sqrt(2*π)) .* exp.( T(-0.5) .* u.^2 )
cosine_kernel(u::Vector{T}) where {T} = T(π / 4) .* cos.( T(π / 2) .* u )
logistic_kernel(u::Vector{T}) where {T} = T(1) ./ (exp.(u) .+ T(2) .+ exp.(-u) )
sigmoid_kernel(u::Vector{T}) where {T} = T(2) ./ ( T(π) .* (exp.(u) .+ exp.(-u) ))