#####################################################
############### Fourier surrogates ##################
#####################################################

# Shift each frequency content by a random phase.
function generate_fourier_surrogate(x::Vector{T}) where {T<:Real}
    F = rfft(x)
    return irfft( F .* exp.( 2*π*im .* rand(length(F)) ), length(x) )
end

"""
    generate_stacked_fourier_surrogates(X::AbstractArray, ns::Int)

Generate ns Fourier surrogates for each line of X.
If X is a CuArray, computation takes place in GPU.
Aformentionned speedup particularly significant if X big or ns large.
"""
function generate_stacked_fourier_surrogates(x::Vector{T}, ns::Int) where {T<:Real}
    S = zeros(T, ns, length(x))
    for i in axes(S, 1)
        S[i, :] = generate_fourier_surrogate(x)
    end
    return S
end

function generate_stacked_fourier_surrogates(X::Matrix{T}, ns::Int) where {T<:Real}
    nx, nt = size(X)
    stacked_surrogates = zeros(T, nx*ns, nt)
    for i in 1:nx
        stacked_surrogates[(i-1)*ns+1:i*ns, :] = generate_fourier_surrogates(X[i, :], ns)
    end
    return stacked_surrogates
end

function generate_stacked_fourier_surrogates(X::CuArray{T, 2}, ns::Int) where {T<:Real}
    nx, nt = size(X)
    Fcuda = repeat( CUDA.CUFFT.rfft( X, 2 ), inner=(ns, 1) )
    stacked_surrogates = CUDA.CUFFT.irfft( Fcuda .* exp.(2*π*im .* CUDA.rand(nx*ns, size(Fcuda,2)) ), nt, 2 )
    return stacked_surrogates
end

#####################################################
################# Ridge regression ##################
#####################################################

"""
    ridge_regression(y::AbstractArray, t::Vector; lambda=0::Real)

Ridge regression with regularization parameter lambda.
Set lambda=0 to recover linear regression.
""" 
function ridge_regression(y::Vector{T}, t::Vector{T}; lambda=0::Real) where {T<:Real}
    t = t .- t[1]
    T_bias_ext = hcat( t, ones(length(t)) )'
    return inv(T_bias_ext * T_bias_ext' + lambda .* I(2)) * T_bias_ext * y
end

function ridge_regression(Y::Matrix{T}, t::Vector{T}; lambda=0::Real) where {T<:Real}
    t = t .- t[1]
    T_bias_ext = hcat( t, ones(length(t)) )'
    return inv(T_bias_ext * T_bias_ext' + lambda .* I(2)) * T_bias_ext * transpose(Y)
end

function ridge_regression_slope(Y::Matrix{T}, t::Vector{T}; lambda=0::Real) where {T<:Real}
    return ridge_regression(Y, t; lambda)[1, :]
end

# GPU version of ridge regression. Y must be nt x ns.
function ridge_regression(Y::CuArray{T, 2}, t::Vector{T}; lambda=0::Real) where {T<:Real}
    t = t .- t[1]
    T_bias_ext = hcat( t, ones(length(t)) )'
    return CuArray( inv(T_bias_ext * T_bias_ext' + lambda .* I(2)) ) * CuArray( T_bias_ext ) * permutedims(Y)
end

function ridge_regression_slope(Y::CuArray{T, 2}, t::Vector{T}; lambda=0::Real) where {T<:Real}
    return ridge_regression(Y, t; lambda)[1, :]
end

#####################################################
################### Kendall tau #####################
#####################################################

function slide_kendall_tau(
    x::Vector{T},
    t::Vector{T},
    hw::Int,
    stride::Int,
) where {T<:Real}

    nx = length(x)
    kt = fill(NaN, nx)
    for i in (hw+1):stride:(nx-hw)
        x_wndwd = centered_wndw( x, i, hw )
        t_wndwd = centered_wndw( t, i, hw )
        kt[i] = corkendall( t_wndwd, x_wndwd )
    end
    return kt
end

function slide_kendall_tau(t::Vector{T}, S::Matrix{T}, hw::Int, stride::Int) where {T<:Real}
    return mapslices( x -> slide_kendall_tau(t, x, hw, stride), S, dims = 2 )
end

#####################################################
################### Kendall tau #####################
#####################################################

function percentile_significance(ref_stat::Vector{T}, sur_stat::Matrix{T}) where {T<:AbstractFloat}
    p = fill(NaN, length(ref_stat))
    for i in eachindex(ref_stat)
        if !isnan(ref_stat[i])
            p[i] = get_percentile( ref_stat[i], sur_stat[:,i] )
        end
    end
    return p
end

function percentile_significance(ref_stat::CuArray{T, 2}, sur_stat::CuArray{T,2}) where {T<:Real}
    return reduce( +, sur_stat .< ref_stat, dims=1 ) ./ size(sur_stat, 1)
end

function get_percentile(ref_stat::T, sur_stat::Vector{T}) where {T<:Real}
    return sum(sur_stat .< ref_stat) / length(sur_stat)
end

function get_percentile(ref_stat::CuArray{T,2}, sur_stat::CuArray{T,2}, ns::Int) where {T<:Real}
    rep_ref = repeat( ref_stat, inner=(ns, 1) )
    return reduce( +, sur_stat .- rep_ref .< T(0), dims=1 ) ./ size(sur_stat, 1)
end

function get_percentile(ref_stat::CuArray{T,2}, sur_stat::CuArray{T,2}, ns::Int, nx::Int) where {T<:Real}
    rep_ref = repeat( ref_stat, inner=(ns, 1) )
    return ( CuArray( kron(I(nx), ones(ns)') ) * (sur_stat .- rep_ref .< T(0))) ./ T(size(sur_stat, 1))
end

#####################################################
############# Slide trend estimation ################
#####################################################

function slide_idtrend(
    X::CuArray{T, 2},
    t::Vector{T},
    p::WindowingParams,
    estimator,
    wndw;
    kwargs...,
) where {T<:Real}

    nl, nt = size(X)
    strided_idx = wndw(p.Nwndw, p.Nstrd, nt)
    nidx = length(strided_idx)

    idtrend = fill(T(NaN), nl, nidx)
    for j1 in eachindex(strided_idx)
        j2 = strided_idx[j1]
        idtrend[:, j1] = Array( estimator( wndw( X, j2, p.Nwndw ), wndw( t, j2, p.Nwndw ); kwargs... ) )
    end
    return CuArray( idtrend )
end

function slide_idtrend(
    X::Matrix{T},
    t::Vector{T},
    p::WindowingParams,
    estimator,
    wndw;
    kwargs...,
) where {T<:Real}

    nl, nt = size(X)
    strided_idx = wndw(p.Nwndw, p.Nstrd, nt)
    nidx = length(strided_idx)

    idtrend = fill(T(NaN), nl, nidx)
    for j1 in eachindex(strided_idx)
        j2 = strided_idx[j1]
        idtrend[:, j1] = estimator( wndw( X, j2, p.Nwndw ), wndw( t, j2, p.Nwndw ); kwargs... )
    end
    return idtrend
end