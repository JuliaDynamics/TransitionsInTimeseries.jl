"""
    permutation_entropy(; m = 3, τ = 1) → f

Return a function that given timeseries `x` it computes the normalized permutation entropy
of order `m` with time delay `τ`.
"""
function permutation_entropy(; m = 3, τ = 1)
    est = SymbolicPermutation(; m, τ)
    return x -> entropy_normalized(est, x)
end
