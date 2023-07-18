"""
    difference_of_means(x::AbstractArray)

Return the absolute difference of the means of the first and second halfs of `x`.
`difference_of_means` can be used as a change metric focused on value differences.
Creating similar statistical differences using other moments instead of `mean` is trivial.
In fact, the source of `difference_of_means` is just:
```julia
# assumes 1-based indexing
x1 = view(x, 1:n÷2)
x2 = view(x, (n÷2 + 1):n)
return abs(mean(x1) - mean(x2))
```
"""
function difference_of_means(x::AbstractArray)
    # assumes 1-based indexing
    x1 = view(x, 1:n÷2)
    x2 = view(x, (n÷2 + 1):n)
    return abs(mean(x1) - mean(x2))
end
