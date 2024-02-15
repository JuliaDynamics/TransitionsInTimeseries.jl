"""
    difference_of_means(x::AbstractArray)

Return the absolute difference of the means of the first and second halfs of `x`.
`difference_of_means` can be used as a change metric focused on value differences.
Creating similar statistical differences using other moments instead of `mean` is trivial.
In fact, the source of `difference_of_means` is just:
```julia
# assumes 1-based indexing
n = length(x)
x1 = view(x, 1:n÷2)
x2 = view(x, (n÷2 + 1):n)
return abs(mean(x1) - mean(x2))
```

`difference_of_means` can also sensibly be used for windows of size 2,
in which case the change metric timeseries is the same as the `abs.(diff(...))`
of the indicator timeseries.
"""
function difference_of_means(x::AbstractArray)
    if length(x) == 2
        return abs(first(x) - last(x))
    end
    # assumes 1-based indexing
    n = length(x)
    x1 = view(x, 1:n÷2)
    x2 = view(x, (n÷2 + 1):n)
    return abs(mean(x1) - mean(x2))
end

function difference_of_maxes(x::AbstractArray)
    if length(x) == 2
        return abs(first(x) - last(x))
    end
    # assumes 1-based indexing
    n = length(x)
    x1 = view(x, 1:n÷2)
    x2 = view(x, (n÷2 + 1):n)
    return abs(maximum(x1) - maximum(x2))
end