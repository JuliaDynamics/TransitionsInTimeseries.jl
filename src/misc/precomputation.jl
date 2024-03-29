"""
    PrecomputableFunction

Supertype of structs containing the necessary field to precompute a `::Function` by:

```julia
precompute(f::PrecomputableFunction, t)
```
"""
abstract type PrecomputableFunction end

"""
    precompute(f, t)

Precompute the function `f` given a time vector `t`.
"""
precompute(f, t) = f

precompute(f::PrecomputableFunction, t) = error("method not implemented")
