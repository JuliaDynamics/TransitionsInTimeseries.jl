#####################################################
# Default choices
#####################################################

const DEFAULT_N_SURROGATES = 10_000
const DEFAULT_SURROGATE_METHOD = RandomFourier()
const DEFAULT_WINDOW_STRIDE = 1

"""
    default_window_width(x)

Return the default window width chosen to be `length(x) ÷ 100` and at least `10`.
"""
function default_window_width(x)
    w = max(length(x)÷100, 10)
    if w < 50
        @warn "The window width chosen by default is w = $w and might be too small!"
    end
    return w
end