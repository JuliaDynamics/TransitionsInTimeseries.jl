#####################################################
# Default choices
#####################################################

default_n_surrogates() = 10_000
default_surrogate_method() = RandomFourier()

function default_window_width(x)
    w = max(length(x)÷100, 10)
    if w < 50
        @warn "The window width chosen by default is w = $w and might be too small!"
    else
        @info "The window width chosen by default is w = $w."
    end
    return w
end

function default_window_stride(x)
    s = 1
    @info "The window stride chosen by default is s = $s."
    return s
end