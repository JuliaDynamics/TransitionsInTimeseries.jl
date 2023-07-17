# TODO: Probably it makes sense to separate surrogate significance
# from computing the change metric? In this way alternative ways for significance
# can be created, such as exceeding a fixed value?

# In this way, we have one step that "calculates change metrics",
# and one step that "quantifies significance in change metrics".
# Does this tie in well with ChangePoints.jl?

"""
    TransitionsSurrogatesConfig(indicators, change_metrics [, surrogate]; kwargs...)

A configuration struct for TransitionsInTimeseries.jl that can be given to
[`estimate_transitions`](@ref). It contains all information necessary to perform the
basic [Workflow](@ref) of the package to detect significant transitions in the input
timeseries, using the method of surrogate testing to quantify significance:

2. Estimate the timeseries of an indicator by sliding a window over the input timeseries.
3. Estimate changes for an indicator by sliding a window over its timeseries.
4. Generate many surrogates that preserve important statistical properties of the original
   input timeseries.
5. Perform steps 1 and 2 for the surrogate timeseries (and for all provided indicators).
6. Estimate when an indicators timeseries shows significant change (trend, jump or anything else)
   when compared to the surrogate timeseries when compared to the surrogate data.
   Significance is estimated from the p-values of the real data vs surrogate data.

## Arguments

- `indicators::AbstractVector{<:Function}` a vector of indicators.
  Some indicators typically used in
  the literature are listed in the documentation section on [indicators](@ref indicators).
  The analysis is performed efficiently for all indicators given.
- `change_metrics` change metrics corresponding to the given indicators. If given
  a function, the same function is used for all indicators. Otherwise, it should be
  a vetor of functions of the same size as `indicators`, using change metric for its
  corresponding indicator. Some change metrics typically used in
  the literature are listed in the documentation section on [change metrics](@ref change_metrics).
- `surrogate::Surrogate` the method to use to generate surrogates of the input timeseries.
  This is an optional argument that defaults to `RandomFourier()`, see [Surrogates](@ref)
  for alternative options.

Both indicators and change metrics are generic Julia functions that input an
`x::AbstractVector` and output an `s::Real`. Any appropriate function may be given and
see [making custom indicators/change metrics](@ref own_indicator) in the documentation
for more information.

## Keyword arguments

- `width_ind::Int=100, stride_ind::Int=1`: width and stride given to the [`WindowViewer`](@ref)
  to compute the indicator from the input timeseries.
- `width_cha::Int=50, stride_cha::Int=1`: width and stride given to the [`WindowViewer`](@ref)
  to compute the change metric timeseries from the indicator timeseries.
- `whichtime = midpoint`: The time vector corresponding to the indicators / change metric
  timeseriesm is obtained from `t` using the keyword `whichtime`. Options include:
    - `last`: use the last timepoint of each window
    - `midpoint`: use the mid timepoint of each time window
    - `first`: use first timepoint of each window
  In fact, the indicators time vector is computed simply via
  ```julia
  t_indicator = windowmap(whichtime, t; width_ind, stride_ind)
  t_change = windowmap(whichtime, t_indicator; width_cha, stride_cha)
  ```
  so any other function of the time window may be given to extract the time point itself,
  such as `mean` or `median`.

- `n_surrogates::Int = 10_000`: how many surrogates to create.
- `rng::AbstractRNG = Random`.default_rng()`: a random number generator for the surrogates.
- `tail::Symbol = :both`: kind of tail test to do (one of `:left, :right, :both`) when
  estimating the p-value from the distribution of surrogate data.

- `T = Float64`: Element type of input timeseries to initialize some computations.

"""
struct TransitionsSurrogatesConfig{F<:Function, G<:Function, S<:Surrogate, W<:Function, R<:AbstractRNG}
    indicators::Vector{F}
    change_metrics::Vector{G}
    surrogate::S
    width_ind::Int
    stride_ind::Int
    width_cha::Int
    stride_cha::Int
    whichtime::W
    n_surrogates::Int
    rng::R
    tail::Symbol
end

function TransitionsSurrogatesConfig(
        indicators, change_metrics, surrogate = RandomFourier();
        width_ind = 100,
        stride_ind = 1,
        width_cha = 50,
        stride_cha = 1,
        whichtime =  midpoint,
        tail = :both,
        n_surrogates = 10_000,
        rng = Random.default_rng(),
        T = Float64,
    )
    # Sanity checks
    if !(indicators isa AbstractVector)
        indicators = [indicators]
    end
    if !(change_metrics isa AbstractVector)
        change_metrics = [change_metrics]
    end
    L = length(indicators)
    if length(change_metrics) ∉ (1, L)
        throw(ArgumentError("The amount of change metrics must be as many as the indicators, or only 1."))
    end
    # Last step: precomputable functions
    indicators = precompute_metrics(indicators, 1:T(width_ind))
    change_metrics = precompute_metrics(change_metrics, 1:T(width_cha))

    return TransitionsSurrogatesConfig(
        indicators, change_metrics, surrogate,
        width_ind, stride_ind, width_cha, stride_cha, whichtime, n_surrogates, rng, tail,
    )
end



"""
    IndicatorsResults

A struct containing the output of [`estimate_transitions`](@ref), which is
the main computational part of TransitionsInTimeseries.jl.
It can be used for analysis, visualization, and further given to [`transition_flags`](@ref).

It has the following fields that the user may access

- `x`: the input timeseries.
- `t`: the time vector of the input timeseries.

- `indicators::Vector{Function}`: indicators used in the processing.
- `x_indicator`, the indicator timeseries (matrix with each column one indicator).
- `t_indicator`, the time vector of the indicator timeseries.

- `change_metrics::Vector{Function}`: change metrics used in the processing.
- `x_change`, the change metric timeseries (matrix with each column one change metric).
- `t_change`, the time vector of the change metric timeseries.
- `pvalues`, the p-value of the change metrics w.r.t. the surrogates.
  It is a 2-dimensional array, where first dimension = time, second dimension = change
  metric. I.e. `pval[:, k]` will give the timeseries of p-value of the `k`-th change metric.
"""
struct IndicatorsResults{TT, T<:Real, X<:Real, XX<:AbstractVector{X},
        F<:Function, Z<:Function, S<:Surrogate}
    t::TT # original time vector; most often it is `Base.OneTo`.
    x::XX

    indicators::Vector{F}
    t_indicator::Vector{T}
    x_indicator::Matrix{X}

    change_metrics::Vector{Z}
    t_change::Vector{T}
    x_change::Matrix{X}
    pvalues::Matrix{X}
    surrogate_method::S
end

function Base.show(io::IO, ::MIME"text/plain", res::IndicatorsResults)
    println(io, "IndicatorsResults")
    descriptors = [
        "input timeseries" => summary(res.x),
        "indicators" => [nameof(i) for i in res.indicators],
        "change metrics" => [nameof(c) for c in res.change_metrics],
        "surrogate" => res.surrogate_method,
    ]
    padlen = maximum(length(d[1]) for d in descriptors) + 3
    for (desc, val) in descriptors
        println(io, rpad(" $(desc): ", padlen), val)
    end
end