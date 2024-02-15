"""
    SlopeChangeSignificance(; moe_slope, moe_offset, pvalue = 0.05)

Test whether the result of [`SlopeChangeResults`](@ref) are statistically significant
by checking whether the _margin of error_ of the fitted parameters `a, b, c, d`
is less than the specified margins of error, for a chosen `pvalue`.
The margin of error is simply half the size of the confidence interval,
also known as radius of the confidence interval.
"""
@kwdef struct SlopeChangeSignificance <: Significance
    moe_slope::Float64
    moe_offset::Float64
    pvalue::Float64 = 0.05
end

function significant_transitions(res::SlopeChangeResults, signif::SlopeChangeSignificance)
    moe = LsqFit.margin_error(res.lsqfit, signif.pvalue)
    flag = (moe[1] ≤ signif.moe_offset &&
        moe[3] ≤ signif.moe_offset &&
        moe[2] ≤ signif.moe_slope &&
        moe[4] ≤ signif.moe_slope)

    return flag
end