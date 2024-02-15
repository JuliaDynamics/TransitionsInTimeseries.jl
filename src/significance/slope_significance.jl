"""
    SlopeChangeSignificance(; moe_slope, moe_offset, slope_diff = moe_slope, pvalue = 0.05)

Test whether the result of [`SlopeChangeResults`](@ref) is statistically significant.

Two tests are done:
1. Check whether the _margin of error_ of the fitted parameters `a, b, c, d`
   of the two linear segments `a + b*t, c + d*t`
   is less than the specified margins of error, for a chosen `pvalue`.
   To disable this test give `Inf` to the `moe` keywords.
2. Test that the two slopes `b, d` have difference greater than `slope_diff`.
   To disable this test give `0` to the `slope_diff` keyword.

The Boolean `&` of the above two is the final test.

The margin of error is simply half the size of the confidence interval,
also known as radius of the confidence interval.
"""
@kwdef struct SlopeChangeSignificance <: Significance
    moe_slope::Float64
    moe_offset::Float64
    slope_diff::Float64 = moe_slope
    pvalue::Float64 = 0.05
end

function significant_transitions(res::SlopeChangeResults, signif::SlopeChangeSignificance)
    moe = LsqFit.margin_error(res.lsqfit, signif.pvalue)
    moeflag = (moe[1] ≤ signif.moe_offset &&
        moe[3] ≤ signif.moe_offset &&
        moe[2] ≤ signif.moe_slope &&
        moe[4] ≤ signif.moe_slope)

    slopeflag = abs(res.fitparams[2] - res.fitparams[4]) > moe_slope
    return moeflag && slopeflag
end