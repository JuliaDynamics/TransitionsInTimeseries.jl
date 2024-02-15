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