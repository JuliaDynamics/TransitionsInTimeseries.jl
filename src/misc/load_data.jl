"""

    load_linear_vs_doublewell()

Load prototypical data from a linear and a double-well model to test some indicators.
"""
function load_linear_vs_doublewell()
    url = "https://raw.githubusercontent.com/JuliaDynamics/JuliaDynamics/master/"*
        "timeseries/linear_vs_doublewell.csv"
    tmp = Downloads.download(url)
    data = readdlm(tmp, ',', Float64, skipstart = 1)
    return view(data, :, 1), view(data, :, 2), view(data, :, 3) # t, xlin, xnlin
end
