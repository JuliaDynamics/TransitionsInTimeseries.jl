using CairoMakie

T = Float32
p = WindowingParams(1f0, 1f0, 1f0, 20, 1)
list = [uniform_kernel, triangular_kernel, parabolic_kernel, biweight_kernel, triweight_kernel,
        tricube_kernel, gaussian_kernel, cosine_kernel, logistic_kernel, sigmoid_kernel]

fig = Figure(resolution = (200, 2000))
axs = [Axis(fig[i,1], title = string(list[i])) for i in eachindex(list)]

for i in eachindex(list)
    lines!( axs[i], scale_x_kernel(T, p), scaled_kernel(list[i], T, p) )
end
fig