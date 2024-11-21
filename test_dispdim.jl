using GLMakie, LinearAlgebra, DSP

include("functions.jl")

mfilter(f::Real) = digitalfilter(Lowpass(Float64(f), fs=1/dt), Butterworth(2))

function gauss(N::Integer)
    x = range(start=-3.0, stop=3.0, length=N)
    return exp.(-x .^ 2)
end

(dt, w0) = loaddata()

w = filtfilt(mfilter(2.0), w0)

L = round(Int, 5.0 / dt)

r1 = zeros(size(w, 1) - L + 1)
r2 = zeros(size(w, 1) - L + 1)
s2 = zeros(size(w, 1) - L + 1)
s3 = zeros(size(w, 1) - L + 1)
for i = eachindex(r1)
    S = svd(w[i:i+L-1, :])
    r1[i] = (S.S[2] - S.S[3]) / (S.S[1] - S.S[3])
    r2[i] = sum(S.S)/maximum(S.S)
    s2[i] = S.S[2]/S.S[1]
    s3[i] = S.S[3]/S.S[1]
end


fig = Figure()
ax = Axis(fig[1, 1])
lines!(w[:, 1] ./ maximum(abs, w), color=:gray60)
lines!(w[:, 2] ./ maximum(abs, w), color=:gray60)
lines!(w[:, 3] ./ maximum(abs, w), color=:gray60)
lines!(L .+ eachindex(s2), s2, color=:blue)
lines!(L .+ eachindex(s3), s3, color=:blue, linestyle=:dash)
# lines!(L .+ eachindex(r1), r1, color=:red, linestyle=:dash)
lines!(L .+ eachindex(r2), r2, color=:red)
fig
