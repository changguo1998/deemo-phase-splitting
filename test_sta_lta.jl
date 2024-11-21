using GLMakie, LinearAlgebra

include("functions.jl")

(dt, w) = loaddata()

fw(x) = log1p(norm(x, 2))

llong = round(Int, 3.0/dt)
lshort = round(Int, 2.0/dt)

rff = zeros(length(w)-llong+1)
rfl = zeros(length(w)-llong-lshort+1)
rll = zeros(length(w)-llong+1)

for i = eachindex(rff)
    rff[i] = fw(w[i+llong-lshort:i+llong-1])/fw(w[i:i+llong-1])
    rll[i] = fw(w[i:i+lshort-1])/fw(w[i:i+llong-1])
end

for i = eachindex(rfl)
    rfl[i] = fw(w[i:i+lshort-1])/fw(w[i+lshort-1:i+lshort+llong-2])
end

tff = eachindex(rff) .+ llong .- 1
tfl = eachindex(rfl) .+ lshort .- 1
tll = eachindex(rll)

pff = argmin(rff) + llong - 1
pfl = argmin(rfl) + lshort - 1
pll = argmin(rll) + llong - lshort

fig = Figure()
ax = Axis(fig[1,1])
lines!(w./maximum(abs, w), color=:gray60)
lines!(tff, rff, color = :black, label="FF")
lines!(tfl, rfl, color = :blue, label="FL")
lines!(tll, rll, color = :red, label="LL")
vlines!([pff, pfl, pll])
text!([pff, pfl, pll], zeros(3), text=["Pff", "Pfl", "Pll"])
axislegend()
fig
