using GLMakie, LinearAlgebra, DSP

include("functions.jl")

(dt, w) = loaddata()

mfilter(f::Real) = digitalfilter(Lowpass(Float64(f), fs=1/dt), Butterworth(2))

freqlist = 0.1:0.1:5.0
L = round(Int, 1.0/dt)

amap = zeros(length(w)-L+1, length(freqlist))
for ifreq = eachindex(freqlist)
    local wf = filtfilt(mfilter(freqlist[ifreq]), w)
    dw = w - wf
    for i = axes(amap, 1)
        amap[i, ifreq] = norm(dw[i:i+L-1])/norm(w[i:i+L-1])
    end
end

wf0 = filtfilt(mfilter(1), w)
dwf0 = w - wf0
a0 = map(i->log1p(norm(dwf0[i:i+L-1])/norm(w[i:i+L-1])), 1:length(w)-L+1)

fig = Figure()
ax = Axis(fig[1,1])
lines!(w./maximum(abs, w), color=:gray60)
lines!(wf0./maximum(abs, w), color=:black)
lines!((w-wf0)./maximum(abs, w), color=:blue)
ax3 = Axis(fig[2,1])
lines!(a0)
ax2 = Axis(fig[3,1])
heatmap!(axes(amap, 1), freqlist, amap, colormap=:jet)
# axislegend()
linkxaxes!(ax, ax2, ax3)
fig
