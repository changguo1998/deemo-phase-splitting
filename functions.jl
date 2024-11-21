using DelimitedFiles, LinearAlgebra

function loaddata()
    dt = 0.01
    w = readdlm("data.txt")
    return (dt, w)
end

# basic

apply_window(func_kernel::Function, w::AbstractVector, L::Integer) = map(i -> func_kernel(w[i:i+L-1]), 1:(length(w)-L+1))

# STA/LTA

@inline kernal_stalta(w::AbstractVector, p::Integer) = log1p(norm(w .- mean(w), p))
kernal_stalta(p::Integer) = w -> kernal_stalta(w, p)

function stalta(w::AbstractVector, Lformer::Integer, Llatter::Integer, p::Integer=2)
    rformer = apply_window(kernal_stalta(2), w, Lformer)
    rlatter = apply_window(kernal_stalta(2), w[Lformer:end], Llatter)
    L = min(length(rformer), length(rlatter))
    r = rformer[1:L] ./ rlatter[1:L]
    return (r, argmin(r) + Lformer)
end

# dimension