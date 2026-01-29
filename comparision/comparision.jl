using OpenHilbert
using LinearAlgebra, SpecialFunctions
using CairoMakie

func(x) = exp(-x^2)
Hfunc(x) = 2 / sqrt(π) * dawson(x)
func_name = "f(x) = exp(-x^2)"

function compare(; pad_rate::Int=2, func_name::String="")
    L0_vec = [30, 60, 120]
    m = length(L0_vec)
    nN = 9

    h_vec = 2.0 .^ ((-nN):-1)
    maxerr_herm_mat = zeros(m, nN)
    maxerr_hann_mat = zeros(m, nN)
    maxerr_trunc_mat = zeros(m, nN)
    for j in 1:m
        # N chosen so that h = 2L0/N matches h_vec = 2^(-nN:-1)
        # This gives N = L0 * 2^(nN+1:-1:2)
        N_vec = (2 .^ ((nN + 1):-1:2)) .* L0_vec[j]
        for k in 1:nN
            println("Running test $j of $m, L0 = $L0_vec[j], N = $N_vec[k]")
            L0 = L0_vec[j]
            N = N_vec[k]
            h = 2L0 / N
            @assert h_vec[k] == h
            x = (-L0) .+ (0:(N - 1)) .* h
            f = func.(x)
            H_exact = Hfunc.(x)
            pad = pad_rate * N
            H_herm = hilbert(f, FFTVecExtrapolation(; n=N÷4, pad=pad, h=h))
            if N < 128
                nδ = 32
            elseif N < 2^12
                nδ = 64
            else
                nδ = div(N÷2, 64)
            end
            H_hann = hilbert(f, FFTVecInterpolation(; δ=nδ, pad=pad))
            H_trunc = hilbert(f, FFTVecTrunc(; pad=pad))
            dH_herm = abs.(H_herm - H_exact)
            dH_hann = abs.(H_hann - H_exact)
            dH_trunc = abs.(H_trunc - H_exact)
            maxerr_herm_mat[j, k] = maximum(dH_herm)
            maxerr_hann_mat[j, k] = maximum(dH_hann)
            maxerr_trunc_mat[j, k] = maximum(dH_trunc)
        end
    end

    fig_trunc = comp_plot("trunc, pad_rate=$pad_rate, $func_name", h_vec, maxerr_trunc_mat)
    fig_hann = comp_plot("hann, pad_rate=$pad_rate, $func_name", h_vec, maxerr_hann_mat)
    fig_herm = comp_plot("herm, pad_rate=$pad_rate, $func_name", h_vec, maxerr_herm_mat)
    return fig_trunc, fig_hann, fig_herm
end

function comp_plot(method, h_vec, maxerr_mat)
    fig = Figure()
    ax = Axis(fig[1, 1];
              title=method,
              xlabel="h",
              ylabel="max error",
              xscale=log10,
              yscale=log10)

    # Color/marker spec for each L0 curve:
    # - L0=30: red diamonds
    # - L0=60: green circles
    # - L0=120: blue triangles
    # Makie `lines!` doesn't support markers; draw lines + scatter markers.
    lines!(ax, h_vec, maxerr_mat[1, :]; color=:red)
    scatter!(ax, h_vec, maxerr_mat[1, :];
             label="L0=30",
             color=:red,
             marker=:diamond,
             markersize=8)

    lines!(ax, h_vec, maxerr_mat[2, :]; color=:green)
    scatter!(ax, h_vec, maxerr_mat[2, :];
             label="L0=60",
             color=:green,
             marker=:circle,
             markersize=8)

    lines!(ax, h_vec, maxerr_mat[3, :]; color=:blue)
    scatter!(ax, h_vec, maxerr_mat[3, :];
             label="L0=120",
             color=:blue,
             marker=:utriangle,
             markersize=8)

    axislegend(ax)
    # Sanitize filename: method may contain '/', spaces, commas, etc. (e.g. "1/(1+x^2)")
    safe_method = replace(string(method), r"[^A-Za-z0-9._-]+" => "_")
    outpath = joinpath(@__DIR__, "comp_$(safe_method).png")
    save(outpath, fig)
    return fig
end

compare(; pad_rate=1000, func_name=func_name)
