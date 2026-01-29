function get_funcname(func_type::SchwartzFunc{T}; details::Bool=false) where {T<:Real}
    !details && return "Schwartz"
    return "exp(-func_type.A * abs((x - func_type.μ) / func_type.σ)^func_type.d)"
end

function get_funcname(func_type::RationalFuncPolesRepresent{T};
                      details::Bool=false) where {T<:Real}
    !details && return "Rational"
    @assert length(func_type.poles) == length(func_type.amplitudes) ==
            length(func_type.Hmask)
    @assert length.(func_type.poles) .>= 1
    res = "real("
    for i in 1:length(func_type.poles)
        for j in 1:length(func_type.poles[i])
            res += "$(round(func_type.amplitudes[i][j], digits=2))/(x - $(round(func_type.poles[i][j], digits=2)))^$(i)"
            i < length(func_type.poles) && (res += " + ")
        end
    end
    res += ")"
    return res
end

function get_funcname(func_type::DRationdlFunc{T}; details::Bool=false) where {T<:Real}
    return "$(round(func_type.d, digits=2))-Rationdl"
end

function get_funcname(func_type::LogRationalFunc{T}; details::Bool=false) where {T<:Real}
    return "LogRational"
end

function get_funcname(func_type::MixedFunc{T}; details::Bool=false) where {T<:Real}
    res = ""
    !isnothing(func_type.swf) &&
        (res += get_funcname(func_type.swf; details=details) + " + ")
    !isnothing(func_type.rtfpr) &&
        (res += get_funcname(func_type.rtfpr; details=details) + " + ")
    !isnothing(func_type.drtf) &&
        (res += get_funcname(func_type.drtf; details=details) + " + ")
    !isnothing(func_type.logrtf) &&
        (res += get_funcname(func_type.logrtf; details=details) + " + ")
    !endswith(res, " + ") && (res = res[1:(end - 1)])
    return res
end

function get_algname(dm::DeModeMethod, pola::PolationMethod,
                     trans::DiscreteTransMethodexport)
    res = ""
    if dm isa AsymptoticDeMode
        res += "ASY DeMode, degree=$(dm.degree) + "
    elseif dm isa AAADeMode
        res += "AAA DeMode, max_degree=$(dm.max_degree) + "
    end

    if pola isa NoPolation
        res += "No Polation + "
    elseif pola isa InterPolation
        res += "Hann Interpolation + "
    elseif pola isa Extrapolation
        res += "Hermitian Extrapolation + "
    end

    if trans isa FFTTrans
        res += "FFT Trans with pad rate = $pad_rate"
    elseif trans isa FIRTrans
        res += "FIR Trans"
    end

    return res
end

function write_setting(func_type::TestFunc{T}, L0_vec::Vector{Real}, grid_gap::Real,
                       points_density::Int;
                       file_place::String=".", dm::DeModeMethod, pola::PolationMethod,
                       trans::DiscreteTransMethodexport) where {T<:Real}
    open(joinpath(file_place, "setting.json"), "w") do f
        write(f, "Used L0: $(round.(L0_vec, digits=2))\n")
        write(f, "grid_gap: $(round(grid_gap, digits=2))\n")
        write(f, "Points density: $(points_density)\n")
        write(f, "func_type: $(get_funcname(func_type; details=true))\n")
        return write(f, "alg: $(get_algname(dm, pola, trans))\n")
    end
end

# ============================== Bench Method =========================================

function loss_bench_plot(L0_vec, maxerr_herm_vec, maxerr_hann_vec, maxerr_trunc_vec,
                         L2relerr_herm_vec, L2relerr_hann_vec, L2relerr_trunc_vec,
                         dm::DeModeMethod, pola::PolationMethod,
                         trans::DiscreteTransMethodexport)
    title1 = "Max error \n $(get_funcname(func_type; details=false)) \n$(get_algname(dm, pola, trans))"
    title2 = "L2 relative error \n $(get_funcname(func_type; details=false)) \n$(get_algname(dm, pola, trans))"

    fig1 = Figure(; size=(800, 600))
    Label(fig1[0, 1], title1; tellwidth=false)
    ax1 = Axis(fig1[1, 1];
               xlabel="L0",
               ylabel="max error",
               title="",
               xscale=log10,
               yscale=log10,)

    lines!(ax1, L0_vec, maxerr_herm_vec; label="Hermite")
    lines!(ax1, L0_vec, maxerr_hann_vec; label="Hann")
    lines!(ax1, L0_vec, maxerr_trunc_vec; label="Trunc")

    axislegend(ax1)

    fig2 = Figure(; size=(800, 600))
    Label(fig2[0, 1], title2; tellwidth=false)

    ax2 = Axis(fig2[1, 1];
               xlabel="L0",
               ylabel="L2 relative error",
               title="",
               xscale=log10,
               yscale=log10,)

    lines!(ax2, L0_vec, L2relerr_herm_vec; label="Hermite")
    lines!(ax2, L0_vec, L2relerr_hann_vec; label="Hann")
    lines!(ax2, L0_vec, L2relerr_trunc_vec; label="Trunc")

    axislegend(ax2)
    return fig1, fig2
end

# The DeModeMethod `dm` here does not provide grid, so we need to generate it ourselves.
function loss_bench_report(funct_type::TestFunc{T}, dm::DeModeMethod;
                           trans::DiscreteTransMethodexport=FIRTrans(); L0_start::Real=2^2,
                           L0_rate::Real=2,
                           test_num::Int=14, point_density::Int=2^5, is_saveset::Bool=false,
                           file_place::String=".")
    @assert test_num >= 1
    @assert point_density >= 1
    h = T(1/32)
    L0_vec = T.(L0_start * L0_rate .^ (0:(test_num - 1)))
    m = length(L0_vec)
    maxerr_herm_vec = zeros(T, m)
    maxerr_hann_vec = zeros(T, m)
    maxerr_trunc_vec = zeros(T, m)
    L2relerr_herm_vec = zeros(T, m)
    L2relerr_hann_vec = zeros(T, m)
    L2relerr_trunc_vec = zeros(T, m)

    for j in 1:m
        L0 = L0_vec[j]
        println("Running test $j of $m, L0 = $L0")
        N = point_num * 2 * L0 + 1
        h = T(2L0 / (N - 1))
        x = T.((-N ÷ 2):(N ÷ 2)) .* h
        local dm1
        if dm isa AsymptoticDeMode
            dm1 = AsymptoticDeMode(; grid=x, degree=dm.degree, n=dm.n, is_print=dm.is_print)
        elseif dm isa AAADeMode
            dm1 = AAADeMode(; grid=x, max_degree=dm.max_degree, show_poles=dm.show_poles)
        else
            error("Unsupported DeModeMethod: $dm")
        end

        f = [origfunc(xi, funct_type) for xi in x]
        H_exact = [Hfunc(xi, funct_type) for xi in x]

        H_trunc = hilbert(f; dm=dm1, pola=NoPolation(), trans=trans)
        H_hann = hilbert(f; dm=dm1, pola=InterPolation(; δ=round(Int, N / 64)), trans=trans)
        H_herm = hilbert(f; dm=dm1, pola=Extrapolation(; n=round(Int, 3 / h), h=h),
                         trans=trans)

        dH_herm = abs.(H_herm - H_exact)
        dH_hann = abs.(H_hann - H_exact)
        dH_trunc = abs.(H_trunc - H_exact)
        err_herm = maximum(dH_herm)
        err_hann = maximum(dH_hann)
        err_trunc = maximum(dH_trunc)
        @show err_herm
        @show err_hann
        @show err_trunc
        maxerr_herm_vec[j] = err_herm
        maxerr_hann_vec[j] = err_hann
        maxerr_trunc_vec[j] = err_trunc

        L2f = sum(f .^ 2) - (f[1]^2 + f[end]^2) / 2
        L2f = sqrt(L2f * h)
        @show L2f

        L2err_herm = sum(dH_herm .^ 2) - (dH_herm[1]^2 + dH_herm[end]^2) / 2
        L2err_herm = sqrt(L2err_herm * h)
        L2relerr_herm_vec[j] = L2err_herm / L2f

        L2err_hann = sum(dH_hann .^ 2) - (dH_hann[1]^2 + dH_hann[end]^2) / 2
        L2err_hann = sqrt(L2err_hann * h)
        L2relerr_hann_vec[j] = L2err_hann / L2f

        L2err_trunc = sum(dH_trunc .^ 2) - (dH_trunc[1]^2 + dH_trunc[end]^2) / 2
        L2err_trunc = sqrt(L2err_trunc * h)
        L2relerr_trunc_vec[j] = L2err_trunc / L2f

        @show L2relerr_herm_vec[j]
        @show L2relerr_hann_vec[j]
        @show L2relerr_trunc_vec[j]
    end

    return loss_bench_plot(L0_vec, maxerr_herm_vec, maxerr_hann_vec, maxerr_trunc_vec,
                           L2relerr_herm_vec, L2relerr_hann_vec, L2relerr_trunc_vec,
                           dm, pola, trans)
end
