function get_funcname(func_type::SchwartzFunc{T}; details::Bool=false) where {T<:Real}
    !details && return "Schwartz"
    return "exp(-$(func_type.A) * abs((x - $(func_type.μ)) / $(func_type.σ))^$(func_type.d))"
end

function get_funcname(func_type::RationalFuncPolesRepresent{T};
                      details::Bool=false) where {T<:Real}
    !details && return "Rational"
    @assert length(func_type.poles) == length(func_type.amplitudes) ==
            length(func_type.Hmask)
    @assert all(length.(func_type.poles) .>= 1)
    res = "real("
    first = true
    for i in 1:length(func_type.poles)
        for j in 1:length(func_type.poles[i])
            if !first
                res *= " + "
            end
            res *= "$(round(func_type.amplitudes[i][j], digits=2))/(x - $(round(func_type.poles[i][j], digits=2)))^$(i)"
            first = false
        end
    end
    res *= ")"
    return res
end

function get_funcname(func_type::DRationdlFunc{T}; details::Bool=false) where {T<:Real}
    !details &&
        return "$(length(func_type.even_order_vec)) orders even + $(length(func_type.odd_order_vec)) orders odd -Rationdl"
    res = "∑ 1/(σ² + (x-s)^2)^(p/2) for (p,s,σ) ∈($([(round(p, digits=2), round(s, digits=2), round(σ, digits=2)) for (p,s,σ) in zip(func_type.even_order_vec, func_type.even_shift_vec, func_type.even_σ_vec)]))"
    !isempty(func_type.odd_order_vec) &&
        (res *= "+ ∑ x/(σ² + (x-s)^2)^((p+1)/2) for (p,s,σ) ∈($([(round(p, digits=2), round(s, digits=2), round(σ, digits=2)) for (p,s,σ) in zip(func_type.odd_order_vec, func_type.odd_shift_vec, func_type.odd_σ_vec)]))")
    return res
end

function get_funcname(func_type::LogRationalFunc; details::Bool=false)
    !details && return "$(round(func_type.d, digits=2))-LogRational"
    return "x^2 / (x^2 + 2)^(1 + $(round(func_type.d, digits=2))/2) / log(x^2 + 2)"
end

function get_funcname(func_type::MixedFunc{T}; details::Bool=false) where {T<:Real}
    res = ""
    !isnothing(func_type.swf) &&
        (res *= get_funcname(func_type.swf; details=details) * " + ")
    !isnothing(func_type.rtfpr) &&
        (res *= get_funcname(func_type.rtfpr; details=details) * " + ")
    !isnothing(func_type.drtf) &&
        (res *= get_funcname(func_type.drtf; details=details) * " + ")
    !isnothing(func_type.logrtf) &&
        (res *= get_funcname(func_type.logrtf; details=details) * " + ")
    # Remove trailing " + " if present
    endswith(res, " + ") && (res = res[1:(end - 3)])
    return res
end

function get_algname(tdm::Union{TestDeMode,Nothing},
                     trans::Union{DiscreteTransMethod,Nothing}; details::Bool=false)
    res = ""
    if !isnothing(tdm)
        if tdm isa TestNoDeMode
            res *= "No DeMode + "
        elseif tdm isa TestAsy
            res *= "ASY DeMode + "

            details &&
                (res *= "ASY DeMode order_vec=$([round(order, digits=2) for order in tdm.order_vec]), ")

            if tdm.mode_length_rate > 0
                res *= "mode_length_rate = $(tdm.mode_length_rate) + \n"
            else # tdm.mode_length >= 0
                res *= "mode_length = $(tdm.mode_length) + \n"
            end

        elseif tdm isa TestVarLsq
            res *= "VarLsq DeMode + "

            details &&
                (res *= "VarLsq DeMode, length expand rate=$(tdm.rate), nseek_vec=$(tdm.nseek_vec), max_deg=$(tdm.max_deg), ")

            if tdm.start_length_rate > 0
                res *= "start_length_rate = $(tdm.start_length_rate) + \n"
            else # tdm.start_length >= 0
                res *= "start_length = $(tdm.start_length) + \n"
            end
        else
            error("Unsupported TestDeMode: $tdm")
        end
    end

    if !isnothing(trans)
        if trans isa OpenHilbert.FFTTrans
            res *= "FFT Trans with pad_rate = $(trans.pad_rate)"
        elseif trans isa OpenHilbert.FIRTrans
            res *= "FIR Trans"
        end
    end

    return res
end

function write_setting(func_type::TestFunc{T}, L0_vec::Vector{<:Real}, grid_gap::Real,
                       points_density::Int;
                       file_place::String=".",
                       tdm::Union{TestDeMode,Nothing}=nothing,
                       trans::Union{DiscreteTransMethod,Nothing}=nothing) where {T<:Real}
    open(joinpath(file_place, "setting.txt"), "w") do f
        write(f, "Used L0: $(round.(L0_vec, digits=2))\n")
        write(f, "grid_gap: $(round(grid_gap, digits=2))\n")
        write(f, "Points density: $(points_density)\n")
        write(f, "func_type: $(get_funcname(func_type; details=true))\n")
        return write(f, "alg: $(get_algname(tdm, trans))\n")
    end
end

# ============================== Bench Method =========================================

function loss_bench_plot(L0_vec, N_vec, maxerr_trunc_vec,
                         L2relerr_trunc_vec;
                         title1::String="", title2::String="")
    N_vec_str = ["10^$(round(log10(N), digits=2))" for N in N_vec][1:2:end]

    fig1 = Figure(; size=(800, 600))
    Label(fig1[0, 1], title1; tellwidth=false)
    ax1 = Axis(fig1[1, 1];
               xlabel="L0",
               ylabel="max error",
               title="",
               xscale=log10,
               yscale=log10,
               yminorticksvisible=true,
               yminorticks=IntervalsBetween(2))

    lines!(ax1, L0_vec, maxerr_trunc_vec; label="Trunc")

    ax1_top = Axis(fig1[1, 1];
                   xlabel="Point Number",
                   xaxisposition=:top,
                   xscale=log10,
                   yticksvisible=false,
                   ylabelvisible=false,
                   ylabel="",
                   yticklabelsvisible=false,
                   bottomspinevisible=false,
                   leftspinevisible=false,
                   rightspinevisible=false,
                   topspinevisible=true)
    linkxaxes!(ax1_top, ax1)
    ax1_top.xticks = (L0_vec[1:2:end], N_vec_str)

    axislegend(ax1)

    fig2 = Figure(; size=(800, 600))
    Label(fig2[0, 1], title2; tellwidth=false)

    ax2 = Axis(fig2[1, 1];
               xlabel="L0",
               ylabel="L2 relative error",
               title="",
               xscale=log10,
               yscale=log10,
               yminorticksvisible=true,
               yminorticks=IntervalsBetween(2))

    lines!(ax2, L0_vec, L2relerr_trunc_vec; label="Trunc")

    ax2_top = Axis(fig2[1, 1];
                   xlabel="Point Number",
                   xaxisposition=:top,
                   xscale=log10,
                   yticksvisible=false,
                   ylabelvisible=false,
                   ylabel="",
                   yticklabelsvisible=false,
                   bottomspinevisible=false,
                   leftspinevisible=false,
                   rightspinevisible=false,
                   topspinevisible=true)
    linkxaxes!(ax2_top, ax2)
    ax2_top.xticks = (L0_vec[1:2:end], N_vec_str)

    axislegend(ax2)
    return fig1, fig2
end

function loss_bench_report(func_type::TestFunc{T}, tdm::TestDeMode,
                           trans::DiscreteTransMethod, L0_start::Real,
                           L0_rate::Real, test_num::Int, point_density::Int,
                           is_saveset::Bool, file_place::String) where {T<:Real}
    @assert test_num >= 1
    @assert point_density >= 1
    h = T(1 / point_density)
    L0_vec = T.(L0_start * L0_rate .^ (0:(test_num - 1)))
    N_vec = [round(Int, point_density * L0) * 2 + 1 for L0 in L0_vec]
    x0 = T.((-N_vec[end] ÷ 2):(N_vec[end] ÷ 2)) .* h
    f0 = [origfunc(xi, func_type) for xi in x0]
    mid = N_vec[end] ÷ 2 + 1

    println("Calculating exact Hilbert transform...")
    if func_type isa LogRationalFunc
        @warn("Hilbert transform of LogRationalFunc does not have a closed form, using numerical approximation")
        # cal_Hlogrtf_nume computes on a larger interval (internally uses L = 10*L0),
        # so we slice its centered part to match the target grid x0 (length N_vec[end]).
        H_log = cal_Hlogrtf_nume(L0_vec[end], point_density, func_type.d)
        mid_log = length(H_log) ÷ 2 + 1
        H_exact0 = H_log[(mid_log - N_vec[end] ÷ 2):(mid_log + N_vec[end] ÷ 2)]
    elseif func_type isa MixedFunc && !isnothing(func_type.logrtf)
        @warn("Hilbert transform of LogRationalFunc does not have a closed form, using numerical approximation")
        func_type_log = MixedFunc(T; swf=func_type.swf, rtfpr=func_type.rtfpr,
                                  drtf=func_type.drtf,
                                  logrtf=nothing)
        H_exact0 = Hfunc(x0, func_type_log)
        H_log = cal_Hlogrtf_nume(L0_vec[end], point_density, func_type.logrtf.d)
        mid_log = length(H_log) ÷ 2 + 1
        H_log0 = H_log[(mid_log - N_vec[end] ÷ 2):(mid_log + N_vec[end] ÷ 2)]
        H_exact0 .+= H_log0
    else
        H_exact0 = Hfunc(x0, func_type)
    end

    println("Calculating Hilbert transform with numerical approximation...")
    m = length(L0_vec)
    maxerr_trunc_vec = zeros(T, m)
    L2relerr_trunc_vec = zeros(T, m)

    for j in 1:m
        L0 = L0_vec[j]
        println("Running test $j of $m, L0 = $L0 \n ================================================================\n")
        N = N_vec[j]
        x = T.((-N ÷ 2):(N ÷ 2)) .* h
        dm = test2demode(x, tdm)

        f = view(f0, (mid - N ÷ 2):(mid + N ÷ 2))
        H_exact = view(H_exact0, (mid - N ÷ 2):(mid + N ÷ 2))

        println("beigin NoPolation")
        H_trunc = hilbert(f; dm=dm, pola=NoPolation(), trans=trans)

        dH_trunc = abs.(H_trunc - H_exact)
        err_trunc = maximum(dH_trunc)
        @show err_trunc
        maxerr_trunc_vec[j] = err_trunc

        L2f = sum(f .^ 2) - (f[1]^2 + f[end]^2) / 2
        L2f = sqrt(L2f * h)
        @show L2f

        L2err_trunc = sum(dH_trunc .^ 2) - (dH_trunc[1]^2 + dH_trunc[end]^2) / 2
        L2err_trunc = sqrt(L2err_trunc * h)
        L2relerr_trunc_vec[j] = L2err_trunc / L2f

        @show L2relerr_trunc_vec[j]
        println("\n\n\n")
    end
    is_saveset &&
        write_setting(func_type, L0_vec, h, point_density; tdm=tdm, trans=trans,
                      file_place=file_place)
    title1 = "Max error \n $(get_funcname(func_type; details=false)) \n$(get_algname(tdm, trans))"
    title2 = "L2 relative error \n $(get_funcname(func_type; details=false)) \n$(get_algname(tdm, trans))"

    return loss_bench_plot(L0_vec, N_vec, maxerr_trunc_vec, L2relerr_trunc_vec;
                           title1=title1, title2=title2)
end

function cal_Hlogrtf_nume(L0::T, point_density::Int, d::Int) where {T<:Real}
    L = 100 * L0
    N = round(Int, point_density * L) * 2 + 1
    h = T(1 / point_density)
    x = T.((-N ÷ 2):(N ÷ 2)) .* h
    func(x) = x / (x^T(2) + T(2))^(T(1) + d / T(2)) / log(x^T(2) + T(2))
    f = [func(xi) for xi in x]
    # IMPORTANT: pass the grid so the transform matches the sampling step h.
    dm = NoDeMode()
    Hf = hilbert(f; dm=dm, pola=NoPolation(), trans=FFTTrans(; pad_rate=0)) .* x
    return Hf
end

function loss_bench_report(func::Function, Hfunc::Function, tdm::TestDeMode,
                           trans::DiscreteTransMethod, L0_start::Real,
                           L0_rate::Real, test_num::Int, point_density::Int;
                           T::Type{<:Real}=Float64)
    @assert test_num >= 1
    @assert point_density >= 1
    h = T(1 / point_density)
    L0_vec = T.(L0_start * L0_rate .^ (0:(test_num - 1)))
    N_vec = [round(Int, point_density * L0) * 2 + 1 for L0 in L0_vec]
    x0 = T.((-N_vec[end] ÷ 2):(N_vec[end] ÷ 2)) .* h
    f0 = [func(xi) for xi in x0]
    mid = N_vec[end] ÷ 2 + 1
    H_exact0 = Hfunc(x0)

    println("Calculating Hilbert transform with numerical approximation...")
    m = length(L0_vec)
    maxerr_trunc_vec = zeros(T, m)
    L2relerr_trunc_vec = zeros(T, m)

    for j in 1:m
        L0 = L0_vec[j]
        println("Running test $j of $m, L0 = $L0")
        N = N_vec[j]
        x = T.((-N ÷ 2):(N ÷ 2)) .* h
        dm = test2demode(x, tdm)

        f = view(f0, (mid - N ÷ 2):(mid + N ÷ 2))
        H_exact = view(H_exact0, (mid - N ÷ 2):(mid + N ÷ 2))

        H_trunc = hilbert(f; dm=dm, pola=NoPolation(), trans=trans)

        dH_trunc = abs.(H_trunc - H_exact)
        err_trunc = maximum(dH_trunc)
        @show err_trunc
        maxerr_trunc_vec[j] = err_trunc

        L2f = sum(f .^ 2) - (f[1]^2 + f[end]^2) / 2
        L2f = sqrt(L2f * h)
        @show L2f

        L2err_trunc = sum(dH_trunc .^ 2) - (dH_trunc[1]^2 + dH_trunc[end]^2) / 2
        L2err_trunc = sqrt(L2err_trunc * h)
        L2relerr_trunc_vec[j] = L2err_trunc / L2f

        @show L2relerr_trunc_vec[j]
    end

    return loss_bench_plot(L0_vec, N_vec, maxerr_trunc_vec, L2relerr_trunc_vec)
end
