struct RationalFuncPolesRepresent{T<:Real} <: TestFunc{T}
    poles::Vector{Vector{Complex{T}}}
    amplitudes::Vector{Vector{T}}
    poles_filter::Function
    Hmask::Vector{Vector{Int}}   # mask for hilbert transform
end

function RationalFuncPolesRepresent(; norder::Int=1, poles_scale::Real=1,
                                    npole::Vector{Int}=[2],
                                    amplitudes::AbstractVector{<:AbstractVector{<:Real}}=[ones(np)
                                                                                          for np in
                                                                                              npole] /
                                                                                         sum(npole),
                                    poles_filter::Function=(x -> !isnan(x) && imag(x) != 0),
                                    imag_gap::Real=1 // 2)
    T = promote_type(eltype(amplitudes[1]), typeof(poles_scale), typeof(imag_gap))
    @assert length(npole) == length(amplitudes) == norder
    @assert norder >= 1
    @assert norm(sum(sum.(amplitudes)) - 1) < eps(real(T))
    poles_vec = Vector{Complex{T}}[]
    Hmask_vec = Vector{Int}[]
    for i in 1:norder
        saved_npole = 0
        poles = Complex{T}[]
        Hmask = Int[]
        while saved_npole < npole[i]
            pole = poles_scale * randn(Complex{T})
            pole += imag_gap * im * ((imag(pole) > 0) ? 1 : -1)

            if poles_filter(pole)
                saved_npole += 1
                push!(poles, pole)
                push!(Hmask, imag(pole) > 0 ? 1 : -1)
            end
        end
        push!(poles_vec, poles)
        push!(Hmask_vec, Hmask)
    end
    amplitudes_converted = [T.(a) for a in amplitudes]
    return RationalFuncPolesRepresent(poles_vec, amplitudes_converted, poles_filter,
                                      Hmask_vec)
end

function origfunc(x::T; rtfpr::RationalFuncPolesRepresent{T}) where {T<:Real}
    res = Complex{T}(0)
    for i in 1:length(rtfpr.poles)
        for j in 1:length(rtfpr.poles[i])
            res += rtfpr.amplitudes[i][j] / (x - rtfpr.poles[i][j])^i
        end
    end
    return real(res)
end

function Hfunc(x::T; rtfpr::RationalFuncPolesRepresent{T}) where {T<:Real}
    res = Complex{T}(0)
    for i in 1:length(rtfpr.poles)
        for j in 1:length(rtfpr.poles[i])
            res += rtfpr.amplitudes[i][j] * rtfpr.Hmask[i][j] * im /
                   (x - rtfpr.poles[i][j])^i
        end
    end
    return real(res)
end
