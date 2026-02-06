# d  means f(x) = x²/(x²+2)^{1+d/2}/ln(x²+2)
struct LogRationalFunc{T<:Real} <: TestFunc{T}
    d::Int
end

function LogRationalFunc(T::Type{<:Real}=Float64; d::Int=1)
    return LogRationalFunc{T}(d)
end

function origfunc(x::T, lrf::LogRationalFunc{T}) where {T<:Real}
    return x^T(2) / (x^T(2) + T(2))^(T(1) + lrf.d / T(2)) / log(x^T(2) + T(2))
end

function Hfunc(x::T, lrf::LogRationalFunc{T})::T where {T<:Real}
    return error("Hfunc for LogRationalFunc is not yet implemented")
end

# ===============================================

drtfunc_even(x::T, p::AbstractFloat) where {T<:Real} = 1 / (x^2 + 1)^(T(p) / 2)
function drtfunc_even(x::T, p::AbstractFloat, shift::Real, σ::Real) where {T<:Real}
    return drtfunc_even((x-T(shift))/T(σ), p) / T(σ)^p
end
drtfunc_odd(x::Real, p::AbstractFloat) = x * drtfunc_even(x, p + 1)
function drtfunc_odd(x::T, p::AbstractFloat, shift::Real, σ::Real) where {T<:Real}
    return drtfunc_odd((x-T(shift))/T(σ), p) / T(σ)^p
end

function Hdrtfunc_even(mesh::AbstractVector{T}, p::AbstractFloat) where {T<:Real}
    res = zero(mesh)
    pre = 2 / beta(T(p / 2), T(1 / 2))
    @inbounds for i in eachindex(mesh)
        res[i] = pre * mesh[i] * pFq((T(1), T(p / 2) + T(0.5)), (T(1.5),), -mesh[i]^2)
    end
    return res
end
function Hdrtfunc_even(mesh::AbstractVector{T}, p::AbstractFloat, shift::Real,
                       σ::Real) where {T<:Real}
    return Hdrtfunc_even((mesh .- T(shift)) ./ T(σ), p) ./ T(σ)^p
end

function Hdrtfunc_odd(mesh::AbstractVector{T}, p::AbstractFloat) where {T<:Real}
    pre = beta(T(p / 2), T(1 / 2)) / T(π)
    return Hdrtfunc_even(mesh, p + 1) .* mesh .- pre
end
function Hdrtfunc_odd(mesh::AbstractVector{T}, p::AbstractFloat, shift::Real,
                      σ::Real) where {T<:Real}
    return Hdrtfunc_odd((mesh .- T(shift)) ./ T(σ), p) ./ T(σ)^p
end

# f(x) = 1/(1+x²)^a where a = d/2
struct DRationdlFunc{T<:Real} <: TestFunc{T}
    even_order_vec::Vector{T}
    even_shift_vec::Vector{T}
    even_σ_vec::Vector{T}
    odd_order_vec::Vector{T}
    odd_shift_vec::Vector{T}
    odd_σ_vec::Vector{T}
end

function DRationdlFunc(even_order_vec::Vector{T}, odd_order_vec::Vector{T};
                       even_shift_vec=T[], even_σ_vec=T[], odd_shift_vec=T[],
                       odd_σ_vec=T[]) where {T<:Real}
    if length(even_shift_vec) == 0 && length(even_σ_vec) == 0
        even_shift_vec = zeros(T, length(even_order_vec))
        even_σ_vec = ones(T, length(even_order_vec))
    end

    if length(odd_shift_vec) == 0 && length(odd_σ_vec) == 0
        odd_shift_vec = zeros(T, length(odd_order_vec))
        odd_σ_vec = ones(T, length(odd_order_vec))
    end

    @assert length(even_order_vec) == length(even_shift_vec) == length(even_σ_vec)
    @assert length(odd_order_vec) == length(odd_shift_vec) == length(odd_σ_vec)
    return DRationdlFunc{T}(even_order_vec, even_shift_vec, even_σ_vec, odd_order_vec,
                            odd_shift_vec, odd_σ_vec)
end

function origfunc(x::T, drf::DRationdlFunc{T}) where {T<:Real}
    res = T(0)
    for i in 1:length(drf.even_order_vec)
        res += drtfunc_even(x, drf.even_order_vec[i], drf.even_shift_vec[i],
                            drf.even_σ_vec[i])
    end
    for i in 1:length(drf.odd_order_vec)
        res += drtfunc_odd(x, drf.odd_order_vec[i], drf.odd_shift_vec[i], drf.odd_σ_vec[i])
    end
    return res
end

function Hfunc(x::T, drf::DRationdlFunc{T}) where {T<:Real}
    res = Complex{T}(0)
    for i in 1:length(drf.even_order_vec)
        res += Hdrtfunc_even(x, drf.even_order_vec[i], drf.even_shift_vec[i],
                             drf.even_σ_vec[i])
    end
    for i in 1:length(drf.odd_order_vec)
        res += Hdrtfunc_odd(x, drf.odd_order_vec[i], drf.odd_shift_vec[i], drf.odd_σ_vec[i])
    end
    return real(res)
end

function Hfunc(x::Vector{T}, drf::DRationdlFunc{T}) where {T<:Real}
    res = zero(x)
    for i in 1:length(drf.even_order_vec)
        res += Hdrtfunc_even(x, drf.even_order_vec[i], drf.even_shift_vec[i],
                             drf.even_σ_vec[i])
    end
    for i in 1:length(drf.odd_order_vec)
        res += Hdrtfunc_odd(x, drf.odd_order_vec[i], drf.odd_shift_vec[i], drf.odd_σ_vec[i])
    end
    return real(res)
end
