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
drtfunc_odd(x::Real, p::AbstractFloat) = x * drtfunc_even(x, p + 1)

function Hdrtfunc_even(x::T, p::AbstractFloat) where {T<:Real}
    return 2 / beta(T(p / 2), T(1 / 2)) * x *
           pFq((T(1), T(p / 2) + T(0.5)), (T(1.5),), -x^2)
end
function Hdrtfunc_even(mesh::AbstractVector{T}, p::AbstractFloat) where {T<:Real}
    res = zero(mesh)
    pre = 2 / beta(T(p / 2), T(1 / 2))
    @inbounds for i in eachindex(mesh)
        res[i] = pre * mesh[i] * pFq((T(1), T(p / 2) + T(0.5)), (T(1.5),), -mesh[i]^2)
    end
    return res
end
function Hdrtfunc_odd(x::T, p::AbstractFloat) where {T<:Real}
    return Hdrtfunc_even(x, p + 1) * x - beta(T(p / 2), T(1 / 2)) / T(π)
end
function Hdrtfunc_odd(mesh::AbstractVector{T}, p::AbstractFloat) where {T<:Real}
    pre = beta(T(p / 2), T(1 / 2)) / T(π)
    return Hdrtfunc_even(mesh, p + 1) .* mesh .- pre
end

# f(x) = 1/(1+x²)^a where a = d/2
struct DRationdlFunc{T<:Real} <: TestFunc{T}
    even_order_vec::Vector{T}
    odd_order_vec::Vector{T}
end

function origfunc(x::T, drf::DRationdlFunc{T}) where {T<:Real}
    res = T(0)
    for i in 1:length(drf.even_order_vec)
        res += drtfunc_even(x, drf.even_order_vec[i])
    end
    for i in 1:length(drf.odd_order_vec)
        res += drtfunc_odd(x, drf.odd_order_vec[i])
    end
    return res
end

function Hfunc(x::T, drf::DRationdlFunc{T}) where {T<:Real}
    res = Complex{T}(0)
    for i in 1:length(drf.even_order_vec)
        res += Hdrtfunc_even(x, drf.even_order_vec[i])
    end
    for i in 1:length(drf.odd_order_vec)
        res += Hdrtfunc_odd(x, drf.odd_order_vec[i])
    end
    return real(res)
end

function Hfunc(x::Vector{T}, drf::DRationdlFunc{T}) where {T<:Real}
    res = zero(x)
    for i in 1:length(drf.even_order_vec)
        res += Hdrtfunc_even(x, drf.even_order_vec[i])
    end
    for i in 1:length(drf.odd_order_vec)
        res += Hdrtfunc_odd(x, drf.odd_order_vec[i])
    end
    return real(res)
end
