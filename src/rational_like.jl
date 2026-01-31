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

# f(x) = 1/(1+x²)^a where a = d/2
struct DRationdlFunc{T<:Real} <: TestFunc{T}
    d::T
    prefactor::T
    function DRationdlFunc{T}(d::T) where {T<:Real}
        @assert d > 0 "d must be > 0"
        return new{T}(d, T(2 / beta(d / 2, T(1 / 2))))
    end
end

function DRationdlFunc(T::Type{<:Real}; d::Real=0.5)
    return DRationdlFunc{T}(T(d))
end

"""
    origfunc(x::T, drf::DRationdlFunc{T}) where {T<:Real}

Compute f(x) = 1/(1+x²)^a where a = d/2.
"""
function origfunc(x::T, drf::DRationdlFunc{T}) where {T<:Real}
    a = drf.d / 2
    return (1 + x^2)^(-a)
end

function Hfunc(x::T, drf::DRationdlFunc{T}) where {T<:Real}
    a = drf.d / 2
    return drf.prefactor * x * pFq((T(1), a + T(0.5)), (T(1.5),), -x^2)
end
