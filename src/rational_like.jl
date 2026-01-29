"""
    LogRationalFunc{T<:Real}

Placeholder for log-rational functions.
TODO: Implement specific log-rational function and its Hilbert transform.
"""
struct LogRationalFunc{T<:Real} <: TestFunc{T} end

function LogRationalFunc(T::Type{<:Real})
    return LogRationalFunc{T}()
end

function origfunc(x::T, lrf::LogRationalFunc{T}) where {T<:Real}
    return error("LogRationalFunc origfunc not yet implemented")
end

function Hfunc(x::T, lrf::LogRationalFunc{T}) where {T<:Real}
    return error("LogRationalFunc Hfunc not yet implemented")
end

"""
    DRationdlFunc{T<:Real}

Represents the function f(x) = 1/(1+x²)^a where a = d/2.
Parameter d must be > 0 (so a > 0).

The Hilbert transform is:
- For 0 < a < 1/2: H[f](x) = x * (1+x²)^(-a) * tan(πa)
- For a = 1/2:     H[f](x) = (2/π) * asinh(x) / √(1+x²)
- For a > 1/2:     H[f](x) = x * (1+x²)^(-a) * Γ(a-1/2) / (√π * Γ(a))
"""
struct DRationdlFunc{T<:Real} <: TestFunc{T}
    d::T
    function DRationdlFunc{T}(d::T) where {T<:Real}
        @assert d > 0 "d must be > 0"
        return new{T}(d)
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

"""
    Hfunc(x::T, drf::DRationdlFunc{T}) where {T<:Real}

Compute the Hilbert transform H[1/(1+x²)^a](x) where a = d/2.

- For 0 < a < 1/2: H[f](x) = x * (1+x²)^(-a) * tan(πa)
- For a = 1/2:     H[f](x) = (2/π) * asinh(x) / √(1+x²)
- For a > 1/2:     H[f](x) = x * (1+x²)^(-a) * Γ(a-1/2) / (√π * Γ(a))
"""
function Hfunc(x::T, drf::DRationdlFunc{T}) where {T<:Real}
    a = drf.d / 2
    if isapprox(drf.d, one(T); atol=eps(T)*10)
        # Special case: a = 1/2, d = 1
        # H[1/√(1+x²)](x) = (2/π) * asinh(x) / √(1+x²)
        return T(2/π) * asinh(x) / sqrt(1 + x^2)
    elseif a < T(0.5)
        # Case: 0 < a < 1/2
        # H[(1+x²)^(-a)](x) = x * (1+x²)^(-a) * tan(πa)
        return x * (1 + x^2)^(-a) * tanpi(a)
    else
        # Case: a > 1/2
        # H[(1+x²)^(-a)](x) = x * (1+x²)^(-a) * Γ(a-1/2) / (√π * Γ(a))
        return x * (1 + x^2)^(-a) * gamma(a - T(0.5)) / (sqrt(T(π)) * gamma(a))
    end
end
