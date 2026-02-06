struct SchwartzFunc{T<:Real} <: TestFunc{T}
    A::T
    μ::T
    σ::T
    d::Int

    function SchwartzFunc{T}(A, μ, σ, d::Int) where {T<:Real}
        @assert A > 0 "A must be positive"
        @assert σ > 0 "σ must be positive"
        @assert d == 2 "d must be 2"
        return new{T}(T(A), T(μ), T(σ), d)
    end
end

function SchwartzFunc(T::Type{<:Real}; A::Real=1.0, μ::Real=0.0, σ::Real=1.0, d::Int=2)
    return SchwartzFunc{T}(T(A), T(μ), T(σ), d)
end

function origfunc(x::T, swf::SchwartzFunc{T}) where {T<:Real}
    return exp(-swf.A * abs((x - swf.μ) / swf.σ)^swf.d)
end

function Hfunc(x::T, swf::SchwartzFunc{T}) where {T<:Real}
    return 2 / T(sqrt(π)) * dawson(sqrt(swf.A) * (x - swf.μ) / swf.σ)
end

function Hfunc(x::Vector{T}, swf::SchwartzFunc{T}) where {T<:Real}
    return [Hfunc(xi, swf) for xi in x]
end
