struct MixedFunc{T<:Real}
    swf::Union{SchwartzFunc{T},Nothing}
    rtfpr::Union{RationalFuncPolesRepresent{T},Nothing}
    drtf::Union{DRationdlFunc{T},Nothing}
    logrtf::Union{LogRationalFunc{T},Nothing}
end

function MixedFunc(T::Type{<:Real}; swf::Union{SchwartzFunc{T},Nothing}=nothing,
                   rtfpr::Union{RationalFuncPolesRepresent{T},Nothing}=nothing,
                   drtf::Union{DRationdlFunc{T},Nothing}=nothing,
                   logrtf::Union{LogRationalFunc{T},Nothing}=nothing)
    return MixedFunc{T}(swf, rtfpr, drtf, logrtf)
end

function origfunc(x::T, mf::MixedFunc{T}) where {T<:Real}
    res = T(0)
    if mf.swf !== nothing
        res += origfunc(x, mf.swf)
    end
    if mf.rtfpr !== nothing
        res += origfunc(x, mf.rtfpr)
    end
    if mf.drtf !== nothing
        res += origfunc(x, mf.drtf)
    end
    if mf.logrtf !== nothing
        res += origfunc(x, mf.logrtf)
    end
    return res
end

function Hfunc(x::T, mf::MixedFunc{T}) where {T<:Real}
    res = T(0)
    if mf.swf !== nothing
        res += Hfunc(x, mf.swf)
    end
    if mf.rtfpr !== nothing
        res += Hfunc(x, mf.rtfpr)
    end
    if mf.drtf !== nothing
        res += Hfunc(x, mf.drtf)
    end
    if mf.logrtf !== nothing
        res += Hfunc(x, mf.logrtf)
    end
    return res
end
