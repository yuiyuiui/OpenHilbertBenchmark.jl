module OpenHilbertBenchmark

using OpenHilbert
using Random, LinearAlgebra, SpecialFunctions, HypergeometricFunctions
using CairoMakie

export origfunc, Hfunc
export SchwartzFunc, RationalFuncPolesRepresent, DRationdlFunc, LogRationalFunc, MixedFunc
export get_funcname, get_algname, write_setting, loss_bench_report

export DeModeMethod, PolationMethod, DiscreteTransMethod, PolationLength
export AsymptoticDeMode, AAADeMode
export NoPolation, InterPolation, ExtraPolation
export FFTTrans, FIRTrans
export hilbert

abstract type TestFunc{T<:Real} end

include("schwartz.jl")
include("rational.jl")
include("rational_like.jl")
include("mix.jl")
include("benchreport.jl")

end
