using Test, LinearAlgebra, Random
using SpecialFunctions
using OpenHilbertBenchmark

Random.seed!(42)

include("schwartz.jl")
include("rational.jl")
#include("rational_like.jl")
include("mix.jl")
include("benchreport.jl")
