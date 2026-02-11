using OpenHilbertBenchmark
using Random, CairoMakie

Random.seed!(42)
T = Float64
norder = 3
even_order_vec = 1//2 * rand(T, norder) + 1//2*[1, 2, 3]
#=
 0.8146725615713044
 1.2251694702980969
 1.738703571716409
 =#
#odd_order_vec = 1//2 * rand(T, norder) + 1//2*[1, 2, 3]
#=
 0.8515649245016007
 1.3366730728197482
 1.582947217396567
=#
#=
even_shift_vec = randn(T, norder)
odd_shift_vec = randn(T, norder)
even_σ_vec = rand(T, norder) .+ 1
odd_σ_vec = rand(T, norder) .+ 1
func_type = DRationdlFunc(even_order_vec, odd_order_vec;
                          even_shift_vec=even_shift_vec, odd_shift_vec=odd_shift_vec,
                          even_σ_vec=even_σ_vec, odd_σ_vec=odd_σ_vec)
=#
func_type = DRationdlFunc(even_order_vec, T[])

L0_start = 2^2
L0_rate = 2
test_num = 16
point_density = 16

tdm = TestVarLsq(; start_length=1, rate=4, is_print=false)
trans = FIRTrans()

dir = "./varlsq/type10-drts/varpro"

fig1, fig2 = loss_bench_report(func_type, tdm, trans,
                               L0_start, L0_rate, test_num, point_density,
                               false, dir)

display(fig1)
display(fig2)

save(joinpath(dir, "maxerr.svg"), fig1)
save(joinpath(dir, "L2relerr.svg"), fig2)
