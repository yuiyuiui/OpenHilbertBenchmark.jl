using OpenHilbertBenchmark
using Random, CairoMakie

Random.seed!(42)
T = Float64
norder = 3
even_order_vec = 1//2 * rand(T, norder) + 1//2*[1, 2, 3]
odd_order_vec = 1//2 * rand(T, norder) + 1//2*[1, 2, 3]
even_shift_vec = randn(T, norder)
odd_shift_vec = randn(T, norder)
even_σ_vec = rand(T, norder) .+ 1
odd_σ_vec = rand(T, norder) .+ 1
func_type = DRationdlFunc(even_order_vec, odd_order_vec;
                          even_shift_vec=even_shift_vec, odd_shift_vec=odd_shift_vec,
                          even_σ_vec=even_σ_vec, odd_σ_vec=odd_σ_vec)

L0_start = 2^2
L0_rate = 2
test_num = 16
point_density = 8

tdm = TestLogAsy(; mode_length_rate=1//5, degree=12, d=1/6, is_print=true)
tp = TestPolation(; hann_length=3, herm_length_rate=1)
trans = FIRTrans()

dir = "./logasy/type10-drts"

fig1, fig2 = loss_bench_report(func_type, tdm, tp, trans,
                               L0_start, L0_rate, test_num, point_density,
                               true, dir)

display(fig1)
display(fig2)

save(joinpath(dir, "maxerr.svg"), fig1)
save(joinpath(dir, "L2relerr.svg"), fig2)
