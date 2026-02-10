using OpenHilbertBenchmark
using Random, CairoMakie

Random.seed!(42)
T = Float64
#=
even_shift_vec = T.([1.0, 2.0, 3.0])
even_σ_vec = T.([4.0, 5.0, 6.0])
func_type = DRationdlFunc(T.([0.5, 0.9, 1.7]), T[];
                          even_shift_vec=even_shift_vec, even_σ_vec=even_σ_vec)
=#
func_type = DRationdlFunc(T.([0.5, 0.9, 1.7]), T[])

L0_start = 2^2
L0_rate = 2
test_num = 16
point_density = 2^4

tdm = TestLogAsy(; mode_length_rate=1 // 5, is_print=true, d=1 / 4, degree=12)
tp = TestPolation(; hann_length=3, herm_length_rate=1)
trans = FIRTrans()

dir = "./logasy/type10-drts/even"

fig1, fig2 = loss_bench_report(func_type, tdm, tp, trans,
                               L0_start, L0_rate, test_num, point_density,
                               true, dir)

display(fig1)
display(fig2)

save(joinpath(dir, "maxerr.svg"), fig1)
save(joinpath(dir, "L2relerr.svg"), fig2)
