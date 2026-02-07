using OpenHilbertBenchmark
using Random, CairoMakie

Random.seed!(42)
T = Float64
func_type = DRationdlFunc(T.([0.5, 1.0, 1.7]), T[])

L0_start = 2^14
L0_rate = 2
test_num = 1
point_density = 2^4

tdm = TestLogLog(; mode_length=200, max_iter=100, is_print=true)
tp = TestPolation(; hann_length=3, herm_length_rate=1)
trans = FIRTrans()

dir = "./logasy/type7-drt"

fig1, fig2 = loss_bench_report(func_type, tdm, tp, trans,
                               L0_start, L0_rate, test_num, point_density,
                               false, dir)
display(fig1)
display(fig2)

save(joinpath(dir, "maxerr.svg"), fig1)
save(joinpath(dir, "L2relerr.svg"), fig2)
