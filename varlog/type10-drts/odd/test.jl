using OpenHilbertBenchmark
using Random, CairoMakie

Random.seed!(42)
T = Float64
func_type = DRationdlFunc(T[], [0.7, 1.2, 3.8])

L0_start = 2^2
L0_rate = 2
test_num = 16
point_density = 2^4

tdm = TestLogLog(; mode_length=10, max_iter=100)
tp = TestPolation(; hann_length=3, herm_length_rate=1)
trans = FIRTrans()

dir = "./loglog/type10-drts/odd"

fig1, fig2 = loss_bench_report(func_type, tdm, tp, trans,
                               L0_start, L0_rate, test_num, point_density,
                               true, dir)

display(fig1)
display(fig2)

save(joinpath(dir, "maxerr.svg"), fig1)
save(joinpath(dir, "L2relerr.svg"), fig2)
