using OpenHilbertBenchmark
using Random, CairoMakie

Random.seed!(42)
T = Float64
func_type = DRationdlFunc([0.5, 0.9, 1.7], [0.7, 1.2, 3.8])

L0_start = 2^2
L0_rate = 2
test_num = 16
point_density = 2^4

tdm = TestNoDeMode()
tp = TestPolation(; hann_length=3, herm_length=3)
trans = FIRTrans()

dir = "./fir/type10-drts"

fig1, fig2 = loss_bench_report(func_type, tdm, tp, trans,
                               L0_start, L0_rate, test_num, point_density,
                               true, dir)

save(joinpath(dir, "maxerr.svg"), fig1)
save(joinpath(dir, "L2relerr.svg"), fig2)
