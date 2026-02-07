using OpenHilbertBenchmark
using Random, CairoMakie

Random.seed!(42)
T = Float64
func_type = RationalFuncPolesRepresent(T)

L0_start = 2^2
L0_rate = 2
test_num = 16
point_density = 12
tdm = TestLogAsy(; mode_length_rate=1//10, degree=12, d=1/6, is_print=true)
tp = TestPolation(; hann_length=3, herm_length_rate=1)
trans = FIRTrans()

dir = "./logasy/type2-rt_1"

fig1, fig2 = loss_bench_report(func_type, tdm, tp, trans,
                               L0_start, L0_rate, test_num, point_density,
                               true, dir)

save(joinpath(dir, "maxerr.svg"), fig1)
save(joinpath(dir, "L2relerr.svg"), fig2)
