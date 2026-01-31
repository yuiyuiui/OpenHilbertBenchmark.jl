using OpenHilbertBenchmark
using Random, CairoMakie

Random.seed!(42)
T = Float64
func_type = DRationdlFunc(T; d=T(1.5))

L0_start = 2^2
L0_rate = 2
test_num = 15
point_density = 2^4
mode_length = 1
dm = AsymptoticDeMode(; degree=5, mode_length=mode_length)

dir = "./asy+fir/1.5_rational"

fig1, fig2 = loss_bench_report(func_type, dm; trans=FIRTrans(), L0_start=L0_start,
                               L0_rate=L0_rate, test_num=test_num,
                               point_density=point_density, is_saveset=true,
                               file_place=dir)

save(joinpath(dir, "maxerr.svg"), fig1)
save(joinpath(dir, "L2relerr.svg"), fig2)
