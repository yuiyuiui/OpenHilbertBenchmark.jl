using OpenHilbertBenchmark
using Random, CairoMakie

Random.seed!(42)
T = Float64
func_type = RationalFuncPolesRepresent(T; norder=2, npole=[2, 1])

L0_start = 2^3
L0_rate = 2
test_num = 14
point_density = 2^5
mode_length = 2
dm = AsymptoticDeMode(; degree=5, mode_length=mode_length)

dir = "./asy/type3-rt_2"

fig1, fig2 = loss_bench_report(func_type, dm; trans=FIRTrans(), L0_start=L0_start,
                               L0_rate=L0_rate, test_num=test_num,
                               point_density=point_density, is_saveset=true, file_place=dir)

save(joinpath(dir, "maxerr.svg"), fig1)
save(joinpath(dir, "L2relerr.svg"), fig2)
