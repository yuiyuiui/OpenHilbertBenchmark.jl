using OpenHilbertBenchmark
using Random, CairoMakie

Random.seed!(42)
T = Float64
func_type = SchwartzFunc(T)

L0_start = 2^2
L0_rate = 2
test_num = 14
point_density = 2^5
dm = AAADeMode(; max_degree=20)

dir = "./aaa+fir/schwartz"

fig1, fig2 = loss_bench_report(func_type, dm; trans=FIRTrans(), L0_start=L0_start,
                               L0_rate=L0_rate, test_num=test_num,
                               point_density=point_density, is_saveset=true, file_place=dir)

display(fig1)
display(fig2)

save(joinpath(dir, "maxerr.svg"), fig1)
save(joinpath(dir, "L2relerr.svg"), fig2)
