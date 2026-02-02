using OpenHilbertBenchmark
using Random, CairoMakie

Random.seed!(42)
T = Float64
func_type = MixedFunc(T; swf=SchwartzFunc(T),
                      rtfpr=RationalFuncPolesRepresent(T; norder=2, npole=[2, 1]),
                      drtf=DRationdlFunc(T; d=T(0.5)))

L0_start = 2^2
L0_rate = 2
test_num = 16
point_density = 2^4
mode_length_rate = 0.1
dm = AsymptoticDeMode(mode_length_rate; degree=10, d=0.5)

dir = "./asy/type8-drt+rt+s"

fig1, fig2 = loss_bench_report(func_type, dm; trans=FIRTrans(), L0_start=L0_start,
                               L0_rate=L0_rate, test_num=test_num,
                               point_density=point_density, is_saveset=true,
                               file_place=dir)

save(joinpath(dir, "maxerr.svg"), fig1)
save(joinpath(dir, "L2relerr.svg"), fig2)
