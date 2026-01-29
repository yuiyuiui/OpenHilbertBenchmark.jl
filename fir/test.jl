include("method.jl")
using .BenchMethod
using CairoMakie

function test_loss_fir(test_num; func_type::String, extra_output::Int, point_num::Int,
                       p=nothing, k=nothing, L0=2^3, L0_rate=2^3, pv=[])
    L0_vec = Int[]
    while test_num > 0
        push!(L0_vec, L0)
        L0 = L0 * L0_rate
        test_num = test_num - 1
    end
    return loss_bench_fir(L0_vec; func_type=func_type, extra_output=extra_output,
                          point_num=point_num, p=p, k=k, pv=pv)
end

func_type = "mixed"
extra_output = 0
L0 = 2^2
L0_rate = 2
test_num = 16
p = 1
k = 1
pv = [1.5]
point_num = 2^5

fig1, fig2 = test_loss_fir(test_num; func_type=func_type, extra_output=extra_output,
                           L0=L0, L0_rate=L0_rate,
                           point_num=point_num, p=p, k=k, pv=pv)
display(fig1)
display(fig2)

save("benchreport/maxerr_fir_$(func_type).svg", fig1)
save("benchreport/L2relerr_fir_$(func_type).svg", fig2)
