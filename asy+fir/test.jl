include("method.jl")
using .BenchMethod
using CairoMakie

function test_loss_asy(test_num; func_type::String, extra_output::Int, point_num::Int,
                       p=nothing, L0=2^3, L0_rate=2^3, pv=[], degree::Int=2,
                       mode_length=2)
    L0_vec = Int[]
    while test_num > 0
        push!(L0_vec, L0)
        L0 = L0 * L0_rate
        test_num = test_num - 1
    end
    return loss_bench_asy(L0_vec; func_type=func_type, extra_output=extra_output,
                          point_num=point_num, p=p, pv=pv, degree=degree,
                          mode_length=mode_length)
end

func_type = "mixed"
extra_output = 0
L0 = 2^3
L0_rate = 2
test_num = 14
p = 1
pv = [1.5]
degree = 10
mode_length = 1
point_num = 2^5

fig1, fig2 = test_loss_asy(test_num; func_type=func_type, extra_output=extra_output,
                           L0=L0, L0_rate=L0_rate,
                           point_num=point_num, p=p, pv=pv, degree=degree,
                           mode_length=mode_length)
display(fig1)
display(fig2)

filename = func_type == "mixed" ? "$(func_type)$(length(pv))" : func_type

save("benchreport/maxerr_asy_$(filename)_$(degree).svg", fig1)
save("benchreport/L2relerr_asy_$(filename)_$(degree).svg", fig2)
