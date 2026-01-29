include("method.jl")
using .BenchMethod
using CairoMakie

function test_loss_aaa(test_num; func_type::String, extra_output::Int, point_num::Int,
                       p=nothing, L0=2^3, L0_rate=2^3, pv=[], max_degree::Int=30)
    L0_vec = Int[]
    while test_num > 0
        push!(L0_vec, L0)
        L0 = L0 * L0_rate
        test_num = test_num - 1
    end
    return loss_bench_aaa(L0_vec; func_type=func_type, extra_output=extra_output,
                          point_num=point_num, p=p, pv=pv, max_degree=max_degree)
end

extra_output = 0
L0 = 2^2
L0_rate = 2
test_num = 14
p = 1
max_degree = 30

# obj_vec = [("schwartz", []), ("mixed", [1]), ("mixed", [1, 2, 3, 4, 5])]
obj_vec = [("mixed", [1.5])]

for obj in obj_vec
    func_type, pv = obj

    fig1, fig2 = test_loss_aaa(test_num; func_type=func_type, extra_output=extra_output,
                               L0=L0, L0_rate=L0_rate,
                               point_num=2^5, p=p, pv=pv, max_degree=max_degree)
    display(fig1)
    display(fig2)

    filename = func_type == "mixed" ? "$(func_type)$(length(pv))" : func_type

    save("benchreport/maxerr_aaa_maxnp_$(max_degree)_$(filename).svg", fig1)
    save("benchreport/L2relerr_aaa_maxnp_$(max_degree)_$(filename).svg", fig2)
end
