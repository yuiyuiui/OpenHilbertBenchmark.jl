include("method.jl")
using .BenchMethod
using CairoMakie

function test_loss_fft(test_num; func_type::String, point_num::Int, p=nothing, L0=2^3,
                       L0_rate=2^3, pv=[], pad_rate=0, k=nothing)
    L0_vec = Int[]
    while test_num > 0
        push!(L0_vec, L0)
        L0 = L0 * L0_rate
        test_num = test_num - 1
    end
    return loss_bench_fft(L0_vec; func_type=func_type,
                          point_num=point_num, p=p, pv=pv, pad_rate=pad_rate, k=k)
end

func_type = "expi_rt"
L0 = 2^2
L0_rate = 2
test_num = 18
p = 1
pv = [1]
pad_rate = 0
point_num = 2^5
k = 1

fig1, fig2 = test_loss_fft(test_num; func_type=func_type, L0=L0, L0_rate=L0_rate,
                           point_num=point_num, p=p, pv=pv, pad_rate=pad_rate, k=k)

display(fig1)
display(fig2)

save("benchreport/maxerr_fft_$(func_type)_pad_rate_$(pad_rate).svg", fig1)
save("benchreport/L2relerr_fft_$(func_type)_pad_rate_$(pad_rate).svg", fig2)
