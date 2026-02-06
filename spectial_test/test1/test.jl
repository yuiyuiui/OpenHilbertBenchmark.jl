using OpenHilbertBenchmark
using Random, CairoMakie

testfunc(x, a) = 1 / sqrt(abs(x) + a)

function Htf(x, a; tol=1e-12)
    abs(x) < tol && return 0
    @assert a > 0 "a must be positive"
    res1 = log(sqrt(abs(x)) / (sqrt(abs(x) + a) + sqrt(a))) / sqrt(abs(x) + a)
    if abs(x) > a + tol
        res2 = acos(sqrt(a / abs(x))) / sqrt(abs(x) - a)
    elseif abs(x) < a - tol
        res2 = acosh(sqrt(a / abs(x))) / sqrt(a - abs(x))
    elseif a < abs(x) && abs(x) < a + tol
        res2 = (1 - (abs(x) - a) / a / 3) / sqrt(a)
    else
        res2 = (1 + 5 * (a - abs(x)) / 12 / a) / sqrt(a)
    end
    return (res1 + res2) * 2 * sign(x) / π
end

function Htf(x::Vector, a::Real)
    return [Htf(xi, a) for xi in x]
end

L0_start = 2^4
L0_rate = 3
test_num = 9
point_density = 2^5

tdm = TestLogLog(; mode_length_rate=1//10, is_print=true)
tp = TestPolation(; hann_length=3, herm_length_rate=1)
trans = FIRTrans()

a = 100

fig1, fig2 = loss_bench_report(x -> testfunc(x, a), x -> Htf(x, a), tdm, tp, trans,
                               L0_start, L0_rate, test_num, point_density)
display(fig1)
display(fig2)
