abstract type TestDeMode end

struct TestNoDeMode <: TestDeMode end

function test2demode(grid::Vector{<:Real}, tdm::TestNoDeMode)
    return NoDeMode()
end

struct TestAsy <: TestDeMode
    mode_length::Real
    mode_length_rate::Real
    order_vec::Vector{<:Real}
    is_print::Bool
end

function TestAsy(; degree::Int=2, mode_length_rate::Real=0, mode_length::Real=0,
                 is_print::Bool=false, d::Real=1, order0::Real=0)
    if mode_length > 0 && mode_length_rate > 0
        error("mode_length and mode_length_rate cannot be both positive")
    end
    S = promote_type(typeof(d), typeof(order0))
    order_vec = [S(order0 + k * d) for k in 1:degree]
    return TestAsy(mode_length, mode_length_rate, order_vec, is_print)
end

function test2demode(grid::Vector{<:Real}, testasy::TestAsy)
    degree = length(testasy.order_vec)
    local d, order0
    if degree == 0
        d = 0
        order0 = 0
    elseif degree == 1
        d = testasy.order_vec[1]
        order0 = 0
    else
        d = testasy.order_vec[2] - testasy.order_vec[1]
        order0 = testasy.order_vec[1] - d
    end
    return AsymptoticDeMode(grid; degree=degree,
                            mode_length=testasy.mode_length,
                            mode_length_rate=testasy.mode_length_rate,
                            is_print=testasy.is_print, d=d, order0=order0)
end

struct TestAAA <: TestDeMode
    minsgl::Real
    tol::Real
    max_degree::Int
    lookaheaad::Int
    pcut::Function
    show_poles::Bool
end

function TestAAA(; minsgl::Real=1e-12, tol::Real=1e-12,
                 max_degree::Int=150,
                 lookaheaad::Int=10, pcut::Function=(x -> !isnan(x) && imag(x) != 0),
                 show_poles::Bool=false)
    return TestAAA(minsgl, tol, max_degree, lookaheaad, pcut, show_poles)
end

function test2demode(grid::Vector{<:Real}, testaaa::TestAAA)
    return AAADeMode(grid; minsgl=testaaa.minsgl, tol=testaaa.tol,
                     max_degree=testaaa.max_degree,
                     lookaheaad=testaaa.lookaheaad, pcut=testaaa.pcut,
                     show_poles=testaaa.show_poles)
end

struct TestLogLog <: TestDeMode
    mode_length::Real
    mode_length_rate::Real
    max_iter::Int
    sign_detect_shorten_rate::Real
    is_print::Bool
end

function TestLogLog(; mode_length_rate::Real=0, mode_length::Real=0,
                    max_iter::Int=10, sign_detect_shorten_rate::Real=0.5,
                    is_print::Bool=false)
    if mode_length > 0 && mode_length_rate > 0
        error("mode_length and mode_length_rate cannot be both positive")
    end
    return TestLogLog(mode_length, mode_length_rate, max_iter,
                      sign_detect_shorten_rate, is_print)
end

function test2demode(grid::Vector{<:Real}, testloglog::TestLogLog)
    return LogLogDeMode(grid; mode_length=testloglog.mode_length,
                        mode_length_rate=testloglog.mode_length_rate,
                        max_iter=testloglog.max_iter,
                        sign_detect_shorten_rate=testloglog.sign_detect_shorten_rate,
                        is_print=testloglog.is_print)
end

struct TestLogAsy <: TestDeMode
    mode_length::Real
    mode_length_rate::Real
    order1_scale::Real
    d::Real
    degree::Int
    sign_detect_shorten_rate::Real
    is_print::Bool
end

function TestLogAsy(; mode_length_rate::Real=0, mode_length::Real=0,
                    order1_scale::Real=0, d::Real=1 / 3,
                    degree::Int=9, sign_detect_shorten_rate::Real=0.9, is_print::Bool=false)
    if mode_length > 0 && mode_length_rate > 0
        error("mode_length and mode_length_rate cannot be both positive")
    end
    return TestLogAsy(mode_length, mode_length_rate, order1_scale, d, degree,
                      sign_detect_shorten_rate, is_print)
end

function test2demode(grid::Vector{<:Real}, testlogasy::TestLogAsy)
    return LogAsyDeMode(grid; mode_length=testlogasy.mode_length,
                        mode_length_rate=testlogasy.mode_length_rate,
                        order1_scale=testlogasy.order1_scale, d=testlogasy.d,
                        degree=testlogasy.degree,
                        sign_detect_shorten_rate=testlogasy.sign_detect_shorten_rate,
                        is_print=testlogasy.is_print)
end

struct TestLsqAsy <: TestDeMode
    mode_length::Real
    mode_length_rate::Real
    sign_detect_shorten_rate::Real
    nseek::Int
    lsq_tol::Real
    start_gap::Real
    d::Real
    degree::Int
    is_print::Bool
end

function TestLsqAsy(; mode_length_rate::Real=0, mode_length::Real=0,
                    sign_detect_shorten_rate::Real=0.9, nseek::Int=6,
                    start_gap::Real=1 / 4, d::Real=1 / 3, degree::Int=6,
                    is_print::Bool=false, lsq_tol=10^(-20 // 3))
    if mode_length > 0 && mode_length_rate > 0
        error("mode_length and mode_length_rate cannot be both positive")
    end
    return TestLsqAsy(mode_length, mode_length_rate, sign_detect_shorten_rate, nseek,
                      lsq_tol, start_gap, d, degree, is_print)
end

function test2demode(grid::Vector{<:Real}, testlsqasy::TestLsqAsy)
    return LsqAsyDeMode(grid; mode_length=testlsqasy.mode_length,
                        mode_length_rate=testlsqasy.mode_length_rate,
                        sign_detect_shorten_rate=testlsqasy.sign_detect_shorten_rate,
                        nseek=testlsqasy.nseek,
                        lsq_tol=testlsqasy.lsq_tol,
                        start_gap=testlsqasy.start_gap, d=testlsqasy.d,
                        degree=testlsqasy.degree, is_print=testlsqasy.is_print)
end

struct TestVarLog <: TestDeMode
    start_length_rate::Real
    start_length::Real
    rate::Real
    max_iter::Int
    max_deg::Real
    pad_rate::Int
    sign_detect_shorten_rate::Real
    is_print::Bool
end

function TestVarLog(; start_length_rate::Real=0, start_length::Real=0,
                    rate::Real=3, max_iter::Int=10, max_deg::Real=4, pad_rate::Int=2,
                    sign_detect_shorten_rate::Real=0.9, is_print::Bool=false)
    if start_length > 0 && start_length_rate > 0
        error("start_length and start_length_rate cannot be both positive")
    end
    @assert rate > 1 "rate must be greater than 1"
    @assert sign_detect_shorten_rate > 0 && sign_detect_shorten_rate <= 1 "sign_detect_shorten_rate must be greater than 0 and less than 1"
    return TestVarLog(start_length_rate, start_length, rate, max_iter, max_deg, pad_rate,
                      sign_detect_shorten_rate, is_print)
end

function test2demode(grid::Vector{<:Real}, testvarlog::TestVarLog)
    return VarLogDeMode(grid; start_length_rate=testvarlog.start_length_rate,
                        start_length=testvarlog.start_length,
                        rate=testvarlog.rate, max_iter=testvarlog.max_iter,
                        pad_rate=testvarlog.pad_rate, is_print=testvarlog.is_print,
                        max_deg=testvarlog.max_deg,
                        sign_detect_shorten_rate=testvarlog.sign_detect_shorten_rate)
end

# ======================== TestPolation ========================

struct TestPolation
    hann_length::Real
    hann_length_rate::Real
    herm_length::Real
    herm_length_rate::Real
end

function TestPolation(; hann_length::Real=0, hann_length_rate::Real=0,
                      herm_length::Real=0, herm_length_rate::Real=0)
    if hann_length > 0 && hann_length_rate > 0
        error("hann_length and hann_length_rate cannot be both positive")
    end
    if herm_length > 0 && herm_length_rate > 0
        error("herm_length and herm_length_rate cannot be both positive")
    end
    return TestPolation(hann_length, hann_length_rate, herm_length, herm_length_rate)
end
