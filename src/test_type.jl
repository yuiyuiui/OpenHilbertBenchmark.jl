abstract type TestDeMode end

struct TestNoDeMode <: TestDeMode end

function test2demode(grid::Vector{<:Real}, tdm::TestNoDeMode)
    return NoDeMode()
end

struct TestAsy <: TestDeMode
    mode_length::Real
    mode_length_rate::Real
    d::Real
    order1::Real
    degree::Int
    order_vec::Vector{<:Real}
    is_print::Bool
    is_rational::Bool
end

function TestAsy(; degree::Int=2, mode_length_rate::Real=0, mode_length::Real=0,
                 is_print::Bool=false, d::Real=1, order1::Real=1, is_rational::Bool=false)
    if mode_length > 0 && mode_length_rate > 0
        error("mode_length and mode_length_rate cannot be both positive")
    end
    S = promote_type(typeof(d), typeof(order1))
    order_vec = [S(order1 + k * d) for k in 0:(degree - 1)]
    return TestAsy(mode_length, mode_length_rate, d, order1, degree, order_vec, is_print,
                   is_rational)
end

function test2demode(grid::Vector{<:Real}, testasy::TestAsy)
    return AsymptoticDeMode(grid; degree=testasy.degree,
                            mode_length=testasy.mode_length,
                            mode_length_rate=testasy.mode_length_rate,
                            is_print=testasy.is_print, d=testasy.d, order1=testasy.order1,
                            is_rational=testasy.is_rational)
end

struct TestVarLsq <: TestDeMode
    start_length_rate::Real
    start_length::Real
    start_n::Int

    rate::Real
    nseek_vec::Vector{Int}
    is_varpro::Bool

    max_deg::Real
    lsq_tol::Real
    fine_rate::Real
    mindeg_gap::Real

    sign_rate::Real
    is_print::Bool
end

function TestVarLsq(; start_length_rate::Real=0, start_length::Real=0,
                    start_n::Int=1, rate::Real=3, nseek_vec::Vector{Int}=[1, 3, 6, 9, 12],
                    is_varpro::Bool=true, max_deg::Real=4, lsq_tol::Real=1e-12,
                    fine_rate::Real=1 // 6, mindeg_gap::Real=1 // 6,
                    sign_rate::Real=0.9, is_print::Bool=false)
    return TestVarLsq(start_length_rate, start_length, start_n, rate, nseek_vec, is_varpro,
                      max_deg, lsq_tol, fine_rate, mindeg_gap, sign_rate,
                      is_print)
end

function test2demode(grid::Vector{<:Real}, testvarlsq::TestVarLsq)
    return VarLsqDeMode(grid; start_length_rate=testvarlsq.start_length_rate,
                        start_length=testvarlsq.start_length,
                        rate=testvarlsq.rate, nseek_vec=testvarlsq.nseek_vec,
                        max_deg=testvarlsq.max_deg, lsq_tol=testvarlsq.lsq_tol,
                        is_varpro=testvarlsq.is_varpro,
                        fine_rate=testvarlsq.fine_rate,
                        mindeg_gap=testvarlsq.mindeg_gap,
                        sign_rate=testvarlsq.sign_rate,
                        is_print=testvarlsq.is_print)
end
