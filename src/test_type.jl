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
    N = length(grid)
    if testasy.mode_length_rate > 0
        n = max(1, round(Int, N * testasy.mode_length_rate))
    elseif testasy.mode_length > 0
        n = max(1, round(Int, testasy.mode_length))
    else
        n = max(1, N ÷ 20)
    end
    return AsymptoticDeMode(n; ndeg=testasy.degree, is_print=testasy.is_print, d=testasy.d,
                            order1=testasy.order1, is_rational=testasy.is_rational)
end
