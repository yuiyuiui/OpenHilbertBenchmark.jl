@testset "RationalFuncPolesRepresent" begin
    @testset "Construction" begin
        T = Float64
        Random.seed!(123)
        # Test basic construction
        rtfpr = RationalFuncPolesRepresent(T; norder=1,
                                           poles_scale=1.0,
                                           npole=[3],
                                           amplitudes=[[1 / 3, 1 / 3, 1 / 3]],
                                           imag_gap=0.5)
        @test rtfpr isa RationalFuncPolesRepresent{T}
        @test length(rtfpr.poles) == 1
        @test length(rtfpr.poles[1]) == 3
        @test length(rtfpr.amplitudes) == 1
        @test length(rtfpr.amplitudes[1]) == 3
        @test length(rtfpr.Hmask) == 1
        @test length(rtfpr.Hmask[1]) == 3

        # All poles should have non-zero imaginary part (due to filter)
        for pole in rtfpr.poles[1]
            @test imag(pole) != 0
        end

        # Hmask should be ±1 based on sign of imaginary part
        for (pole, mask) in zip(rtfpr.poles[1], rtfpr.Hmask[1])
            @test mask == (imag(pole) > 0 ? 1 : -1)
        end
    end

    @testset "origfunc correctness" begin
        T = Float64
        Random.seed!(456)
        # Create a simple case with known poles
        # For a single simple pole p with amplitude a:
        # f(x) = a / (x - p)
        # This is complex, but real(f(x)) should match origfunc

        rtfpr = RationalFuncPolesRepresent(T; norder=1,
                                           poles_scale=0.5,
                                           npole=[2],
                                           amplitudes=[[0.5, 0.5]],
                                           imag_gap=1.0)
        @test rtfpr isa RationalFuncPolesRepresent{T}
        # Test that origfunc returns real values
        for x in [-2.0, -1.0, 0.0, 1.0, 2.0]
            result = origfunc(x, rtfpr)
            @test result isa Real
            @test isfinite(result)
        end

        # Verify the formula manually
        x = 1.5
        expected = 0.0 + 0.0im
        for i in 1:length(rtfpr.poles)
            for j in 1:length(rtfpr.poles[i])
                expected += rtfpr.amplitudes[i][j] / (x - rtfpr.poles[i][j])^i
            end
        end
        @test origfunc(x, rtfpr) ≈ real(expected)
    end

    @testset "Hfunc correctness" begin
        T = Float64
        Random.seed!(789)
        rtfpr = RationalFuncPolesRepresent(T; norder=1,
                                           poles_scale=0.5,
                                           npole=[2],
                                           amplitudes=[[0.5, 0.5]],
                                           imag_gap=1.0)
        @test rtfpr isa RationalFuncPolesRepresent{T}
        # Test that Hfunc returns real values
        for x in [-2.0, -1.0, 0.0, 1.0, 2.0]
            result = Hfunc(x, rtfpr)
            @test result isa Real
            @test isfinite(result)
        end

        # Verify the formula manually
        x = 1.5
        expected = 0.0 + 0.0im
        for i in 1:length(rtfpr.poles)
            for j in 1:length(rtfpr.poles[i])
                expected += rtfpr.amplitudes[i][j] * rtfpr.Hmask[i][j] * im /
                            (x - rtfpr.poles[i][j])^i
            end
        end
        @test Hfunc(x, rtfpr) ≈ real(expected)
    end

    @testset "Type stability" begin
        T = Float64
        Random.seed!(101)
        rtfpr = RationalFuncPolesRepresent(T; norder=1,
                                           poles_scale=1.0,
                                           npole=[2],
                                           amplitudes=[[0.5, 0.5]],
                                           imag_gap=0.5)
        @test rtfpr isa RationalFuncPolesRepresent{T}
        # Test origfunc type stability
        @test @inferred(origfunc(1.0, rtfpr)) isa Float64

        # Test Hfunc type stability
        @test @inferred(Hfunc(1.0, rtfpr)) isa Float64
    end

    @testset "Multi-order construction" begin
        T = Float64
        Random.seed!(202)
        # Test construction with multiple orders
        rtfpr = RationalFuncPolesRepresent(T; norder=2,
                                           poles_scale=1.0,
                                           npole=[2, 3],
                                           amplitudes=[[0.2, 0.2], [0.2, 0.2, 0.2]],
                                           imag_gap=0.5)
        @test rtfpr isa RationalFuncPolesRepresent{T}
        @test length(rtfpr.poles) == 2
        @test length(rtfpr.poles[1]) == 2
        @test length(rtfpr.poles[2]) == 3
        @test length(rtfpr.amplitudes) == 2
        @test length(rtfpr.amplitudes[1]) == 2
        @test length(rtfpr.amplitudes[2]) == 3
    end

    @testset "Amplitude normalization" begin
        # Test that amplitudes must sum to 1
        T = Float64
        @test_throws AssertionError RationalFuncPolesRepresent(T; norder=1,
                                                               poles_scale=1.0,
                                                               npole=[2],
                                                               amplitudes=[[0.3, 0.3]],  # sum = 0.6 ≠ 1
                                                               imag_gap=0.5)
    end
end
