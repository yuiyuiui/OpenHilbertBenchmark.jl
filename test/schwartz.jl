@testset "SchwartzFunc" begin
    @testset "Construction" begin
        # Test basic construction with Float64
        swf = SchwartzFunc(Float64; A=1.0, μ=0.0, σ=1.0, d=2)
        @test swf.A == 1.0
        @test swf.μ == 0.0
        @test swf.σ == 1.0
        @test swf.d == 2

        # Test construction with Float32
        swf32 = SchwartzFunc(Float32; A=2.0f0, μ=1.0f0, σ=0.5f0, d=1)
        @test swf32.A isa Float32
        @test swf32.μ isa Float32
        @test swf32.σ isa Float32
        @test swf32.d isa Int

        # Test assertion errors
        @test_throws AssertionError SchwartzFunc(Float64; A=-1.0, μ=0.0, σ=1.0, d=2)  # A must be positive
        @test_throws AssertionError SchwartzFunc(Float64; A=1.0, μ=0.0, σ=-1.0, d=2)  # σ must be positive
        @test_throws AssertionError SchwartzFunc(Float64; A=1.0, μ=0.0, σ=1.0, d=3)   # d must be 1 or 2
    end

    @testset "origfunc correctness" begin
        # Test d=2 case: exp(-A * ((x-μ)/σ)^2)
        swf2 = SchwartzFunc(Float64; A=1.0, μ=0.0, σ=1.0, d=2)
        @test origfunc(0.0, swf2) ≈ 1.0  # exp(0) = 1
        @test origfunc(1.0, swf2) ≈ exp(-1.0)
        @test origfunc(-1.0, swf2) ≈ exp(-1.0)
        @test origfunc(2.0, swf2) ≈ exp(-4.0)

        # Test with different parameters
        swf2_shifted = SchwartzFunc(Float64; A=2.0, μ=1.0, σ=0.5, d=2)
        @test origfunc(1.0, swf2_shifted) ≈ 1.0  # at μ
        @test origfunc(1.5, swf2_shifted) ≈ exp(-2.0 * ((1.5 - 1.0) / 0.5)^2)

        # Test d=1 case: exp(-A * |x-μ|/σ)
        swf1 = SchwartzFunc(Float64; A=1.0, μ=0.0, σ=1.0, d=1)
        @test origfunc(0.0, swf1) ≈ 1.0
        @test origfunc(1.0, swf1) ≈ exp(-1.0)
        @test origfunc(-1.0, swf1) ≈ exp(-1.0)
    end

    @testset "Hfunc correctness" begin
        # Test d=2 case: uses Dawson function
        # H[exp(-x^2)] = (2/√π) * dawson(x)
        swf2 = SchwartzFunc(Float64; A=1.0, μ=0.0, σ=1.0, d=2)
        @test Hfunc(0.0, swf2) ≈ 0.0 atol = 1e-14  # dawson(0) = 0
        @test Hfunc(1.0, swf2) ≈ 2 / sqrt(π) * dawson(1.0)
        @test Hfunc(-1.0, swf2) ≈ 2 / sqrt(π) * dawson(-1.0)

        # Test d=1 case: sign(x-μ) * exp(-A/σ * |x-μ|)
        swf1 = SchwartzFunc(Float64; A=1.0, μ=0.0, σ=1.0, d=1)
        @test Hfunc(1.0, swf1) ≈ exp(-1.0)
        @test Hfunc(-1.0, swf1) ≈ -exp(-1.0)
        @test Hfunc(0.0, swf1) ≈ 0.0 atol = 1e-14
    end

    @testset "Type stability" begin
        swf64 = SchwartzFunc(Float64; A=1.0, μ=0.0, σ=1.0, d=2)
        swf32 = SchwartzFunc(Float32; A=1.0f0, μ=0.0f0, σ=1.0f0, d=2)
        swf1_64 = SchwartzFunc(Float64; A=1.0, μ=0.0, σ=1.0, d=1)
        swf1_32 = SchwartzFunc(Float32; A=1.0f0, μ=0.0f0, σ=1.0f0, d=1)

        # Test origfunc type stability
        @test @inferred(origfunc(1.0, swf64)) isa Float64
        @test @inferred(origfunc(1.0f0, swf32)) isa Float32

        # Test Hfunc type stability for d=2
        @test @inferred(Hfunc(1.0, swf64)) isa Float64
        @test @inferred(Hfunc(1.0f0, swf32)) isa Float32

        # Test Hfunc type stability for d=1
        @test @inferred(Hfunc(1.0, swf1_64)) isa Float64
        @test @inferred(Hfunc(1.0f0, swf1_32)) isa Float32
    end

    @testset "Symmetry properties" begin
        swf2 = SchwartzFunc(Float64; A=1.0, μ=0.0, σ=1.0, d=2)
        swf1 = SchwartzFunc(Float64; A=1.0, μ=0.0, σ=1.0, d=1)

        # origfunc should be even (symmetric) for μ=0
        for x in [0.5, 1.0, 2.0, 3.0]
            @test origfunc(x, swf2) ≈ origfunc(-x, swf2)
            @test origfunc(x, swf1) ≈ origfunc(-x, swf1)
        end

        # Hfunc should be odd (antisymmetric) for μ=0
        for x in [0.5, 1.0, 2.0, 3.0]
            @test Hfunc(x, swf2) ≈ -Hfunc(-x, swf2)
            @test Hfunc(x, swf1) ≈ -Hfunc(-x, swf1)
        end
    end
end
