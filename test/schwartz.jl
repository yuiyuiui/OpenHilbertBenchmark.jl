@testset "SchwartzFunc" begin
    @testset "Construction" begin
        # Test basic construction with Float64
        swf = SchwartzFunc(Float64; A=1.0, μ=0.0, σ=1.0, d=2)
        @test swf.A == 1.0
        @test swf.μ == 0.0
        @test swf.σ == 1.0
        @test swf.d == 2

        # Test construction with Float32
        swf32 = SchwartzFunc(Float32; A=2.0f0, μ=0.0f0, σ=1.0f0, d=2)
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
        swf = SchwartzFunc(Float64; A=1.0, μ=0.0, σ=1.0, d=2)
        @test origfunc(0.0, swf) ≈ 1.0  # exp(0) = 1
        @test origfunc(1.0, swf) ≈ exp(-1.0)
        @test origfunc(-1.0, swf) ≈ exp(-1.0)
        @test origfunc(2.0, swf) ≈ exp(-4.0)

        # Test with different parameters
        swf_shifted = SchwartzFunc(Float64; A=2.0, μ=1.0, σ=0.5, d=2)
        @test origfunc(1.0, swf_shifted) ≈ 1.0  # at μ
        @test origfunc(1.5, swf_shifted) ≈ exp(-2.0 * ((1.5 - 1.0) / 0.5)^2)
    end

    @testset "Hfunc correctness" begin
        # Test d=2 case: uses Dawson function
        # H[exp(-x^2)] = (2/√π) * dawson(x)
        swf = SchwartzFunc(Float64; A=1.0, μ=0.0, σ=1.0, d=2)
        @test Hfunc(0.0, swf) ≈ 0.0 atol = 1e-14  # dawson(0) = 0
        @test Hfunc(1.0, swf) ≈ 2 / sqrt(π) * dawson(1.0)
        @test Hfunc(-1.0, swf) ≈ 2 / sqrt(π) * dawson(-1.0)
    end

    @testset "Type stability" begin
        swf64 = SchwartzFunc(Float64; A=1.0, μ=0.0, σ=1.0, d=2)
        swf32 = SchwartzFunc(Float32; A=1.0f0, μ=0.0f0, σ=1.0f0, d=2)

        # Test origfunc type stability
        @test @inferred(origfunc(1.0, swf64)) isa Float64
        @test @inferred(origfunc(1.0f0, swf32)) isa Float32

        # Test Hfunc type stability for d=2
        @test @inferred(Hfunc(1.0, swf64)) isa Float64
        @test @inferred(Hfunc(1.0f0, swf32)) isa Float32
    end

    @testset "Symmetry properties" begin
        swf = SchwartzFunc(Float64; A=1.0, μ=0.0, σ=1.0, d=2)

        # origfunc should be even (symmetric) for μ=0
        for x in [0.5, 1.0, 2.0, 3.0]
            @test origfunc(x, swf) ≈ origfunc(-x, swf)
        end

        # Hfunc should be odd (antisymmetric) for μ=0
        for x in [0.5, 1.0, 2.0, 3.0]
            @test Hfunc(x, swf) ≈ -Hfunc(-x, swf)
        end
    end
end
