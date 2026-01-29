@testset "MixedFunc" begin
    @testset "Construction" begin
        T = Float64
        # Test construction with only SchwartzFunc
        swf = SchwartzFunc(T; A=1.0, μ=0.0, σ=1.0, d=2)
        mf1 = MixedFunc(T; swf=swf)
        @test mf1.swf === swf
        @test mf1.rtfpr === nothing
        @test mf1.drtf === nothing
        @test mf1.logrtf === nothing

        # Test construction with only DRationdlFunc
        drf = DRationdlFunc(T; d=0.5)
        mf2 = MixedFunc(T; drtf=drf)
        @test mf2.swf === nothing
        @test mf2.drtf === drf

        # Test construction with multiple components
        Random.seed!(42)
        rtfpr = RationalFuncPolesRepresent(T; norder=1,
                                           poles_scale=1.0,
                                           npole=[2],
                                           amplitudes=[[0.5, 0.5]],
                                           imag_gap=0.5)
        mf3 = MixedFunc(T; swf=swf, rtfpr=rtfpr, drtf=drf)
        @test mf3.swf === swf
        @test mf3.rtfpr === rtfpr
        @test mf3.drtf === drf
        @test mf3.logrtf === nothing
    end

    @testset "origfunc correctness" begin
        # Test that origfunc is additive
        T = Float64
        swf = SchwartzFunc(T; A=1.0, μ=0.0, σ=1.0, d=2)
        drf = DRationdlFunc(T; d=0.5)

        mf = MixedFunc(T; swf=swf, drtf=drf)

        for x in [-2.0, -1.0, 0.0, 1.0, 2.0]
            expected = origfunc(x, swf) + origfunc(x, drf)
            @test origfunc(x, mf) ≈ expected
        end

        # Test with single component
        mf_single = MixedFunc(T; swf=swf)
        for x in [-1.0, 0.0, 1.0]
            @test origfunc(x, mf_single) ≈ origfunc(x, swf)
        end
    end

    @testset "Hfunc correctness" begin
        # Test that Hfunc is additive (linearity of Hilbert transform)
        T = Float64
        swf = SchwartzFunc(T; A=1.0, μ=0.0, σ=1.0, d=2)
        drf = DRationdlFunc(T; d=0.5)

        mf = MixedFunc(T; swf=swf, drtf=drf)

        for x in [-2.0, -1.0, 0.5, 1.0, 2.0]
            expected = Hfunc(x, swf) + Hfunc(x, drf)
            @test Hfunc(x, mf) ≈ expected
        end

        # Test with single component
        mf_single = MixedFunc(T; drtf=drf)
        for x in [-1.0, 0.5, 1.0]
            @test Hfunc(x, mf_single) ≈ Hfunc(x, drf)
        end
    end

    @testset "Type stability" begin
        T1 = Float64
        T2 = Float32
        swf64 = SchwartzFunc(T1; A=1.0, μ=0.0, σ=1.0, d=2)
        swf32 = SchwartzFunc(T2; A=1.0f0, μ=0.0f0, σ=1.0f0, d=2)

        mf64 = MixedFunc(T1; swf=swf64)
        mf32 = MixedFunc(T2; swf=swf32)

        # Test origfunc type stability
        @test @inferred(origfunc(1.0, mf64)) isa T1
        @test @inferred(origfunc(1.0f0, mf32)) isa T2

        # Test Hfunc type stability
        @test @inferred(Hfunc(1.0, mf64)) isa T1
        @test @inferred(Hfunc(1.0f0, mf32)) isa T2
    end

    @testset "Symmetry with mixed components" begin
        # If all components have symmetric origfunc, mixed should too
        T = Float64
        swf = SchwartzFunc(T; A=1.0, μ=0.0, σ=1.0, d=2)  # even
        drf = DRationdlFunc(T; d=0.5)  # even

        mf = MixedFunc(T; swf=swf, drtf=drf)

        # origfunc should be even
        for x in [0.5, 1.0, 2.0]
            @test origfunc(x, mf) ≈ origfunc(-x, mf)
        end

        # Hfunc should be odd
        for x in [0.5, 1.0, 2.0]
            @test Hfunc(x, mf) ≈ -Hfunc(-x, mf)
        end
    end

    @testset "Empty MixedFunc" begin
        # All components are nothing
        T = Float64
        mf_empty = MixedFunc{T}(nothing, nothing, nothing, nothing)

        # Should return 0 for both functions
        @test origfunc(1.0, mf_empty) ≈ T(0)
        @test Hfunc(1.0, mf_empty) ≈ T(0)
    end

    @testset "With RationalFuncPolesRepresent" begin
        Random.seed!(123)
        T = Float64
        rtfpr = RationalFuncPolesRepresent(T; norder=1,
                                           poles_scale=0.5,
                                           npole=[2],
                                           amplitudes=[[0.5, 0.5]],
                                           imag_gap=1.0)
        swf = SchwartzFunc(T; A=1.0, μ=0.0, σ=1.0, d=2)

        mf = MixedFunc(T; swf=swf, rtfpr=rtfpr)

        for x in [-1.0, 0.0, 1.0, 2.0]
            expected_orig = origfunc(x, swf) + origfunc(x, rtfpr)
            expected_H = Hfunc(x, swf) + Hfunc(x, rtfpr)
            @test origfunc(x, mf) ≈ expected_orig
            @test Hfunc(x, mf) ≈ expected_H
        end
    end
end
