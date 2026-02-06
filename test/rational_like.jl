@testset "DRationdlFunc" begin
    @testset "Construction" begin
        # Test basic construction with Float64
        drf = DRationdlFunc([0.5], Float64[])
        @test drf isa DRationdlFunc{Float64}
        @test drf.even_order_vec isa Vector{Float64}
        @test drf.odd_order_vec isa Vector{Float64}
        @test drf.even_order_vec == [0.5]
        @test drf.odd_order_vec == Float64[]

        # Test construction with Float32
        drf32 = DRationdlFunc([0.5f0], Float32[])
        @test drf32 isa DRationdlFunc{Float32}
        @test drf32.even_order_vec isa Vector{Float32}
        @test drf32.odd_order_vec isa Vector{Float32}
        @test drf32.even_order_vec == [0.5f0]
        @test drf32.odd_order_vec == Float32[]

        # Test boundary values
        drf_min = DRationdlFunc([0.01], Float64[])
        @test drf_min.even_order_vec == [0.01]
        @test drf_min.odd_order_vec == Float64[]

        drf_max = DRationdlFunc([1.0], Float64[])
        @test drf_max.even_order_vec == [1.0]
        @test drf_max.odd_order_vec == Float64[]

        # Test large d values (a > 1/2)
        drf_large = DRationdlFunc([2.0], Float64[])
        @test drf_large.even_order_vec == [2.0]
        @test drf_large.odd_order_vec == Float64[]

        # Test mixed even and odd orders
        drf_mixed = DRationdlFunc([1.0, 2.0], [0.5, 1.5])
        @test length(drf_mixed.even_order_vec) == 2
        @test length(drf_mixed.odd_order_vec) == 2
    end

    @testset "origfunc correctness" begin
        # Test d=0.5 case: f(x) = 1/(1+x²)^0.25
        drf = DRationdlFunc([0.5], Float64[])
        @test origfunc(0.0, drf) ≈ 1.0  # (1+0)^(-0.25) = 1
        @test origfunc(1.0, drf) ≈ 2.0^(-0.25)
        @test origfunc(-1.0, drf) ≈ 2.0^(-0.25)
        @test origfunc(2.0, drf) ≈ 5.0^(-0.25)

        # Test d=1.0 case: f(x) = 1/√(1+x²)
        drf1 = DRationdlFunc([1.0], Float64[])
        @test origfunc(0.0, drf1) ≈ 1.0
        @test origfunc(1.0, drf1) ≈ 1 / sqrt(2.0)
        @test origfunc(-1.0, drf1) ≈ 1 / sqrt(2.0)
        @test origfunc(3.0, drf1) ≈ 1 / sqrt(10.0)

        # Test d=2.0 case (a=1): f(x) = 1/(1+x²)
        drf2 = DRationdlFunc([2.0], Float64[])
        @test origfunc(0.0, drf2) ≈ 1.0
        @test origfunc(1.0, drf2) ≈ 0.5
        @test origfunc(2.0, drf2) ≈ 0.2

        # Test general formula: f(x) = (1+x²)^(-d/2) for single even order
        for d in [0.1, 0.3, 0.7, 0.9, 1.5, 2.0, 3.0]
            drf_test = DRationdlFunc([d], Float64[])
            for x in [-2.0, -1.0, 0.0, 1.0, 2.0]
                expected = (1 + x^2)^(-d / 2)
                @test origfunc(x, drf_test) ≈ expected
            end
        end
    end

    @testset "Hfunc correctness (numerical verification)" begin
        T = Float64
        function hilbert_numerical(x::Real, a::Real; rtol=1e-10)
            f(t) = (1 + t^2)^(-a)

            integrand(u) = (f(x - u) - f(x + u)) / u

            val, err = quadgk(integrand, 0, Inf; rtol=rtol)

            return val / π
        end

        function hilbert_classic(x::Real, a::Real)
            if isapprox(a, 1.0; atol=1e-12)
                return x / (1 + x^2)
            elseif isapprox(a, 0.5; atol=1e-12)
                return (2 / π) * asinh(x) / sqrt(1 + x^2)
            else
                return NaN
            end
        end

        a_vals = [0.5, 1.0, 2.0, 3.5]
        x_vals = [0.2, 1.0, 2.0, 5.0]

        for a in a_vals
            for x in x_vals
                val_num = hilbert_numerical(x, a)
                val_cl = Hfunc(x, DRationdlFunc([2 * a], Float64[]))
                val_ref = hilbert_classic(x, a)

                println("x = $(x) | numerical = $(val_num) | closed = $(val_cl)")

                if !isnan(val_ref)
                    println(" | classic = $(val_ref)")
                end

                # numerical vs closed form
                @test isapprox(val_cl, val_num; rtol=1e-6)

                # when analytic solution is known
                if !isnan(val_ref)
                    @test isapprox(val_cl, val_ref; rtol=1e-10)
                end
            end
        end
    end

    @testset "Type stability" begin
        drf64 = DRationdlFunc([0.5], Float64[])
        drf32 = DRationdlFunc([0.5f0], Float32[])
        drf64_1 = DRationdlFunc([1.0], Float64[])
        drf32_1 = DRationdlFunc([1.0f0], Float32[])
        drf64_2 = DRationdlFunc([2.0], Float64[])  # a > 1/2
        drf32_2 = DRationdlFunc([2.0f0], Float32[])

        # Test origfunc type stability
        @test @inferred(origfunc(1.0, drf64)) isa Float64
        @test @inferred(origfunc(1.0f0, drf32)) isa Float32
        @test @inferred(origfunc(1.0, drf64_1)) isa Float64
        @test @inferred(origfunc(1.0f0, drf32_1)) isa Float32
        @test @inferred(origfunc(1.0, drf64_2)) isa Float64
        @test @inferred(origfunc(1.0f0, drf32_2)) isa Float32

        # Test Hfunc type stability (case: 0 < a < 1/2)
        @test @inferred(Hfunc(1.0, drf64)) isa Float64
        @test @inferred(Hfunc(1.0f0, drf32)) isa Float32

        # Test Hfunc type stability (special case d=1, a=1/2)
        @test @inferred(Hfunc(1.0, drf64_1)) isa Float64
        @test @inferred(Hfunc(1.0f0, drf32_1)) isa Float32

        # Test Hfunc type stability (case: a > 1/2)
        @test @inferred(Hfunc(1.0, drf64_2)) isa Float64
        @test @inferred(Hfunc(1.0f0, drf32_2)) isa Float32
    end

    @testset "Symmetry properties" begin
        # origfunc should be even (symmetric) for single even order
        for d in [0.2, 0.5, 0.8, 1.0, 1.5, 2.0, 3.0]
            drf = DRationdlFunc([d], Float64[])
            for x in [0.5, 1.0, 2.0, 3.0]
                @test origfunc(x, drf) ≈ origfunc(-x, drf)
            end
        end

        # Hfunc should be odd (antisymmetric) for single even order
        for d in [0.2, 0.5, 0.8, 1.0, 1.5, 2.0, 3.0]
            drf = DRationdlFunc([d], Float64[])
            for x in [0.5, 1.0, 2.0, 3.0]
                @test Hfunc(x, drf) ≈ -Hfunc(-x, drf)
            end
        end
    end
end

@testset "LogRationalFunc" begin
    @testset "Construction" begin
        # Test basic construction with Float64
        lrf = LogRationalFunc(Float64; d=1)
        @test lrf.d == 1
        @test lrf isa LogRationalFunc{Float64}

        # Test construction with Float32
        lrf32 = LogRationalFunc(Float32; d=1)
        @test lrf32 isa LogRationalFunc{Float32}
        @test lrf32.d == 1

        # Test default construction (defaults to Float64, d=1)
        lrf_default = LogRationalFunc()
        @test lrf_default isa LogRationalFunc{Float64}
        @test lrf_default.d == 1

        # Test different d values
        for d in [1, 2, 3, 5]
            lrf_test = LogRationalFunc(Float64; d=d)
            @test lrf_test.d == d
        end
    end

    @testset "origfunc correctness" begin
        # f(x) = x² / (x² + 2)^(1 + d/2) / ln(x² + 2)

        # Test d=1 case
        lrf = LogRationalFunc(Float64; d=1)
        # At x=0: f(0) = 0 / (2)^(1.5) / ln(2) = 0
        @test origfunc(0.0, lrf) ≈ 0.0 atol = 1e-15

        # At x=1: f(1) = 1 / 3^(1.5) / ln(3)
        @test origfunc(1.0, lrf) ≈ 1.0 / 3.0^1.5 / log(3.0)
        @test origfunc(-1.0, lrf) ≈ 1.0 / 3.0^1.5 / log(3.0)

        # At x=2: f(2) = 4 / 6^(1.5) / ln(6)
        @test origfunc(2.0, lrf) ≈ 4.0 / 6.0^1.5 / log(6.0)

        # Test d=2 case: f(x) = x² / (x² + 2)^2 / ln(x² + 2)
        lrf2 = LogRationalFunc(Float64; d=2)
        @test origfunc(0.0, lrf2) ≈ 0.0 atol = 1e-15
        @test origfunc(1.0, lrf2) ≈ 1.0 / 9.0 / log(3.0)

        # Test general formula: f(x) = x² / (x² + 2)^(1 + d/2) / ln(x² + 2)
        for d in [1, 2, 3]
            lrf_test = LogRationalFunc(Float64; d=d)
            for x in [-2.0, -1.0, 0.0, 1.0, 2.0]
                expected = x^2 / (x^2 + 2.0)^(1.0 + d / 2.0) / log(x^2 + 2.0)
                @test origfunc(x, lrf_test) ≈ expected
            end
        end
    end

    @testset "Hfunc not implemented" begin
        lrf = LogRationalFunc(Float64; d=1)
        @test_throws ErrorException Hfunc(1.0, lrf)
    end

    @testset "Type stability" begin
        lrf64 = LogRationalFunc(Float64; d=1)
        lrf32 = LogRationalFunc(Float32; d=1)
        lrf64_2 = LogRationalFunc(Float64; d=2)
        lrf32_2 = LogRationalFunc(Float32; d=2)

        # Test origfunc type stability
        @test @inferred(origfunc(1.0, lrf64)) isa Float64
        @test @inferred(origfunc(1.0f0, lrf32)) isa Float32
        @test @inferred(origfunc(1.0, lrf64_2)) isa Float64
        @test @inferred(origfunc(1.0f0, lrf32_2)) isa Float32
    end

    @testset "Symmetry properties" begin
        # origfunc should be even (symmetric): f(x) = f(-x)
        for d in [1, 2, 3]
            lrf = LogRationalFunc(Float64; d=d)
            for x in [0.5, 1.0, 2.0, 3.0]
                @test origfunc(x, lrf) ≈ origfunc(-x, lrf)
            end
        end

        # origfunc should be non-negative for all x
        for d in [1, 2, 3]
            lrf = LogRationalFunc(Float64; d=d)
            for x in [-3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0]
                @test origfunc(x, lrf) >= 0.0
            end
        end
    end

    @testset "Decay behavior" begin
        # f(x) should decay to 0 as |x| → ∞
        lrf = LogRationalFunc(Float64; d=1)
        vals = [origfunc(Float64(x), lrf) for x in [10, 100, 1000]]
        # Each successive value should be smaller
        @test vals[1] > vals[2] > vals[3]
        @test vals[3] < 1e-4
    end
end
