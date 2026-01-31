@testset "DRationdlFunc" begin
    @testset "Construction" begin
        # Test basic construction with Float64
        drf = DRationdlFunc(Float64; d=0.5)
        @test drf.d == 0.5

        # Test construction with Float32
        drf32 = DRationdlFunc(Float32; d=0.5f0)
        @test drf32.d isa Float32
        @test drf32.d == 0.5f0

        # Test boundary values
        drf_min = DRationdlFunc(Float64; d=0.01)
        @test drf_min.d == 0.01

        drf_max = DRationdlFunc(Float64; d=1.0)
        @test drf_max.d == 1.0

        # Test large d values (a > 1/2)
        drf_large = DRationdlFunc(Float64; d=2.0)
        @test drf_large.d == 2.0

        drf_larger = DRationdlFunc(Float64; d=4.0)
        @test drf_larger.d == 4.0

        # Test assertion errors
        @test_throws AssertionError DRationdlFunc(Float64; d=0.0)   # d must be > 0
        @test_throws AssertionError DRationdlFunc(Float64; d=-0.5)  # d must be > 0
    end

    @testset "origfunc correctness" begin
        # Test d=0.5 case: f(x) = 1/(1+x²)^0.25
        drf = DRationdlFunc(Float64; d=0.5)
        @test origfunc(0.0, drf) ≈ 1.0  # (1+0)^(-0.25) = 1
        @test origfunc(1.0, drf) ≈ 2.0^(-0.25)
        @test origfunc(-1.0, drf) ≈ 2.0^(-0.25)
        @test origfunc(2.0, drf) ≈ 5.0^(-0.25)

        # Test d=1.0 case: f(x) = 1/√(1+x²)
        drf1 = DRationdlFunc(Float64; d=1.0)
        @test origfunc(0.0, drf1) ≈ 1.0
        @test origfunc(1.0, drf1) ≈ 1 / sqrt(2.0)
        @test origfunc(-1.0, drf1) ≈ 1 / sqrt(2.0)
        @test origfunc(3.0, drf1) ≈ 1 / sqrt(10.0)

        # Test d=2.0 case (a=1): f(x) = 1/(1+x²)
        drf2 = DRationdlFunc(Float64; d=2.0)
        @test origfunc(0.0, drf2) ≈ 1.0
        @test origfunc(1.0, drf2) ≈ 0.5
        @test origfunc(2.0, drf2) ≈ 0.2

        # Test general formula: f(x) = (1+x²)^(-d/2)
        for d in [0.1, 0.3, 0.7, 0.9, 1.5, 2.0, 3.0]
            drf_test = DRationdlFunc(Float64; d=d)
            for x in [-2.0, -1.0, 0.0, 1.0, 2.0]
                expected = (1 + x^2)^(-d / 2)
                @test origfunc(x, drf_test) ≈ expected
            end
        end
    end

    @testset "Hfunc correctness" begin
        # Test d=1.0 case (a=1/2): H[1/√(1+x²)](x) = (2/π) * asinh(x) / √(1+x²)
        drf1 = DRationdlFunc(Float64; d=1.0)
        @test Hfunc(0.0, drf1) ≈ 0.0 atol = 1e-14  # asinh(0) = 0
        @test Hfunc(1.0, drf1) ≈ (2 / π) * asinh(1.0) / sqrt(2.0)
        @test Hfunc(-1.0, drf1) ≈ (2 / π) * asinh(-1.0) / sqrt(2.0)
        @test Hfunc(2.0, drf1) ≈ (2 / π) * asinh(2.0) / sqrt(5.0)

        # Test d=2.0 case (a=1): H[1/(1+x²)](x) = x/(1+x²)
        # This is a well-known result!
        drf2 = DRationdlFunc(Float64; d=2.0)
        @test Hfunc(0.0, drf2) ≈ 0.0 atol = 1e-14
        @test Hfunc(1.0, drf2) ≈ 1.0 / 2.0
        @test Hfunc(-1.0, drf2) ≈ -1.0 / 2.0
        @test Hfunc(2.0, drf2) ≈ 2.0 / 5.0
        @test Hfunc(3.0, drf2) ≈ 3.0 / 10.0

        # Test case: 0 < a < 1/2
        # H[(1+x²)^(-a)](x) = x * (1+x²)^(-a) * tan(πa)
        for d in [0.1, 0.3, 0.5, 0.7, 0.9]
            drf_test = DRationdlFunc(Float64; d=d)
            a = d / 2
            for x in [-2.0, -1.0, 0.5, 1.0, 2.0]
                expected = x * (1 + x^2)^(-a) * tanpi(a)
                @test Hfunc(x, drf_test) ≈ expected rtol = 1e-10
            end
        end

        # Test case: a > 1/2
        # H[(1+x²)^(-a)](x) = x * (1+x²)^(-a) * Γ(a-1/2) / (√π * Γ(a))
        for d in [1.5, 2.0, 3.0, 4.0]
            drf_test = DRationdlFunc(Float64; d=d)
            a = d / 2
            for x in [-2.0, -1.0, 0.5, 1.0, 2.0]
                expected = x * (1 + x^2)^(-a) * gamma(a - 0.5) / (sqrt(π) * gamma(a))
                @test Hfunc(x, drf_test) ≈ expected rtol = 1e-10
            end
        end

        # Hfunc(0, drf) should be 0 for any d (odd function)
        for d in [0.1, 0.5, 0.9, 1.0, 1.5, 2.0, 3.0]
            drf_test = DRationdlFunc(Float64; d=d)
            @test Hfunc(0.0, drf_test) ≈ 0.0 atol = 1e-14
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
                val_cl = Hfunc(x, DRationdlFunc(T; d=2 * a))
                val_ref = hilbert_classic(x, a)

                @printf("x = %-5.2f | numerical = %+12.8e | closed = %+12.8e",
                        x, val_num, val_cl)

                if !isnan(val_ref)
                    @printf(" | classic = %+12.8e", val_ref)
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
        drf64 = DRationdlFunc(Float64; d=0.5)
        drf32 = DRationdlFunc(Float32; d=0.5f0)
        drf64_1 = DRationdlFunc(Float64; d=1.0)
        drf32_1 = DRationdlFunc(Float32; d=1.0f0)
        drf64_2 = DRationdlFunc(Float64; d=2.0)  # a > 1/2
        drf32_2 = DRationdlFunc(Float32; d=2.0f0)

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
        # origfunc should be even (symmetric) for all d
        for d in [0.2, 0.5, 0.8, 1.0, 1.5, 2.0, 3.0]
            drf = DRationdlFunc(Float64; d=d)
            for x in [0.5, 1.0, 2.0, 3.0]
                @test origfunc(x, drf) ≈ origfunc(-x, drf)
            end
        end

        # Hfunc should be odd (antisymmetric) for all d
        for d in [0.2, 0.5, 0.8, 1.0, 1.5, 2.0, 3.0]
            drf = DRationdlFunc(Float64; d=d)
            for x in [0.5, 1.0, 2.0, 3.0]
                @test Hfunc(x, drf) ≈ -Hfunc(-x, drf)
            end
        end
    end

    @testset "Asymptotic behavior" begin
        # For large |x|, origfunc should decay as |x|^(-d)
        x_large = 100.0

        for d in [0.5, 1.0, 2.0, 3.0]
            drf = DRationdlFunc(Float64; d=d)
            # (1 + x²)^(-d/2) ≈ x^(-d) for large x
            @test origfunc(x_large, drf) ≈ x_large^(-d) rtol = 0.01
        end

        # For large |x| with a < 1/2, Hfunc ≈ x^(1-d) * tan(πa)
        drf = DRationdlFunc(Float64; d=0.5)
        a = 0.25
        @test Hfunc(x_large, drf) ≈ x_large^(1 - 0.5) * tanpi(a) rtol = 0.01

        # For large |x| with a > 1/2, Hfunc ≈ x^(1-d) * Γ(a-1/2)/(√π*Γ(a))
        drf2 = DRationdlFunc(Float64; d=2.0)
        a2 = 1.0
        coef = gamma(a2 - 0.5) / (sqrt(π) * gamma(a2))
        @test Hfunc(x_large, drf2) ≈ x_large^(1 - 2.0) * coef rtol = 0.01
    end

    @testset "Consistency checks" begin
        # Check that Hfunc satisfies expected properties

        # For small x, Hfunc ≈ x * coefficient (to leading order)
        drf = DRationdlFunc(Float64; d=0.5)
        a = 0.25
        x_small = 0.01
        @test Hfunc(x_small, drf) ≈ x_small * tanpi(a) rtol = 0.01

        # For a < 1/2: ratio Hfunc(x)/origfunc(x) = x*tan(πa)
        for d in [0.2, 0.4, 0.6, 0.8]
            drf_test = DRationdlFunc(Float64; d=d)
            a_test = d / 2
            for x in [0.5, 1.0, 2.0]
                ratio = Hfunc(x, drf_test) / origfunc(x, drf_test)
                expected_ratio = x * tanpi(a_test)
                @test ratio ≈ expected_ratio rtol = 1e-10
            end
        end

        # For a > 1/2: ratio Hfunc(x)/origfunc(x) = x*Γ(a-1/2)/(√π*Γ(a))
        for d in [1.5, 2.0, 3.0, 4.0]
            drf_test = DRationdlFunc(Float64; d=d)
            a_test = d / 2
            coef = gamma(a_test - 0.5) / (sqrt(π) * gamma(a_test))
            for x in [0.5, 1.0, 2.0]
                ratio = Hfunc(x, drf_test) / origfunc(x, drf_test)
                expected_ratio = x * coef
                @test ratio ≈ expected_ratio rtol = 1e-10
            end
        end

        # For d=1 (a=1/2), check the special formula
        drf1 = DRationdlFunc(Float64; d=1.0)
        for x in [0.5, 1.0, 2.0]
            expected = (2 / π) * asinh(x) / sqrt(1 + x^2)
            @test Hfunc(x, drf1) ≈ expected rtol = 1e-10
        end

        # Verify well-known result: H[1/(1+x²)](x) = x/(1+x²)
        drf2 = DRationdlFunc(Float64; d=2.0)
        for x in [0.5, 1.0, 2.0, 3.0]
            @test Hfunc(x, drf2) ≈ x / (1 + x^2) rtol = 1e-10
        end
    end

    @testset "Limiting behavior" begin
        # As d → 0, origfunc → 1 and Hfunc → 0 (approximately)
        drf_small = DRationdlFunc(Float64; d=0.01)
        @test origfunc(1.0, drf_small) ≈ 1.0 rtol = 0.01
        @test abs(Hfunc(1.0, drf_small)) < 0.1  # small value

        # As x → ∞, the ratio |Hfunc|/|origfunc| → |x * coefficient|
        x_large = 1000.0

        # Case a < 1/2
        drf = DRationdlFunc(Float64; d=0.5)
        a = 0.25
        ratio = abs(Hfunc(x_large, drf) / origfunc(x_large, drf))
        @test ratio ≈ abs(x_large * tanpi(a)) rtol = 0.001

        # Case a > 1/2
        drf2 = DRationdlFunc(Float64; d=2.0)
        a2 = 1.0
        coef = gamma(a2 - 0.5) / (sqrt(π) * gamma(a2))
        ratio2 = abs(Hfunc(x_large, drf2) / origfunc(x_large, drf2))
        @test ratio2 ≈ abs(x_large * coef) rtol = 0.001
    end

    @testset "Known special values" begin
        # H[1/(1+x²)](x) = x/(1+x²) for a=1 (d=2)
        drf_a1 = DRationdlFunc(Float64; d=2.0)
        for x in [-3.0, -1.0, 0.5, 1.0, 2.0, 5.0]
            @test Hfunc(x, drf_a1) ≈ x / (1 + x^2)
        end

        # H[1/(1+x²)²](x) = x(3+x²)/(2(1+x²)²) for a=2 (d=4)
        # Coefficient: Γ(3/2)/(√π*Γ(2)) = (√π/2)/(√π*1) = 1/2
        # So H = x * (1+x²)^(-2) * 1/2 ... but the actual formula should be verified
        drf_a2 = DRationdlFunc(Float64; d=4.0)
        a = 2.0
        coef = gamma(a - 0.5) / (sqrt(π) * gamma(a))  # = Γ(1.5)/(√π*Γ(2)) = (√π/2)/(√π*1) = 0.5
        @test coef ≈ 0.5
        for x in [0.5, 1.0, 2.0]
            expected = x * (1 + x^2)^(-2) * coef
            @test Hfunc(x, drf_a2) ≈ expected
        end
    end
end
