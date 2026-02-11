@testset "Benchreport" begin
    @testset "get_funcname for SchwartzFunc" begin
        swf = SchwartzFunc(Float64; A=1.0, μ=0.0, σ=1.0, d=2)

        # Without details
        @test get_funcname(swf) == "Schwartz"
        @test get_funcname(swf; details=false) == "Schwartz"

        # With details
        detailed = get_funcname(swf; details=true)
        @test occursin("exp", detailed)
    end

    @testset "get_funcname for RationalFuncPolesRepresent" begin
        Random.seed!(42)
        T = Float64
        rtfpr = RationalFuncPolesRepresent(T; norder=1,
                                           poles_scale=1.0,
                                           npole=[2],
                                           amplitudes=[[0.5, 0.5]],
                                           imag_gap=0.5)

        # Without details
        @test get_funcname(rtfpr) == "Rational"
        @test get_funcname(rtfpr; details=false) == "Rational"

        # With details
        detailed = get_funcname(rtfpr; details=true)
        @test occursin("real(", detailed)

        @test rtfpr isa RationalFuncPolesRepresent{T}
    end

    @testset "get_funcname for DRationdlFunc" begin
        drf = DRationdlFunc([0.5], Float64[])

        name = get_funcname(drf)
        @test occursin("Rationdl", name)
        @test occursin("even", name)

        drf2 = DRationdlFunc([2.0], Float64[])
        name2 = get_funcname(drf2)
        @test occursin("Rationdl", name2)
    end

    @testset "get_funcname for MixedFunc" begin
        swf = SchwartzFunc(Float64; A=1.0, μ=0.0, σ=1.0, d=2)
        drf = DRationdlFunc([0.5], Float64[])

        mf = MixedFunc(Float64; swf=swf, drtf=drf)

        name = get_funcname(mf)
        @test occursin("Schwartz", name)
        @test occursin("Rationdl", name)

        # Single component
        mf_single = MixedFunc(Float64; swf=swf)
        name_single = get_funcname(mf_single)
        @test occursin("Schwartz", name_single)
    end

    @testset "get_algname" begin
        # Test with nothing values
        @test get_algname(nothing, nothing) == ""

        # Test returns a string (actual content depends on OpenHilbert types)
        @test get_algname(nothing, nothing) isa String
    end

    @testset "write_setting" begin
        # Create a temporary directory for testing
        tmpdir = mktempdir()

        try
            swf = SchwartzFunc(Float64; A=1.0, μ=0.0, σ=1.0, d=2)
            L0_vec = [4.0, 8.0, 16.0]
            grid_gap = 0.1
            points_density = 32

            # Call write_setting with nothing for tdm, trans
            write_setting(swf, L0_vec, grid_gap, points_density;
                          file_place=tmpdir, tdm=nothing, trans=nothing)

            # Check that file was created
            setting_file = joinpath(tmpdir, "setting.txt")
            @test isfile(setting_file)

            # Read and verify content
            content = read(setting_file, String)
            @test occursin("Used L0:", content)
            @test occursin("4.0", content)
            @test occursin("8.0", content)
            @test occursin("16.0", content)
            @test occursin("grid_gap:", content)
            @test occursin("0.1", content)
            @test occursin("Points density:", content)
            @test occursin("32", content)
            @test occursin("func_type:", content)
            @test occursin("exp", content)  # SchwartzFunc details contain "exp"
            @test occursin("alg:", content)

            # Test with MixedFunc
            drf = DRationdlFunc([0.5], Float64[])
            mf = MixedFunc(Float64; swf=swf, drtf=drf)

            write_setting(mf, L0_vec, grid_gap, points_density;
                          file_place=tmpdir, tdm=nothing, trans=nothing)

            content2 = read(setting_file, String)
            # With details=true, MixedFunc returns detailed expressions
            @test occursin("exp", content2)  # SchwartzFunc detailed expression
            @test occursin("func_type:", content2)
        finally
            # Cleanup
            rm(tmpdir; recursive=true)
        end
    end

    @testset "cal_Hlogrtf_nume (LogRationalFunc numerical approximation solution) basic properties" begin
        T = Float64
        L0 = T(2.0)
        point_density = 8
        d = 1

        H_big = OpenHilbertBenchmark.cal_Hlogrtf_nume(L0, point_density, d)
        N_big = round(Int, point_density * (100 * L0)) * 2 + 1

        @test eltype(H_big) == T
        @test length(H_big) == N_big
        @test all(isfinite, H_big)

        # when f(x) is an even function, the Hilbert transform should be an odd function: H(-x) = -H(x)
        mid_big = length(H_big) ÷ 2 + 1
        @test abs(H_big[mid_big]) ≤ 1e-6
        @test maximum(abs.(H_big .+ reverse(H_big))) ≤ 1e-5

        # take the center slice to match the target grid length N0 (corresponding to [-L0, L0])
        N0 = round(Int, point_density * L0) * 2 + 1
        H0 = H_big[(mid_big - N0 ÷ 2):(mid_big + N0 ÷ 2)]
        @test length(H0) == N0
        @test abs(H0[N0 ÷ 2 + 1]) ≤ 1e-6
        @test maximum(abs.(H0 .+ reverse(H0))) ≤ 1e-5

        # d change should result in numerical result change (basic sanity check)
        H_big_d2 = OpenHilbertBenchmark.cal_Hlogrtf_nume(L0, point_density, 2)
        @test length(H_big_d2) == length(H_big)
        @test norm(H_big_d2 - H_big) > 0
    end

    @testset "MixedFunc + LogRationalFunc numerical approximation solution concatenation (non-log part + log numerical slice)" begin
        T = Float64
        L0 = T(2.0)
        point_density = 8
        h = T(1 / point_density)
        N0 = round(Int, point_density * L0) * 2 + 1
        x0 = T.((-N0 ÷ 2):(N0 ÷ 2)) .* h

        # non-log part: select a component with closed-form (Schwartz d=2)
        swf = SchwartzFunc(T; A=1.0, μ=0.0, σ=1.0, d=2)
        mf_nonlog = MixedFunc(T; swf=swf, logrtf=nothing)
        H_nonlog = [Hfunc(xi, mf_nonlog) for xi in x0]
        @test length(H_nonlog) == N0
        @test all(isfinite, H_nonlog)

        # log part: use the center slice of cal_Hlogrtf_nume to align with x0
        H_log_big = OpenHilbertBenchmark.cal_Hlogrtf_nume(L0, point_density, 1)
        mid_log = length(H_log_big) ÷ 2 + 1
        H_log0 = H_log_big[(mid_log - N0 ÷ 2):(mid_log + N0 ÷ 2)]
        @test length(H_log0) == N0

        # the concatenated solution should have the same length and be finite
        H_ref = H_nonlog .+ H_log0
        @test length(H_ref) == N0
        @test all(isfinite, H_ref)
    end
end
