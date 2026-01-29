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
        rtfpr = RationalFuncPolesRepresent(norder=1,
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
    end

    @testset "get_funcname for DRationdlFunc" begin
        drf = DRationdlFunc(Float64; d=0.5)

        name = get_funcname(drf)
        @test occursin("Rationdl", name)
        @test occursin("0.5", name)

        drf2 = DRationdlFunc(Float64; d=2.0)
        name2 = get_funcname(drf2)
        @test occursin("2.0", name2)
    end

    @testset "get_funcname for MixedFunc" begin
        swf = SchwartzFunc(Float64; A=1.0, μ=0.0, σ=1.0, d=2)
        drf = DRationdlFunc(Float64; d=0.5)

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
        @test get_algname(nothing, nothing, nothing) == ""

        # Test returns a string (actual content depends on OpenHilbert types)
        @test get_algname(nothing, nothing, nothing) isa String
    end

    @testset "write_setting" begin
        # Create a temporary directory for testing
        tmpdir = mktempdir()

        try
            swf = SchwartzFunc(Float64; A=1.0, μ=0.0, σ=1.0, d=2)
            L0_vec = [4.0, 8.0, 16.0]
            grid_gap = 0.1
            points_density = 32

            # Call write_setting with nothing for dm, pola, trans
            write_setting(swf, L0_vec, grid_gap, points_density;
                          file_place=tmpdir, dm=nothing, pola=nothing, trans=nothing)

            # Check that file was created
            setting_file = joinpath(tmpdir, "setting.json")
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
            drf = DRationdlFunc(Float64; d=0.5)
            mf = MixedFunc(Float64; swf=swf, drtf=drf)

            write_setting(mf, L0_vec, grid_gap, points_density;
                          file_place=tmpdir, dm=nothing, pola=nothing, trans=nothing)

            content2 = read(setting_file, String)
            # With details=true, MixedFunc returns detailed expressions
            # Check for components: SchwartzFunc details contain "exp", DRationdlFunc contains "Rationdl"
            @test occursin("exp", content2)  # SchwartzFunc detailed expression
            @test occursin("Rationdl", content2)  # DRationdlFunc name
            @test occursin("func_type:", content2)
        finally
            # Cleanup
            rm(tmpdir; recursive=true)
        end
    end
end
