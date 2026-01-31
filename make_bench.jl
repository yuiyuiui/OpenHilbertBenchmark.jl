using Pkg
Pkg.activate(".")

dir0_vec = ["./aaa+fir/", "./asy+fir/", "./fir/", "./fft"]

# Collect all .txt and .svg files to delete
files_to_delete = String[]
test_files = String[]

for dir0 in dir0_vec
    if !isdir(dir0)
        @warn "Directory $dir0 does not exist, skipping..."
        continue
    end

    # Walk through all subdirectories
    for (root, dirs, files) in walkdir(dir0)
        for file in files
            filepath = joinpath(root, file)
            if endswith(file, ".txt") || endswith(file, ".svg")
                push!(files_to_delete, filepath)
            elseif file == "test.jl"
                push!(test_files, filepath)
            end
        end
    end
end

# Delete all .txt and .svg files
println("Deleting $(length(files_to_delete)) files (.txt and .svg)...")
for file in files_to_delete
    try
        rm(file)
        println("  Deleted: $file")
    catch e
        @warn "Failed to delete $file: $e"
    end
end

# Run all test.jl files
println("\nRunning $(length(test_files)) test.jl files...")
# Get the project root directory (where make_bench.jl is located)
project_root = dirname(@__FILE__)
for (idx, test_file) in enumerate(test_files)
    println("\n[$(idx)/$(length(test_files))] Running: $test_file")
    try
        # Convert to absolute path and ensure directory exists
        test_file_abs = abspath(project_root, test_file)
        test_dir = dirname(test_file_abs)
        mkpath(test_dir)  # Ensure directory exists for write_setting

        # Run from project root to ensure relative paths work correctly
        cd(project_root) do
            return include(test_file_abs)
        end
        println("  ✓ Completed: $test_file")
    catch e
        @error "Failed to run $test_file: $e"
        rethrow(e)
    finally
        # Force garbage collection to free memory between tests
        GC.gc()
    end
end

println("\n✓ All done!")
