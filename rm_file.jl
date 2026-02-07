using Pkg
Pkg.activate(".")

dir0_vec = ["./aaa/", "./asy/", "./fir/", "./fft", "./loglog", "./logasy"]

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
