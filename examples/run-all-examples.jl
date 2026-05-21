using WhereTheWaterFlows

const EXAMPLES_DIR = @__DIR__

function collect_examples(root::AbstractString)
    files = String[]
    for (dir, _, fs) in walkdir(root)
        for f in fs
            endswith(f, ".jl") || continue
            fl = joinpath(dir, f)
            basename(fl) == "run-all-examples.jl" && continue
            push!(files, fl)
        end
    end
    sort!(files)
    return files
end

function run_example(path::AbstractString)
    modname = Symbol("AllExamples_", hash(path))
    m = Module(modname)
    Core.eval(m, :(using WhereTheWaterFlows))
    Base.include(m, path)
    return nothing
end

files = collect_examples(EXAMPLES_DIR)

println("Running all examples from: ", EXAMPLES_DIR)
println("Found ", length(files), " scripts")

for fl in files
    println("\n==> ", relpath(fl, EXAMPLES_DIR))
    run_example(fl)
end

println("\nDone.")
