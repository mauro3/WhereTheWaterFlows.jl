using CSV
using DataFrames
using Dates
using Printf

benchdir = @__DIR__
outroot = joinpath(benchdir, "outputs")
outfile = joinpath(benchdir, "results.md")

function read_csv_if_exists(path::String)
    return isfile(path) ? CSV.read(path, DataFrame) : nothing
end

function detect_cpu_model()
    if Sys.islinux()
        cpuinfo = "/proc/cpuinfo"
        if isfile(cpuinfo)
            for ln in eachline(cpuinfo)
                if startswith(ln, "model name")
                    parts = split(ln, ":"; limit = 2)
                    if length(parts) == 2
                        return strip(parts[2])
                    end
                end
            end
        end
    elseif Sys.isapple()
        try
            return strip(read(`sysctl -n machdep.cpu.brand_string`, String))
        catch
        end
    end
    return "unknown"
end

function detect_cpu_info()
    return (
        platform = string(Sys.KERNEL),
        cpu_model = detect_cpu_model(),
        logical_threads = Sys.CPU_THREADS,
    )
end

function fmt_float(x)
    if ismissing(x)
        return ""
    end
    if x isa Real && isnan(Float64(x))
        return ""
    end
    return @sprintf("%.3f", Float64(x))
end

function dem_size_label(npixels)
    if ismissing(npixels)
        return ""
    end
    n = Int(npixels)
    if n == 2000 * 2000
        return "2000x2000"
    elseif n == 8000 * 8000
        return "8000x8000"
    end
    return ""
end

function infer_dataset(row)
    if hasproperty(row, :dataset)
        v = getproperty(row, :dataset)
        if !ismissing(v)
            return string(v)
        end
    end

    if hasproperty(row, :features)
        m = match(r"\((?:real|mock),\s*(small|large)\)", String(getproperty(row, :features)))
        if m !== nothing
            return m.captures[1]
        end
    end

    if hasproperty(row, :npixels)
        n = Int(getproperty(row, :npixels))
        if n == 2000 * 2000
            return "small"
        elseif n == 8000 * 8000
            return "large"
        end
    end

    return ""
end

function dem_sort_key(row)
    if hasproperty(row, :npixels)
        v = getproperty(row, :npixels)
        if !ismissing(v)
            return Int(v)
        end
    end

    ds = infer_dataset(row)
    if ds == "small"
        return 2000 * 2000
    elseif ds == "large"
        return 8000 * 8000
    end

    return typemax(Int)
end

function nruns_sort_key(row)
    if hasproperty(row, :nruns)
        v = getproperty(row, :nruns)
        if !ismissing(v)
            return Int(v)
        end
    end
    return typemax(Int)
end

function nthreads_sort_key(row)
    if hasproperty(row, :nthreads)
        v = getproperty(row, :nthreads)
        if !ismissing(v)
            return Int(v)
        end
    end
    return typemax(Int)
end

function sort_for_report(df::DataFrame)
    idx = collect(1:nrow(df))
    perm = sortperm(idx; by = i -> (
        dem_sort_key(df[i, :]),
        nruns_sort_key(df[i, :]),
        nthreads_sort_key(df[i, :]),
    ))
    return df[perm, :]
end

function filter_real_runs(df::DataFrame)
    if "features" in names(df)
        return filter(row -> occursin(r"\breal\b", String(row.features)), df)
    end
    return df
end

function add_table!(lines::Vector{String}, df::DataFrame, specs)
    headers = [h for (h, _) in specs]
    cols = [c for (_, c) in specs]

    push!(lines, "| " * join(headers, " | ") * " |")
    push!(lines, "| " * join(fill("---", length(headers)), " | ") * " |")

    for row in eachrow(df)
        vals = String[]
        for c in cols
            v = hasproperty(row, c) ? getproperty(row, c) : missing
            if c in (:runtime_mean_s, :runtime_std_s)
                push!(vals, fmt_float(v))
            elseif c == :dem_size
                npix = hasproperty(row, :npixels) ? getproperty(row, :npixels) : missing
                push!(vals, dem_size_label(npix))
            elseif c == :dataset
                push!(vals, infer_dataset(row))
            else
                push!(vals, string(v))
            end
        end
        push!(lines, "| " * join(vals, " | ") * " |")
    end
end

lines = String[]
push!(lines, "# Benchmark Results")
push!(lines, "")
push!(lines, "Generated: $(Dates.format(now(), DateFormat("yyyy-mm-dd HH:MM:SS")))")
cpu = detect_cpu_info()
push!(lines, "Platform: $(cpu.platform)")
push!(lines, "CPU: $(cpu.cpu_model)")
push!(lines, "Logical threads: $(cpu.logical_threads)")
push!(lines, "")

wwf_csv = joinpath(outroot, "wwf", "benchmark_wwf.csv")
whitebox_csv = joinpath(outroot, "whitebox", "benchmark_whitebox.csv")
topotoolbox_csv = joinpath(outroot, "topotoolbox", "benchmark_topotoolbox.csv")
grass_csv = joinpath(outroot, "grass", "benchmark_grass.csv")

wwf = read_csv_if_exists(wwf_csv)
if wwf === nothing || nrow(wwf) == 0
    push!(lines, "## WWF")
    push!(lines, "")
    push!(lines, "No results found in `outputs/wwf/benchmark_wwf.csv`.")
    push!(lines, "")
else
    push!(lines, "## WWF")
    push!(lines, "")
    wwf = filter_real_runs(wwf)
    if nrow(wwf) == 0
        push!(lines, "No real-run results found in `outputs/wwf/benchmark_wwf.csv`.")
    else
        wwf = sort_for_report(wwf)
        add_table!(
            lines,
            wwf,
            [
                ("method", :method),
                ("dataset", :dataset),
                ("dem_size", :dem_size),
                ("features", :features),
                ("nthreads", :nthreads),
                ("nruns", :nruns),
                ("runtime_mean_s", :runtime_mean_s),
                ("runtime_std_s", :runtime_std_s),
            ],
        )
    end
    push!(lines, "")
end

grass = read_csv_if_exists(grass_csv)
if grass === nothing || nrow(grass) == 0
    push!(lines, "## GRASS")
    push!(lines, "")
    push!(lines, "No results found in `outputs/grass/benchmark_grass.csv`.")
    push!(lines, "")
else
    push!(lines, "## GRASS")
    push!(lines, "")
    grass = filter_real_runs(grass)
    if nrow(grass) == 0
        push!(lines, "No real-run results found in `outputs/grass/benchmark_grass.csv`.")
    else
        grass = sort_for_report(grass)
        add_table!(
            lines,
            grass,
            [
                ("method", :method),
                ("dataset", :dataset),
                ("dem_size", :dem_size),
                ("features", :features),
                ("nthreads", :nthreads),
                ("nruns", :nruns),
                ("runtime_mean_s", :runtime_mean_s),
                ("runtime_std_s", :runtime_std_s),
            ],
        )
    end
    push!(lines, "")
end

whitebox = read_csv_if_exists(whitebox_csv)
if whitebox === nothing || nrow(whitebox) == 0
    push!(lines, "## Whitebox")
    push!(lines, "")
    push!(lines, "No results found in `outputs/whitebox/benchmark_whitebox.csv`.")
    push!(lines, "")
else
    push!(lines, "## Whitebox")
    push!(lines, "")
    whitebox = filter_real_runs(whitebox)
    if nrow(whitebox) == 0
        push!(lines, "No real-run results found in `outputs/whitebox/benchmark_whitebox.csv`.")
    else
        whitebox = sort_for_report(whitebox)
        add_table!(
            lines,
            whitebox,
            [
                ("method", :method),
                ("dataset", :dataset),
                ("dem_size", :dem_size),
                ("features", :features),
                ("nthreads", :nthreads),
                ("nruns", :nruns),
                ("runtime_mean_s", :runtime_mean_s),
                ("runtime_std_s", :runtime_std_s),
            ],
        )
    end
    push!(lines, "")
end

topotoolbox = read_csv_if_exists(topotoolbox_csv)
if topotoolbox === nothing || nrow(topotoolbox) == 0
    push!(lines, "## TopoToolbox")
    push!(lines, "")
    push!(lines, "No results found in `outputs/topotoolbox/benchmark_topotoolbox.csv`.")
    push!(lines, "")
else
    topotoolbox = sort_for_report(topotoolbox)
    push!(lines, "## TopoToolbox")
    push!(lines, "")
    add_table!(
        lines,
        topotoolbox,
        [
            ("method", :method),
            ("dataset", :dataset),
            ("dem_size", :dem_size),
            ("features", :features),
            ("nthreads", :nthreads),
            ("nruns", :nruns),
            ("runtime_mean_s", :runtime_mean_s),
            ("runtime_std_s", :runtime_std_s),
        ],
    )
    push!(lines, "")
end

open(outfile, "w") do io
    write(io, join(lines, "\n"))
end

println("Wrote benchmark results: ", outfile)
