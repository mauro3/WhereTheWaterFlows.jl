import Downloads
import Random

function print_usage(scriptname::String, defaults_text::String)
    println("Usage: julia --project=benchmarks benchmarks/$(scriptname) [--mode real|mock] [--dataset small|large] [--runs N]")
    println("Defaults: $(defaults_text)")
end

function parse_common_args(
    args::Vector{String};
    scriptname::String,
    defaults_text::String,
    default_mode::String = "real",
    default_dataset::String = "small",
    default_runs::Union{Nothing, Int} = 6,
    mock_default_runs::Union{Nothing, Int} = 1,
)
    mode = default_mode
    dataset = default_dataset
    runs = nothing

    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--help" || arg == "-h"
            print_usage(scriptname, defaults_text)
            exit(0)
        elseif arg == "--mode"
            i += 1
            i <= length(args) || error("Missing value for --mode")
            mode = lowercase(args[i])
        elseif arg == "--dataset"
            i += 1
            i <= length(args) || error("Missing value for --dataset")
            dataset = lowercase(args[i])
        elseif arg == "--runs"
            i += 1
            i <= length(args) || error("Missing value for --runs")
            runs = parse(Int, args[i])
        else
            error("Unknown argument: $arg")
        end
        i += 1
    end

    mode in ("real", "mock") || error("--mode must be real or mock")
    dataset in ("small", "large") || error("--dataset must be small or large")

    if isnothing(runs)
        if mode == "mock" && !isnothing(mock_default_runs)
            runs = mock_default_runs
        else
            runs = default_runs
        end
    end

    isnothing(runs) && error("No default runs configured for this mode")
    runs >= 1 || error("--runs must be >= 1")

    return (; mode, dataset, runs)
end

function create_mock_dem(dataset::String)
    n = dataset == "small" ? 128 : 256
    Random.seed!(42)
    x = range(-pi, pi, length=n)
    y = range(-pi, pi, length=n)
    return Float32.(sin.(x) .* cos.(y') .+ 0.05 .* randn(length(x), length(y)))
end

function ensure_small_dem(datadir::String)
    dem_small = joinpath(datadir, "swissalti3d_tile.tif")

    if !isfile(dem_small)
        stac_url = "https://data.geo.admin.ch/api/stac/v1/collections/ch.swisstopo.swissalti3d/items?bbox=7.7,46.5,7.8,46.6&limit=20"
        txt = Downloads.download(stac_url) |> read |> String

        urls = String.(m.match for m in eachmatch(r"https://[^\" ]+\.tif", txt))
        candidates = filter(u -> occursin("_0.5_2056_", u), urls)
        isempty(candidates) && error("No swissALTI3D 0.5 m GeoTIFF found.")

        url = first(candidates)
        println("Downloading DEM:")
        println(url)
        Downloads.download(url, dem_small)
    end

    return dem_small
end
