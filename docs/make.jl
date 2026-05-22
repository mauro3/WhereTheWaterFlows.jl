using Documenter, WhereTheWaterFlows, CairoMakie

makedocs(;
    modules = [WhereTheWaterFlows],
    format  = Documenter.HTML(),
    pages   = [
        "Home"     => "index.md",
        "Tutorial" => "tutorial.md",
        "Feedback" => "feedback.md",
        "Randomly" => "randomly.md",
        "Subglacially" => "subglacially.md",
        "Examples" => "examples.md",
        "API"      => "api.md",
    ],
    repo     = "https://github.com/mauro3/WhereTheWaterFlows.jl/blob/{commit}{path}#L{line}",
    sitename = "WhereTheWaterFlows.jl",
    authors  = "Mauro Werder <mauro3@runbox.com>",
)

deploydocs(;
    repo         = "github.com/mauro3/WhereTheWaterFlows.jl.git",
    push_preview = true,
)
