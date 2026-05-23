using Documenter, WhereTheWaterFlows, CairoMakie

makedocs(;
    modules = [WhereTheWaterFlows],
    format  = Documenter.HTML(),
    pages   = [
        "Home"     => "index.md",
        "Tutorial" => "tutorial.md",
        "Subglacial flow `Subglacially`" => "subglacially.md",
        "Monte Carlo `Randomly`" => "randomly.md",
        "Worked Example" => "worked-example.md",
        "Example scipts" => "examples.md",
        "Feedback function" => "feedback.md",
        "API"      => "api.md",
        "References"      => "references.md",
    ],
    repo     = "https://github.com/mauro3/WhereTheWaterFlows.jl/blob/{commit}{path}#L{line}",
    sitename = "WhereTheWaterFlows.jl",
    authors  = "Mauro Werder <mauro3@runbox.com>",
)

deploydocs(;
    repo         = "github.com/mauro3/WhereTheWaterFlows.jl.git",
    push_preview = true,
)
