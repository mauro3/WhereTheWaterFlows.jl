using Documenter, WhereTheWaterFlows

makedocs(;
    modules=[WhereTheWaterFlows],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/mauro3/WhereTheWaterFlows.jl/blob/{commit}{path}#L{line}",
    sitename="WhereTheWaterFlows.jl",
    authors="Mauro Werder <mauro3@runbox.com>",
    assets=String[],
)

deploydocs(;
    repo="github.com/mauro3/WhereTheWaterFlows.jl",
)
