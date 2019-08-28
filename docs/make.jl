using Documenter, SeisNetwork

makedocs(;
    modules=[SeisNetwork],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/kura-okubo/SeisNetwork.jl/blob/{commit}{path}#L{line}",
    sitename="SeisNetwork.jl",
    authors="kurama",
    assets=String[],
)

deploydocs(;
    repo="github.com/kura-okubo/SeisNetwork.jl",
)
