using Pkg; Pkg.activate("docs")
using Documenter
using EcotoxSystems

makedocs(
    sitename = "EcotoxSystems.jl", 
    format = Documenter.HTML()
    )

deploydocs(
    repo = "github.com/SimonHansul/EcotoxSystems.jl.git",
)