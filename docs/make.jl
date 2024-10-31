using Pkg; Pkg.activate("docs")
using Documenter
using EcotoxSystems

makedocs(
    sitename = "EcotoxSystems.jl", 
    format = Documenter.HTML()
    )