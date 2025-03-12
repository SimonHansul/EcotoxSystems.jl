using Pkg; Pkg.activate("docs")
using Documenter
using EcotoxSystems

makedocs(
    sitename = "EcotoxSystems.jl", 
    format = Documenter.HTML()
    )

deploydocs(
    repo = "github.com/USER_NAME/PACKAGE_NAME.jl.git",
)