using Pkg; Pkg.activate("test")
using Aqua
using Revise, EcotoxSystems

Aqua.test_all(EcotoxSystems)

