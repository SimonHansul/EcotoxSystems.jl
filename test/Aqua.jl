using Pkg; Pkg.activate("test")
using Aqua
using Revise, EcotoxSystems

Aqua.test_stale_deps(EcotoxSystems)

Aqua.test_persistent_tasks(EcotoxSystems, tmax = 10)


# Aqua.find_persistent_tasks_deps(EcotoxSystems) # this returns no positive results

Aqua.test_all(EcotoxSystems)

