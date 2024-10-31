using Pkg; Pkg.activate("test")
using Aqua
using Revise, EcotoxSystems

Aqua.test_stale_deps(EcotoxSystems)

# where is the source of the persistant task?
#   - utils.jl : yes
#       - reformat_labels : no
#       - getclonames: no
#       - extract_colnames: no
#   
#   - drcfuncts.jl : 

# seems to be in robustmin or robustmean


Aqua.test_persistent_tasks(EcotoxSystems)


# Aqua.find_persistent_tasks_deps(EcotoxSystems) # this returns no positive results

Aqua.test_all(EcotoxSystems)

