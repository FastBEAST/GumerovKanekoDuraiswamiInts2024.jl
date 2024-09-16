module GumerovKanekoDuraiswamiInts2024

#importing used modules
using LinearAlgebra
using CompScienceMeshes
using StaticArrays

#defeine globla constants and variables
const zerotol = 3e-16 #2e-16 gives error in determining small h4 = 2.482534153247273e-16 as non zero when it seems to be zero
const testMode::Int = 0 #change this into 1 will enable many print parts within the code.
global case = 0

include("GalerkinLaplaceTriGS.jl")
include("GSorthogonalization.jl")
include("cases.jl")
include("I0.jl")
include("I1.jl")
include("I2s.jl")
include("I2t.jl")
include("I3.jl")
include("I0m.jl")
include("I1m.jl")
include("I2sm.jl")
include("I2tm.jl")
include("I3m.jl")
include("IntegrateMesh.jl")

# export GSorthogonalization_expan
export GalerkinLaplaceTriGS
export IntegrateMesh

end