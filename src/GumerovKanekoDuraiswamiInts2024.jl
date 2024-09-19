module GumerovKanekoDuraiswamiInts2024

#importing used modules
using LinearAlgebra
using CompScienceMeshes
using StaticArrays

#define global constants and variables
const zerotol = 3e-16

include("GalerkinLaplaceTriGS.jl")
include("GSorthogonalization.jl")
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
