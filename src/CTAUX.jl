module CTAUX
export ctaux_main
using Random
using FFTW
using Dierckx
using LinearAlgebra
using InteractiveUtils

include("mesh.jl")
include("parameterstructs.jl")
include("greenfunctions.jl")

include("qmc.jl")
include("remove.jl")
include("insert.jl")
# Write your package code here.

end
