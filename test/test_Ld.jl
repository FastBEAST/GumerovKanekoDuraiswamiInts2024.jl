#test transposed_doublelayer
using CompScienceMeshes
using GumerovKanekoDuraiswamiInts2024
using LinearAlgebra
using Test
#an example of mesh: Mesh{3, 3, Float64}(SVector{3, Float64}[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], SVector{3, Int64}[[1, 2, 3]], Dict{SVector{3, Int64}, Int64}())

#make a sphere
Γ1 = meshsphere(; radius=1.0, h=1.0)
Γ2 = meshsphere(; radius=1.0, h=1.0)
#here for testing

#get result from GL
@time begin
    Ag = IntegrateMesh(Γ1, Γ2; operator="doublelayer_transposed")
end
Ab = IntegrateMesh(Γ1, Γ2; operator="doublelayer")'
#checking how close the values are
diff1 = Ag .- Ab
diff_norm1 = norm(diff1)
perc_diff1 = norm(diff1) / norm(Ab)

@test perc_diff1 < 1e-6
@test diff_norm1 < 1e-6
