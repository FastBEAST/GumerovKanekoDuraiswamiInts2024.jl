#testing of two shared vertices, compare section 7.2 in Gumerov et al., 2023
using GumerovKanekoDuraiswamiInts2024
using LinearAlgebra
using Test

x1 = [0.0, 0.0, 0.0]
x2 = [1.0, 0.0, 0.0]
x3 = [1 / 2, sqrt(3) / 2, 0.0]
y1 = [0.0, 0.0, 0.0]
y2 = [1.0, 0.0, 0.0]
y3 = [1 / 2, 0.0, sqrt(3) / 2]

L, M, Ld, Md = GalerkinLaplaceTriGS(x1, x2, x3, y1, y2, y3; L_output=1, Md_output=1);

@test L ≈ 0.415922738854561
@test M ≈ 0.706739910625218
@test Md ≈ 2.85747144125268
