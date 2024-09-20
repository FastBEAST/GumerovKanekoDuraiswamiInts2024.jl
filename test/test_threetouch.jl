#testing of self interactions, compare section 7.2 in Gumerov et al., 2023
using GumerovKanekoDuraiswamiInts2024
using LinearAlgebra
using Test

x1 = y1 = [0.0, 0.0, 0.0]
x2 = y2 = [1.0, 0.0, 0.0]
x3 = y3 = [1 / 2, sqrt(3) / 2, 0.0]

L, M, Ld, Md = GalerkinLaplaceTriGS(x1, x2, x3, y1, y2, y3; L_output=1, Md_output=1);

@test L ≈ 3 / 4 * log(3)
@test M ≈ 0
@test Ld ≈ 0
@test Md ≈ 6 * log(3)
