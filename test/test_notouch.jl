#testing of the no touch cases, compare table 1 in Gumerov et al., 2023
using GumerovKanekoDuraiswamiInts2024
using LinearAlgebra
using Test

x1 = [0.0, 0.0, 0.0]
x2 = [1.0, 0.0, 0.0]
x3 = [1 / 2, sqrt(3) / 2, 0.0]
y1 = [1.0, 0.0, 1.0]
y2 = [0.0, 0.0, 1.0]

#test case 1
y3 = [1 / 2, 0.0, 1 + sqrt(3) / 2]

L, M, Ld, Md = GalerkinLaplaceTriGS(x1, x2, x3, y1, y2, y3; L_output=1, Md_output=1);

@test L ≈ 0.139757030669707
@test M ≈ 0.099860729206614
@test Ld ≈ 0.022035244796804
@test Md ≈ 0.046564310284965

#test case 2
y3 = [1 / 2, sqrt(6) / 4, 1 + sqrt(6) / 4]

L, M, Ld, Md = GalerkinLaplaceTriGS(x1, x2, x3, y1, y2, y3; L_output=1, Md_output=1);

ny = cross(y3 - y1, y3 - y2) / norm(cross(y3 - y1, y3 - y2))
@test L ≈ 0.149630247150535
@test M ≈ 0.114715727210190
@test Ld ≈ ny ⋅ [0.0, 0.010953212167802, -0.114715727210190]
@test Md ≈ 0.137859073743097

#test case 3
y3 = [1 / 2, -sqrt(3) / 2, 1]

L, M, Ld, Md = GalerkinLaplaceTriGS(x1, x2, x3, y1, y2, y3; L_output=1, Md_output=1);
@test L ≈ 0.156068357679434
@test M ≈ 0.111863573921226
@test Ld ≈ -0.111863573921226
@test Md ≈ -0.138417139905960
