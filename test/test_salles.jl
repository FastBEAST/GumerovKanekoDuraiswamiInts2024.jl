
using GumerovKanekoDuraiswamiInts2024
using LinearAlgebra
using Test

x1 = [0.0, 0.0, 0.0]
x2 = [0.0, 1.0, 0.0]
x3 = [1.0, 0.0, 0.0]
y1 = [0.0, 0.0, 0.01]
y2 = [0.0, 1.0, 0.01]
y3 = [-1.0, 0.0, 0.01]

L, M, Ld, Md = GalerkinLaplaceTriGS(x1, x2, x3, y1, y2, y3; L_output=1, Md_output=1);

@test L ≈ 0.4150963397038614

y1 = [-2.0, 0.5, 0.01]
y2 = [-1.0, 1.0, 0.01]
y3 = [-1.0, 0.0, 0.01]

ny = cross(y3 - y1, y3 - y2) / norm(cross(y3 - y1, y3 - y2))

L, M, Ld, Md = GalerkinLaplaceTriGS(y1, y2, y3, x1, x2, x3; L_output=1, Md_output=1);

@test Ld ≈ -0.0006289369951278
