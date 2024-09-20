# GumerovKanekoDuraiswamiInts2024

[![Build Status](https://github.com/FastBEAST/GumerovKanekoDuraiswamiInts2024.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/FastBEAST/GumerovKanekoDuraiswamiInts2024.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/FastBEAST/GumerovKanekoDuraiswamiInts2024.jl/graph/badge.svg?token=RDRQTBWQS3)](https://codecov.io/gh/FastBEAST/GumerovKanekoDuraiswamiInts2024.jl)

This repository contains a Julia implementation of the analytical solutions of the double surface integrals for all layer potential kernels associated with the Laplace Green's function in the Galerkin boundary element method in $\boldsymbol{\mathbb{R}^3}$ for constant basis functions on flat triangles as presented in [1].

The implementation also took inspiration from the already existing Matlab code published with the paper which can be found [here](https://github.com/pirl-lab/analytical-quadrature-laplace-galerkin).

## Available Integrals
For the Laplace kernel
$$G(\vec r,\vec r')\frac{1}{4\pi |\vec r - \vec r'|}\,,$$
the four integrals for the single layer, double layer, transposed double layer and hypersingular operator are available:
$$L = \int_{S_y}\int_{S_x}G(\vec r, \vec r')dS(\vec r)dS(\vec r')\,,$$
$$M = \int_{S_y}\int_{S_x}\vec n_{\vec r} \cdot\nabla_{\vec r} G(\vec r, \vec r')dS(\vec r)dS(\vec r')\,,$$
$$L' = \int_{S_y}\vec n_{\vec r'} \cdot\nabla_{\vec r'}\int_{S_x} G(\vec r, \vec r')dS(\vec r)dS(\vec r')\,,$$
$$M' = \int_{S_y}\vec n_{\vec r'} \cdot\nabla_{\vec r'}\int_{S_x} \vec n_{\vec r} \cdot\nabla_{\vec r}G(\vec r, \vec r')dS(\vec r)dS(\vec r')\,.$$
Here, $\vec n_{\vec r}$ and $\vec n_{\vec r'}$ are the normal vectors of the flat triangles, and $\nabla_{\vec r}$ and $\nabla_{\vec r'}$ are the gradients with respect to position $\vec r$ and $\vec r'$.

## Example of Usage
The package can be installed by

```@julia
import Pkg
Pkg.add(url="https://github.com/FastBEAST/GumerovKanekoDuraiswamiInts2024.jl.git")
```
Given two triangles with vertices $x_1$, $x_2$, $x_3$ and $y_1$, $y_2$, $y_3$, the results for all four integrals can be calculated as follows:
```@julia
using GumerovKanekoDuraiswamiInts2024

x1 = [0.0, 0.0, 0.0]
x2 = [1.0, 0.0, 0.0]
x3 = [1 / 2, sqrt(3) / 2, 0.0]
y1 = [1.0, 0.0, 1.0]
y2 = [0.0, 0.0, 1.0]
y3 = [1 / 2, 0.0, 1 + sqrt(3) / 2]
L, M, Ld, Md = GalerkinLaplaceTriGS(x1, x2, x3, y1, y2, y3; L_output=1, Md_output=1);
```
By setting Md_output to zero, the calculation of the hypersingular part is supressed, by setting L_output to zero, the calculation of the other three operators does not take place. These three are always calculated together because they require the same fundamental calculations.

If an entire boundary element method type matrix for one of the operators has to be assembled, the function "IntegrateMesh" can be used:
```@julia
using GumerovKanekoDuraiswamiInts2024
using CompScienceMeshes
using StaticArrays

Mesh1 = meshsphere(1.0,0.3)
Mesh2 = CompScienceMeshes.translate(Mesh1, SVector(2.0,0.0,0.0))

L = IntegrateMesh(Γ1, Γ2; operator="singlelayer")
```
The available operator keywords are "singlelayer", "doublelayer", "doublelayer_transposed" and "hypersingular".

## References
[1] Gumerov, Nail A., Shoken Kaneko, and Ramani Duraiswami. 2024. “Analytical Galerkin Boundary Integrals of Laplace Kernel Layer Potentials in $\boldsymbol{\mathbb{R}^3}$.” SIAM Journal on Scientific Computing 46 (2): A974–97. https://doi.org/10.1137/23M1547688.
