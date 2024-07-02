using BEAST
using CompScienceMeshes
using LinearAlgebra
using PlotlyJS
using GumerovKanekoDuraiswamiInts2024
using IterativeSolvers
o, x, y, z = euclidianbasis(3)

Γ = meshsphere(radius=1.0, h=1.0)
# plot(wireframe(Γ))

# X = raviartthomas(Γ)
X = lagrangecxd0(Γ)

ui = Helmholtz3D.linearpotential(direction = z)
e = strace(ui,Γ)

κ = 0.0
# t = Helmholtz3D.singlelayer(wavenumber=κ)
t = Helmholtz3D.singlelayer(;)

# A = assemble(t, X, X)
A = IntegrateMesh(Γ, Γ)
b = assemble(e, X)
u = IterativeSolvers.gmres(A,b)
fcr, geo = facecurrents(u,X)

pt =PlotlyJS.plot(patch(geo,real.(fcr)))
display(pt)