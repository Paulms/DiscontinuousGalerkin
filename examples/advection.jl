# Simple Advection problem

# build mesh
include("../src/DiscontinuousGalerkin.jl")
using DiscontinuousGalerkin
N = 10
mesh = Uniform1DMesh(N,-2.0,2.0)
#define initial condition
f0(x) = sin(2*π*x)

#define flux function
f(u) = u

#define max wave speed
max_w_speed(u) = 1

#Setup problem
tspan = (0.0,1.0)
cfl = 0.5
problem = DGScalar1DProblem(f0, f, max_w_speed, mesh, tspan, cfl)

#Generate basis for Discontinuous Galerkin Scheme
basis=legendre_basis(3)

#Solve problem
u = solve(problem, DiscontinuousGalerkinScheme(basis, advection_solver))

#Plot
using Plots
gr()
plot(u)
