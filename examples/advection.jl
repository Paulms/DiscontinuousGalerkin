# Simple Advection problem

# build mesh
include("../src/DiscontinuousGalerkin.jl")
using DiscontinuousGalerkin

const CFL = 0.5
const Tend = 1.0

f(::Type{Val{:jac}},u::Vector) = eye(size(u,1))
f(u::Vector) = u

#define max wave speed
max_w_speed(u, f) = 1

f0(x) = sin(2*π*x)

function get_problem(N)
  mesh = Uniform1DFVMesh(N,-2.0,2.0,:PERIODIC,:PERIODIC)
  ConservationLawsProblem(f0,f,CFL,Tend,mesh)
end
#Compile
prob = get_problem(10)

#Run
prob = get_problem(200)

#Generate basis for Discontinuous Galerkin Scheme
basis=legendre_basis(3)

#Optional dgLimiter
limiter! = DGLimiter(prob, basis, Linear_MUSCL_Limiter())

#Solve problem
@time sol = solve(prob, DiscontinuousGalerkinScheme(basis, advection_num_flux; max_w_speed = max_w_speed); TimeIntegrator = SSPRK22(limiter!))

#Plot
using Plots
plot(sol)
