# Simple Advection problem

# build mesh
include("../src/DiscontinuousGalerkin.jl")
using DiscontinuousGalerkin
N = 10
mesh = Uniform1DMesh(N,-2.0,2.0)
#define initial condition
f0(x) = sin(2*π*x)

#Setup problem
problem = DGScalar1DProblem(f0)

#Assign Initial values
mybasis=legendre_basis(3)
u0 = zeros(mybasis.order+1,N)
for i = 1:N
  value = project_function(problem.f0,mybasis,(mesh.cell_faces[i],mesh.cell_faces[i+1]))
  u0[:,i] = value.param
end

#Add ghost cells
u = hcat(zeros(u0[:,1]),u0,zeros(u0[:,1]))

#build inverse of mass matrix
M_inv = get_local_inv_mass_matrix(mybasis, mesh)

#Reconstruct first u
uₕ = mybasis.φₕ*u
