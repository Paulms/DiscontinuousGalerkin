# Simple Advection problem

# build mesh
include("../src/DiscontinuousGalerkin.jl")
using DiscontinuousGalerkin
N = 2
mesh = Uniform1DMesh(N,-1.0,1.0)
#set initial condition
f0(x) = sin(x-π/3)^3+1
mybasis=legendre_basis(3)
u0 = zeros(mybasis.order+1,N)
for i = 1:N
  value = project_function(f0,mybasis,(mesh.cell_faces[i],mesh.cell_faces[i+1]))
  u0[:,i] = value.param
end
