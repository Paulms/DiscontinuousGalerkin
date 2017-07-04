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

#Reconstruct u in finite space
uₕ = mybasis.φₕ*u

#Reconstruct u on faces
uₛ = mybasis.ψₕ*u

"Calculates residual for DG method
  residual = M_inv*(Q+F)
         where M_inv is the inverse mass matrix
              Q is the edge fluxes
              F is the interior flux"
function residual{T}(basis::Basis{T})
  F = uₕ
  # Integrate interior fluxes
  F = A_mul_B!(F, basis.φₕ.*basis.weights, F)

  # Evaluate edge fluxes
  q = zeros(T,size(uₛ,2)-1)
  for i = 1:(size(uₛ,2)-1)
    ul = @view uₛ[2,i]; ur = @view uₛ[1,i+1]
    q[i] = advection_solver(ul,ur)
  end

end
