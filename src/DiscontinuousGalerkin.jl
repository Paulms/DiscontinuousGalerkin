__precompile__()
module DiscontinuousGalerkin
  using Polynomials
  using FastGaussQuadrature
  using LsqFit

  include("mesh.jl")
  include("problem.jl")
  include("riemann_solvers.jl")
  include("basis.jl")

  export Uniform1DMesh
  export legendre_basis, project_function
  export get_local_mass_matrix, get_local_inv_mass_matrix
  export DGScalar1DProblem
end
