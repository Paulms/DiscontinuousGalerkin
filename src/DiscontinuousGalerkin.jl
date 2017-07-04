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
end
