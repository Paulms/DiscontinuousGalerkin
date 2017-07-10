__precompile__()
module DiscontinuousGalerkin
  using Polynomials
  using FastGaussQuadrature
  using LsqFit
  using Parameters, Reexport
  using DiffEqBase
  @reexport using OrdinaryDiffEq
  using RecipesBase

  import DiffEqBase: solve, @def

  include("mesh.jl")
  include("problem.jl")
  include("riemann_solvers.jl")
  include("basis.jl")
  include("solve.jl")
  include("solution.jl")
  include("plots.jl")

  export Uniform1DMesh
  export legendre_basis
  export DG1DProblem
  export advection_solver, rusanov_solver
  export solve, DiscontinuousGalerkinScheme
end
