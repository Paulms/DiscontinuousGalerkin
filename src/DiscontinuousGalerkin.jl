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
  import Base: show

  abstract type AbstractConservationLawProblem{islinear,isstochastic,MeshType} <: DEProblem end
  abstract type AbstractFEAlgorithm <: DEAlgorithm end
  abstract type AbstractDGLimiter end
  abstract type AbstractFVMesh end

  include("mesh.jl")
  include("problem.jl")
  include("riemann_solvers.jl")
  include("basis.jl")
  include("solve.jl")
  include("solution.jl")
  include("plots.jl")
  include("limiter.jl")

  export Uniform1DFVMesh
  export legendre_basis
  export cell_faces
  export cell_centers, cell_volume, cell_indices, numcells
  export ConservationLawsProblem
  export advection_num_flux, rusanov_euler_num_flux, glf_num_flux
  export legendre_basis, PolynomialBasis, poly_jacobi, poly_legendre, reference_to_interval
  export DGLimiter, Linear_MUSCL_Limiter, WENO_Limiter
  export solve, DiscontinuousGalerkinScheme
  export fluxρ, myblock, apply_boundary
end
