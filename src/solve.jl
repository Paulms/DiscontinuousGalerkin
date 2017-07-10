# Define Scheme type for solve function dispatch
struct DiscontinuousGalerkinScheme
  TimeAlgorithm::OrdinaryDiffEqAlgorithm
  basis::Basis
  riemann_solver::Function
end

"DiscontinuousGalerkinScheme constructor"
function DiscontinuousGalerkinScheme(basis, riemann_solver, TimeIntegrator = SSPRK22())
  DiscontinuousGalerkinScheme(TimeIntegrator, basis, riemann_solver)
end

"Solve scalar 1D conservation laws problems"
function solve{MeshType}(
  problem::DG1DProblem{MeshType},
  alg::DiscontinuousGalerkinScheme;  kwargs...)
  # Unpack some useful variables
  @unpack basis, TimeAlgorithm, riemann_solver = alg
  mesh = problem.mesh
  N = mesh.N
  NC = size(problem.f0(mesh.cell_faces[1]),1)
  NN = basis.order+1
  #Assign Initial values (u0 = φₕ⋅u0ₘ)
  u0ₘ = zeros(NN*NC, N)
  for i = 1:N
    for j = 1:NC
      value = project_function(problem.f0,basis,(mesh.cell_faces[i],mesh.cell_faces[i+1]); component = j)
      u0ₘ[(j-1)*NN+1:NN*j,i] = value.param
    end
  end

  #build inverse of mass matrix
  M_inv = get_local_inv_mass_matrix(basis, mesh)

  #Time loop
  #First dt
  u0ₕ = reconstruct_u(u0ₘ, basis.φₕ, NC)
  dt = update_dt(u0ₕ, problem, basis.order)
  # Setup time integrator
  semidiscretef(t,u,du) = residual!(du, u, basis, problem, riemann_solver, M_inv,NC)
  prob = ODEProblem(semidiscretef, u0ₘ, problem.tspan)
  timeIntegrator = init(prob, TimeAlgorithm;dt=dt, kwargs...)
  @inbounds for i in timeIntegrator
    uₕ = reconstruct_u(timeIntegrator.u, basis.φₕ, NC)
    dt = update_dt(uₕ, problem, basis.order)
    set_proposed_dt!(timeIntegrator, dt)
  end
  if timeIntegrator.sol.t[end] != problem.tspan[end]
    savevalues!(timeIntegrator)
  end
  return build_solution(timeIntegrator.sol,basis,problem, NC)
end

"Reconstruc solution from basis space"
function reconstruct_u(u::Matrix, φ::Matrix, NC::Int)
  uh = myblock(φ,NC)*u
  NN = size(φ,1); Nx = size(u,2)
  uₕ = zeros(eltype(u), NN*Nx,NC)
  for j in 1:NC
    uₕ[:,j] = uh[(j-1)*NN+1:j*NN,:][:]
  end
  return uₕ
end

"Update dt based on CFL condition"
function update_dt(u, problem::DG1DProblem, order)
  ν = problem.max_w_speed(u)
  cfl = problem.cfl
  dx = max(problem.mesh.cell_dx)
  dx * cfl / (ν * (2 * order + 1))
end

"Apply boundary conditions on scalar problems"
function apply_boundary(u, problem::DG1DProblem)
  if problem.left_boundary == :Periodic
      u[:,1] = u[:,end-1]
  elseif  problem.left_boundary == :ZeroFlux
      u[:,1] = u[:,2]
  end
  if problem.right_boundary == :Periodic
      u[:,end] = u[:,2]
  elseif  problem.left_boundary == :ZeroFlux
      u[:,end] = u[:,end-1]
  end
end

"build a block diagonal matrix by repeating a matrix N times"
function myblock(A::Matrix,N::Int)
  M = size(A,1)
  Q = size(A,2)
  B = zeros(eltype(A),M*N,Q*N)
  for i = 1:N
    B[(i-1)*M+1:i*M,(i-1)*Q+1:i*Q] = A
  end
  B
end

"Calculates residual for DG method
  Inputs:
    H = matrix used to store residuals
    u = coefficients of current finite solution approx.
  residual = M_inv*(Q+F)
         where M_inv is the inverse mass matrix
              Q is the edge fluxes
              F is the interior flux"
@def scalar_1D_residual_common begin
  #Add ghost cells
  uₘ = hcat(zeros(u[:,1]),u,zeros(u[:,1]))

  #Apply boundary conditions TODO: Other boundary types
  apply_boundary(uₘ, problem)

  #Reconstruct u in finite space: uₕ(ξ)
  uₕ = myblock(basis.φₕ,NC)*uₘ
  F = zeros(uₕ)
  Fₕ = zeros(uₕ)
  NN = basis.order+1
  for k in 1:size(uₕ,2)
    for j in 1:NN
      Fₕ[j:NN:size(uₕ,1),k] = problem.f(uₕ[j:NN:size(uₕ,1),k])
    end
  end
  # Integrate interior fluxes ∫f(uₕ)φ'(ξ)dξ
  F = A_mul_B!(F,myblock(basis.dφₕ.*basis.weights,NC)',Fₕ)

  # Evaluate edge fluxes
  uₛ = myblock(basis.ψₕ,NC)*uₘ
  q = zeros(T,NC,size(uₛ,2)-1)
  for i = 1:(size(uₛ,2)-1)
    ul = uₛ[2:2:size(uₛ,1),i]; ur = uₛ[1:2:size(uₛ,1),i+1]
    q[:,i] = riemann_solver(ul,ur)
  end
  Q = zeros(F)
  NN = size(basis.φₕ,2)
  for l in 1:NN
    for j in 1:NC
      Q[(j-1)*NN+l,2:end-1] = q[j,2:end] + (-1)^l*q[j,1:end-1]
    end
  end
end
function residual!{T}(H, u, basis::Basis{T}, problem::DGProblem, riemann_solver, M_inv,NC)
  @scalar_1D_residual_common
  H[:,:] = F[:,2:(end-1)]-Q[:,2:(end-1)]
  #Calculate residual
  for k in 1:mesh.N
    H[:,k] = myblock(M_inv[k],NC)*H[:,k]
  end
  H
end

"Efficient residual computation for uniform problems"
function residual!{T, MeshType<:DGU1DMesh}(H, u, basis::Basis{T}, problem::DG1DProblem{MeshType}, riemann_solver, M_inv,NC)
  @scalar_1D_residual_common
  #Calculate residual
  A_mul_B!(H,myblock(M_inv,NC),F[:,2:(end-1)]-Q[:,2:(end-1)])
end
