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
  problem::DGScalar1DProblem{MeshType},
  alg::DiscontinuousGalerkinScheme;  kwargs...)
  # Unpack some useful variables
  @unpack basis, TimeAlgorithm, riemann_solver = alg
  mesh = problem.mesh
  N = mesh.N
  NC = size(problem.f0(mesh.cell_faces[1]),1)
  #Assign Initial values (u0 = φₕ⋅u0ₘ)
  u0ₘ = zeros(basis.order+1, N, NC)
  for i = 1:N
    for j in 1:NC
      value = project_function(problem.f0,basis,(mesh.cell_faces[i],mesh.cell_faces[i+1]);component = j)
      u0ₘ[:,i,j] = value.param
    end
  end

  #build inverse of mass matrix
  M_inv = get_local_inv_mass_matrix(basis, mesh)

  #Time loop
  #First dt
  dt = update_dt(u0ₘ, problem, basis.order)
  # Setup time integrator
  semidiscretef(t,u,du) = residual!(du, u, basis, problem, riemann_solver, M_inv)
  prob = ODEProblem(semidiscretef, u0ₘ, problem.tspan)
  timeIntegrator = init(prob, TimeAlgorithm;dt=dt, kwargs...)
  @inbounds for i in timeIntegrator
    dt = update_dt(timeIntegrator.u, problem, basis.order)
    set_proposed_dt!(timeIntegrator, dt)
  end
  if timeIntegrator.sol.t[end] != problem.tspan[end]
    savevalues!(timeIntegrator)
  end
  return build_solution(timeIntegrator.sol,basis,problem)
end

"Update dt based on CFL condition"
function update_dt(u, problem::DGScalar1DProblem, order)
  ν = problem.max_w_speed(u)
  cfl = problem.cfl
  dx = max(problem.mesh.cell_dx)
  dx * cfl / (ν * (2 * order + 1))
end

"Apply boundary conditions on scalar problems"
function apply_boundary(u, problem::DGScalar1DProblem)
  if problem.left_boundary == :Periodic
      u[:,1] = u[:,end-1]
  end
  if problem.right_boundary == :Periodic
      u[:,end] = u[:,2]
  end
end

"Calculates residual for DG method
  Inputs:
    H = matrix used to store residuals
    u = coefficients of current finite solution approx.
  residual = M_inv*(Q+F)
         where M_inv is the inverse mass matrix
              Q is the edge fluxes
              F is the interior flux"

"Efficient residual computation for uniform 1D problems"
function residual!{T, MeshType<:DGU1DMesh}(H, u, basis::Basis{T}, problem::DGScalar1DProblem{MeshType}, riemann_solver, M_inv)
  #Add ghost cells
  NC = size(u,3)
  H = zeros(u)
  for j in 1:NC
    uₘ = zeros(T, size(u,1),size(u,2)+2)
    uₘ = hcat(zeros(u[:,1,j]),u[:,:,j],zeros(u[:,1,j]))

    #Apply boundary conditions TODO: Other boundary types
    apply_boundary(uₘ, problem)

    #Reconstruct u in finite space: uₕ(ξ)
    uₕ = basis.φₕ*uₘ

    #Reconstruct u on faces
    uₛ = basis.ψₕ*uₘ

    F = zeros(uₕ)
    Fₕ = zeros(uₕ)
    for k in 1:size(uₕ,2)
      Fₕ[:,k] = problem.f(uₕ[:,k])
    end
    # Integrate interior fluxes ∫f(uₕ)φ'(ξ)dξ
    F = A_mul_B!(F,basis.dφₕ',Fₕ.*basis.weights)

    # Evaluate edge fluxes
    q = zeros(T,size(uₛ,2)-1)
    for i = 1:(size(uₛ,2)-1)
      ul = uₛ[2,i]; ur = uₛ[1,i+1]
      q[i] = riemann_solver(ul,ur)
    end
    Q = zeros(F)
    for l in 1:size(Q,1)
      Q[l,2:end-1] = q[2:end] + (-1)^l*q[1:end-1]
    end
    #Calculate residual
    Hc = @view H[:,:,j]
    A_mul_B!(Hc,M_inv,F[:,2:(end-1)]-Q[:,2:(end-1)])
  end
end
