struct DGSolution{T,N,uType,usType,xType,xfType,uType2,DType,tType,rateType,P,B,A,IType} <: AbstractTimeseriesSolution{T,N}
  u::uType
  uₛ::usType
  nodes::xType
  face_nodes::xfType
  NC::Int
  u_analytic::uType2
  errors::DType
  t::tType
  k::rateType
  prob::P
  basis::B
  alg::A
  interp::IType
  dense::Bool
  tslocation::Int
  retcode::Symbol
end

function reconstruct_u(u::Matrix, φ::Matrix, NC::Int)
  uh = myblock(φ,NC)*u
  NN = size(φ,1); Nx = size(u,2)
  uₕ = zeros(eltype(u), NN*Nx,NC)
  for j in 1:NC
    uₕ[:,j] = uh[(j-1)*NN+1:j*NN,:][:]
  end
  return uₕ
end

function build_solution{T,N,uType,uType2,DType,tType,
rateType,P,A,IType}(
  ode_sol::ODESolution{T,N,uType,uType2,DType,tType,
  rateType,P,A,IType}, basis, problem, NC)
  @unpack u_analytic, errors, t, k,alg,interp,
  dense,tslocation,retcode = ode_sol
  mesh = problem.mesh
  #Use first value to infere types
  k = 1
  NN = size(basis.φₕ,1)
  NM = size(basis.ψₕ,1)
  Nx = size(ode_sol.u[k],2)
  uₕ = reconstruct_u(ode_sol.u[k], basis.φₕ, NC)
  u = Vector{typeof(uₕ)}(0)
  push!(u, copy(uₕ))

  uf = reconstruct_u(ode_sol.u[k], basis.ψₕ, NC)
  uₛ = Vector{typeof(uf)}(0)
  push!(uₛ, copy(uf))

  # Save the rest of the iterations
  for k in 2:size(ode_sol.u,1)
    uₕ = reconstruct_u(ode_sol.u[k], basis.φₕ, NC)
    push!(u, copy(uₕ))
    uf = reconstruct_u(ode_sol.u[k], basis.ψₕ, NC)
    push!(uₛ, copy(uf))
  end
  xg = zeros(NN,Nx)
  xf = zeros(NM,Nx)
  for i in 1:size(xg,2)
    b = mesh.cell_faces[i+1]; a=mesh.cell_faces[i]
    xg[:,i] = 0.5 * (b - a) * basis.nodes + 0.5 * (b + a)
    xf[:,i] = [a,b]
  end
  xg = xg[:]
  xf = xf[:]
  #TODO: Actually compute errors
  DGSolution{T,N,typeof(u),typeof(uₛ),typeof(xg),typeof(xf),typeof(u_analytic),typeof(errors),
  typeof(t),typeof(k),typeof(problem),typeof(basis),typeof(alg),typeof(interp)}(u, uₛ, xg, xf,
  NC,u_analytic, errors, t, k,problem,basis,alg,interp,dense,tslocation,retcode)
end
@def interp_common begin
  NN = size(sol.basis.φₕ,1)
  if typeof(uₘ) <: Vector
    u = Vector{Matrix{eltype(uₘ[1])}}(0)
    for M in uₘ
      uₕ = reconstruct_u(M, sol.basis.φₕ, sol.NC)
      push!(u, copy(uₕ))
    end
    return u
  else
    uₕ = reconstruct_u(uₘ, sol.basis.φₕ, sol.NC)
    return uₕ
  end
end
function (sol::DGSolution)(t,deriv::Type=Val{0};idxs=nothing)
  uₘ = sol.interp(t,idxs,deriv)
  @interp_common
end
function (sol::DGSolution)(v,t,deriv::Type=Val{0};idxs=nothing)
  uₘ = sol.interp(v,t,idxs,deriv)
  @interp_common
end
