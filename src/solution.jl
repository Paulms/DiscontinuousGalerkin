struct DGSolution{T,N,uType,usType,xType,xfType,uType2,DType,tType,rateType,P,B,A,IType} <: AbstractTimeseriesSolution{T,N}
  u::uType
  uₛ::usType
  nodes::xType
  face_nodes::xfType
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

function build_solution{T,N,uType,uType2,DType,tType,
rateType,P,A,IType}(
  ode_sol::ODESolution{T,N,uType,uType2,DType,tType,
  rateType,P,A,IType}, basis, problem)
  @unpack u_analytic, errors, t, k,alg,interp,
  dense,tslocation,retcode = ode_sol
  mesh = problem.mesh
  #Use first value to infere types
  NC = size(ode_sol.u[1],3)

  u = Vector{Matrix{eltype(ode_sol.u[1])}}(0)
  uₛ = Vector{Matrix{eltype(ode_sol.u[1])}}(0)

  for k in 1:size(ode_sol.u,1)
    uh = zeros(eltype(ode_sol.u[k]),size(basis.φₕ,1)*size(ode_sol.u[k],2), NC)
    uf = zeros(eltype(ode_sol.u[k]),size(basis.ψₕ,1)*size(ode_sol.u[k],2), NC)
    for j in 1:NC
      uₕ = basis.φₕ*ode_sol.u[k]
      us = basis.ψₕ*ode_sol.u[k]
      uh[:,j] = uₕ
      uf[:,j] = us
    end
    push!(u, uₕ)
    push!(uₛ, uf)
  end
  xg = zeros(ode_sol.u[1][:,:,1])
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
  u_analytic, errors, t, k,problem,basis,alg,interp,dense,tslocation,retcode)
end
function (sol::DGSolution)(t,deriv::Type=Val{0};idxs=nothing)
  uₘ = sol.interp(t,idxs,deriv)
  if typeof(uₘ) <: Vector
    NC = size(uₘ[1],3)
    u = Vector{Matrix{eltype(uₘ[1])}}(0)
    for k in 1:size(uₘ,1)
      uh = zeros(eltype(ode_sol.u[k]),size(basis.φₕ,1)*size(ode_sol.u[k],2), NC)
      uf = zeros(eltype(ode_sol.u[k]),size(basis.ψₕ,1)*size(ode_sol.u[k],2), NC)
      for j in 1:NC
        uₕ = basis.φₕ*ode_sol.u[k]
        us = basis.ψₕ*ode_sol.u[k]
        uh[:,j] = uₕ
        uf[:,j] = us
      end
      push!(u, uₕ)
      push!(uₛ, uf)
    end
    return [(sol.basis.φₕ*M[:,:,j])[:] for j in size(M,3) for M in uₘ]
  else
    return (sol.basis.φₕ*uₘ)[:]
  end
end
function (sol::DGSolution)(v,t,deriv::Type=Val{0};idxs=nothing)
  uₘ = sol.interp(v,t,idxs,deriv)
  if typeof(uₘ) <: Vector
    return [(sol.basis.φₕ*M)[:] for M in uₘ]
  else
    return (sol.basis.φₕ*uₘ)[:]
  end
end
