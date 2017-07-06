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
  k = 1
  uₕ = basis.φₕ*ode_sol.u[k]
  uₕ = uₕ[:]
  u = Vector{typeof(uₕ)}(0)
  push!(u, uₕ)
  uf = basis.ψₕ*ode_sol.u[end]
  xf = zeros(uf)
  uf = uf[:]
  uₛ = Vector{typeof(uf)}(0)
  push!(uₛ, uf)
  # Save the rest of the iterations
  for k in 2:size(ode_sol.u,1)
    uₕ = basis.φₕ*ode_sol.u[k]
    push!(u, uₕ[:])
    uf = basis.ψₕ*ode_sol.u[k]
    push!(uₛ, uf[:])
  end
  xg = zeros(ode_sol.u[1])
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
    return [(sol.basis.φₕ*M)[:] for M in uₘ]
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
