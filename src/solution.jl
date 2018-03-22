struct DGSolution{T,N,uType,xType,uType2,DType,tType,rateType,P,B,A,IType} <: AbstractTimeseriesSolution{T,N}
  u::uType
  nodes::xType
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

function build_solution{T,N,uType,uType2,DType,tType,
rateType,P,A,IType}(
  ode_sol::ODESolution{T,N,uType,uType2,DType,tType,
  rateType,P,A,IType}, xg, basis, problem, NC)
  @unpack u_analytic, errors, t, k,alg,interp,
  dense,tslocation,retcode = ode_sol
  mesh = problem.mesh
  #Use first value to infere types
  k = 1
  NN = size(basis.φ,1)
  Nx = size(ode_sol.u[k],2)
  uh = flat_u(ode_sol.u[k], basis.order, NC)
  u = Vector{typeof(uh)}(0)
  push!(u, copy(uh))

  # Save the rest of the iterations
  for k in 2:size(ode_sol.u,1)
    uh = flat_u(ode_sol.u[k], basis.order, NC)
    push!(u, copy(uh))
  end
  xg = xg[:]
  #TODO: Actually compute errors
  DGSolution{T,N,typeof(u),typeof(xg),typeof(u_analytic),typeof(errors),
  typeof(t),typeof(k),typeof(problem),typeof(basis),typeof(alg),typeof(interp)}(u, xg,
  NC,u_analytic, errors, t, k,problem,basis,alg,interp,dense,tslocation,retcode)
end

function interp_common(sol, u)
  NN = size(sol.basis.φ,1)
  if typeof(u) <: Vector
    u = Vector{Matrix{eltype(uₘ[1])}}(0)
    for M in u
      uₕ = flat_u(M, sol.basis.order, sol.NC)
      push!(u, copy(uₕ))
    end
    return u
  else
    uₕ = flat_u(u, sol.basis.order, sol.NC)
    return uₕ
  end
end
function (sol::DGSolution)(t,deriv::Type=Val{0};idxs=nothing)
  u = sol.interp(t,idxs,deriv)
  return interp_common(sol, u)
end
function (sol::DGSolution)(v,t,deriv::Type=Val{0};idxs=nothing)
  u = sol.interp(v,t,idxs,deriv)
  return interp_common(sol, u)
end

function save_csv(sol::DGSolution, file_name::String; idx = -1)
  if !endswith(file_name,".csv")
    file_name = "$file_name.csv"
  end
  if idx == -1
    writedlm(file_name, hcat(cell_centers(sol.prob.mesh),sol.u[end]), ',')
  else
    writedlm(file_name, hcat(cell_centers(sol.prob.mesh),sol.u[idx]), ',')
  end
end
