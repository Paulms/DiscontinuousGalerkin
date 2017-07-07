abstract type DGProblem end
struct DG1DProblem{MeshType,tType,T1}  <: DGProblem
  f0::Function
  f::Function
  max_w_speed::Function
  mesh::MeshType
  tspan::Tuple{tType,tType}
  cfl::T1
  left_boundary::Symbol
  right_boundary::Symbol
end

function DG1DProblem(f0,f,max_w_speed, mesh, tspan,cfl; left_b = :Periodic, right_b = :Periodic)
  DG1DProblem{typeof(mesh), eltype(tspan), typeof(cfl)}(f0,f, max_w_speed, mesh, tspan,cfl,left_b, right_b)
end
