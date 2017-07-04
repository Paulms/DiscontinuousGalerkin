abstract type DGProblem end
struct DGScalar1DProblem  <: DGProblem
  f0::Function
  left_boundary::Symbol
  right_boundary::Symbol
end

function DGScalar1DProblem(f0; left_b = :Periodic, right_b = :Periodic)
  DGScalar1DProblem(f0, left_b, right_b)
end
