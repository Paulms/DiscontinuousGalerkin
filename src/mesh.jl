struct DG1DMesh{T}
  N ::Int
  cell_dx :: Vector{T}
  cell_center :: Vector{T}
  cell_faces :: Vector{T}
end

"Build an uniform 1D mesh"
function Uniform1DMesh{T<:Number}(N::Int,xinit::T,xend::T)
  L = xend - xinit
  dx = L/N
  xx = [i*dx+dx/2+xinit for i in 0:(N-1)]
  faces = [xinit + dx*i for i in 0:N]
  DG1DMesh{T}(N,ones(T,N)*dx,xx,faces)
end
