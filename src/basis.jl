struct Basis{T}
  order::Int
  nodes::Vector{T}
  weights::Vector{T}
  polynomials::Vector{Poly}
  φₕ::Matrix{T}
  ψₕ::Matrix{T}
  dφₕ::Matrix{T}
end

function legendre{T<:Number}(n, ::Type{T}=Float64, var=:x)
    if n==0
        return Poly{T}([one(T)], var)
    elseif n==1
        return Poly{T}([zero(T), one(T)], var)
    end
    px = Poly{T}([zero(T), one(T)], var)
    p0 = Poly{T}([one(T)], var)
    p1 = px
    for i = 2:n
        p2 = ( (2i-1)*px*p1 - (i-1)*p0 ) / i
        p0 = p1
        p1 = p2
    end
    return p1
end

function legendre_basis{T<:Number}(order, ::Type{T}=Float64)
  nodes, weights = gausslegendre(order+1)
  φₕ = zeros(T,order+1,order+1)
  dφₕ = zeros(T,order+1,order+1)
  # TODO: # of faces depend on dimensions
  ψₕ = zeros(T,2,order+1)
  polynomials = Vector{Poly}(order+1)
  for n = 0:order
    p = legendre(n, T)
    dp = polyder(p)
    polynomials[n+1] = p
    # Eval interior nodes
    φₕ[:,n+1] = polyval(p, nodes)
    dφₕ[:,n+1] = polyval(dp, nodes)
    # Eval faces nodes
    ψₕ[:,n+1] = polyval(p, [-1.0,1.0])
  end
  Basis{T}(order,nodes,weights,polynomials,φₕ,ψₕ,dφₕ)
end

function quadrature_to_interval(x,a::Tuple)
   0.5*(a[2]-a[1])*x + 0.5*(a[2]+a[1])
end

function project_function(f, basis, interval::Tuple)
  nodes = quadrature_to_interval(basis.nodes, interval)
  f_val = f.(nodes)
  function model(x,p)
    result = zeros(x)
    for i in 1:size(p,1)
      result.+=p[i]*polyval(basis.polynomials[i], x)
    end
    result
  end
  p0 = zeros(eltype(basis.nodes),basis.order+1)
  curve_fit(model, basis.nodes, f_val, p0)
end

#TODO: Dispatch on different basis types
"Get mass matrix"
function get_local_mass_matrix{T}(basis::Basis{T}, mesh)
  diagnal = zeros(T, basis.order+1)
  diagnal[:] = 2.0/(2*(0:basis.order)+1)
  M = Vector{T}(basis.order+1)
  m = diagm(diagnal)
  for k in 1:(basis.order+1)
    M[k] = mesh.cell_dx[k]/2.0*m
  end
  return M
end

"Get mass matrix inverse"
function get_local_inv_mass_matrix{T}(basis::Basis{T}, mesh)
  diagnal = zeros(T, basis.order+1)
  diagnal[:] = (2*(0:basis.order)+1) / 2.0
  M_inv = Vector{Matrix{T}}(basis.order+1)
  M = diagm(diagnal)
  for k in 1:(basis.order+1)
    M_inv[k] = 2.0./mesh.cell_dx[k]*M
  end
  return M_inv
end
