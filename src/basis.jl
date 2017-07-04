struct basis{T}
  order::Int
  nodes::Vector{T}
  weights::Vector{T}
  polynomials::Vector{Poly}
  φₕ::Matrix{T}
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
  polynomials = Vector{Poly}(order+1)
  for n = 0:order
    p = legendre(n, T)
    dp = polyder(p)
    polynomials[n+1] = p
    φₕ[:,n+1] = polyval(p, nodes)
    dφₕ[:,n+1] = polyval(dp, nodes)
  end
  basis{T}(order,nodes,weights,polynomials,φₕ,dφₕ)
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
