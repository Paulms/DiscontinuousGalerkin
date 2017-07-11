"Apply slope limiter πᴺ to u assuming u an Nth polynomial"
function SlopLimitN(uₘ::Matrix, basis::PolynomialBasis)
  # Compute cell averages
  uh = basis.L2M*uₘ
  uh = uh ./ [factorial(i) for i in 0:basis.order]
  ū = [0.5*(polyval(polyint(Poly(uh[:,i])),1)-polyval(polyint(Poly(uh[:,i])),-1)) for i in size(uh,2)]

  ulimit = u
  # find end values of each element
  ue1 = uₘ[1,:];ue2=uₘ[end,:]

  #find cell averages
  vk = ū; vkm1=[ū[1],ū[1:end-1]];vkp1=[ū[2:end],ū[end]]
end
