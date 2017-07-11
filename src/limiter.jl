function minmod(a,b,c)
  if (a > 0 && b > 0 && c > 0)
    min(a,b,c)
  elseif (a < 0 && b < 0 && c < 0)
    max(a,b,c)
  else
    zero(a)
  end
end

"Apply slope limiter πᴺ to u assuming u an Nth polynomial"
function SlopLimitN(uₘ::Matrix, basis::DiscontinuousGalerkin.PolynomialBasis, mesh)
  # Compute cell averages
  uh = basis.L2M*uₘ
  uh = uh ./ [factorial(i) for i in 0:basis.order]
  ū = [0.5*(polyval(polyint(Poly(uh[:,i])),1)-polyval(polyint(Poly(uh[:,i])),-1)) for i in size(uh,2)]

  ulimit = uₘ
  # find end values of each element
  ue1 = uₘ[1,:];ue2=uₘ[end,:]

  #find cell averages
  vk = ū; vkm1=[ū[1],ū[1:end-1]];vkp1=[ū[2:end],ū[end]]

  #apply reconstruction to find elements in need of limiting
  ve1 = vk - minmod.(vk-ue1, vk -vkm1,vkp1-vk)
  ve2 = vk + minmod.(ue2-vk,vk-vkm1,vkp1-vk)
  tol = 1e-8
  idx = (abs.(ve1-ue1)>tol)|(abs.(ve2-ue2)>tol)
  if (!isempty(idx))
    h2 = 2/diff(mesh.cell_faces)
    x1 = zeros(basis.order+1,mesh.N)
    for k in 1:mesh.N
      x1[:,k] = reference_to_interval(basis.nodes, (mesh.cell_faces[k],mesh.cell_faces[k+1]))
    end
    x0 = mesh.cell_center
    ux = 0.5*h2.*((basis.dφₕ*uₘ)'*basis.weights)
    for j in 1:size(ulimit,1)
      ulimit[j,idx] = vk[idx]+(x1[j,idx]-x0[idx])*
      minmod.(ux[idx],h2[idx].*(vkp1[idx]-vk[idx]),h2[idx].*(vk[idx]-vkm1[idx]))
    end
  end
  return ulimit
end
