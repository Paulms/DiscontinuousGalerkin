# sod shock tube problem
# build mesh
include("../src/DiscontinuousGalerkin.jl")
using DiscontinuousGalerkin
N = 20
mesh = Uniform1DMesh(N,-1.0,1.0)
#Parameters
xdiafragm = 0.0
const γ = 1.4

#define initial condition
function f0(x::Real)
  #u_fis = [ρ, u, p]
  ufis_l = [1.0, 0.0,	1.0]
  ufis_r = [0.125 0.0 0.1]
  if x <= xdiafragm
    ρ,u,p = ufis_l
  else
    ρ,u,p = ufis_r
  end
  return [ρ, ρ*u, 1.0/(γ - 1.0)*p + 0.5*ρ*u*u]
end

#define flux function (generic for euler systems)
function f(u)
  # Fisical variables
  ρ = u[1];v = u[2]/ρ;E = u[3]
  p = (γ - 1)*(E - 0.5*ρ*v*v)
  # Flux
  return [u[2], ρ*v*v + p, (E + p)*v]
end

#define max wave speed (generic euler system)
function max_w_speed(u)
  ρ = u[:,1];v = u[:,2]./ρ;E = u[:,3]
  p = zeros(ρ); wave_speed = zeros(ρ)
  @. p = (γ-1)*(E-0.5*ρ*v*v)
  @. wave_speed = abs(v) + sqrt(γ*p/ρ)
  return maximum(wave_speed)
end

#Setup problem
tspan = (0.0,0.2)
cfl = 0.5
problem = DG1DProblem(f0, f, max_w_speed, mesh, tspan, cfl;
          left_b = :ZeroFlux, right_b = :ZeroFlux)

#Generate basis for Discontinuous Galerkin Scheme
basis=legendre_basis(3)

#Solve problem
u = solve(problem, DiscontinuousGalerkinScheme(basis, rusanov_solver))

#Plot
using Plots
gr()
plot(u)
