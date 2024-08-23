using Distributed
using DistributedArrays
using DistributedArrays.SPMD
@everywhere using Distributed
@everywhere using DistributedArrays
@everywhere using DistributedArrays.SPMD
using StaticArrays
using Base.Threads
using Printf

includet("aux_functions_poisson_divergence.jl")

includet("aux_functions_particles.jl")

includet("aux_functions_E-B.jl")

includet("aux_functions_shapes.jl")

includet("aux_functions_Interpolate.jl")

includet("aux_functions_Integrate.jl")

includet("aux_functions_density_charge.jl")

includet("aux_functions_grid.jl")

includet("aux_functions_RHS.jl")

includet("aux_functions_IO.jl")

includet("choques_utils.jl")

includet("aux_functions_RHS.jl")


function compare_electric_field_constraints(v,j,par_grid, par_evolv, run_name, save_plots)
  N, L, J, dx, order = par_grid
  t_i, t_f, M, M_g, dt = par_evolv
  ρ_f = zeros(J)
  E_f = zeros(J)
  E_i = v[2N+1:end,1]
  ϕ_f = zeros(J)
  #S_f = zeros(J)

  @assert j <= M_g


  get_density!(v[:,j], ρ_f, par_grid,0.0)
  get_ϕ!(ϕ_f, ρ_f .+ 1, 2π/L)
  get_E_from_ϕ!(ϕ_f,E_f,dx)

  if !remote_server || plots
      plt = plot(x,E_f,label="from final density", ls=:dash, lw=2)
      plot!(x,E_i,label="E_initial")
      plot!(x,v[2N+1:end,j], label="E_final"
      )
      if save_plots
      png("Images/" * run_name * "_Efield")
      end
  end
return plt
end



function get_temperature(u,N;m=1)
  return m*var(u[N+1:2N])
end

function get_temperature_rel(u,N;m=1)
  return m*var(p2v.(u[N+1:2N]))
end

"""
Chequeada en ini_dat_v2.ipynb
"""
function get_temperature_rel_D(u,N,D;m=1)
  sv = zeros(D)
  sv2 = zeros(D)
  for d in 1:D
    sv[d] = sum(u[D+d:2D:N*2D])
    sv2[d] = sum(u[D+d:2D:N*2D].^2)
  end
  return m*(sum(sv2) - sv'*sv)/N/D
end



function get_energy_rel(u,p)
  L, N, J = p
  dx = L/J
  n0 = N/L
  energy_K = 0.0
  energy_E = 0.0
  for i in 1:N
    get_momenta!
    energy_K = energy_K + (sqrt(1+u[N+i]^2) - 1)
  end
  for j in 1:J
    energy_E = energy_E + u[2N+j]^2
  end
  # return energy_K,  dx * energy_E / 2 * n0
  return energy_K / n0,  dx * energy_E / 2 # normalized version
end

function get_energy_rel(u,par; m=1)
  Box, N, J= par
  dx = differentials(Box,J)
  D = length(dx)
  n0 = N /volume(Box)
  energy_K = 0.0
  energy_E = 0.0
  p = zeros(D)
  for i in 1:N
    get_momenta!(p,i,u)
    energy_K = energy_K + (sqrt(m^2+p'*p/m^2) - m)
  end
  
    #Fu = view(u,4N+1, 4N+3*prod(J))

    energy_E = u[4N+1:end]'*u[4N+1:end]

    #E = get_E(u,N,J)
    #B = get_B(u,N,J)

    #energy_E = sum.(E'*E+B'*B)
  
  # return energy_K,  dx * energy_E / 2 * n0
  return energy_K / n0,  prod(dx) * energy_E / 2 # normalized version
end



