using Distributed
using DistributedArrays
using DistributedArrays.SPMD
@everywhere using Distributed
@everywhere using DistributedArrays
@everywhere using DistributedArrays.SPMD
using StaticArrays
using Base.Threads
using Printf

include("aux_functions_poisson_divergence.jl")

include("aux_functions_particles.jl")

include("aux_functions_E-B.jl")

include("aux_functions_shapes.jl")

include("aux_functions_Interpolate.jl")

include("aux_functions_Integrate.jl")

include("aux_functions_density_charge.jl")

include("aux_functions_grid.jl")

include("aux_functions_RHS.jl")


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
  n0 = N/volume(Box)
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


#=
################################################################################################
Some functions to handle data
################################################################################################
=#
function get_averages(v,par_grid,par_evolv, par_f)
  (N, L, J, dx, order) = par_grid
  (t_i, t_f, M, M_g, dt) = par_evolv
  (θ, nm, κ) = par_f
  Energy_K = zeros(M_g)
  Energy_E = zeros(M_g)
  EField_T = zeros(M_g)
  p_T = zeros(M_g)
  Q_T = zeros(M_g)
  S_T = zeros(M_g)
  #E_E = 0.0
  E_mode = zeros(M_g)
  T = zeros(M_g)
  #P = zeros(J)
  ρ = zeros(J)
  S = zeros(J)

  for j in 1:M_g
      (Energy_K[j],Energy_E[j]) = get_energy_rel(v[:,j],(L,N,J))
      EField_T[j] = sum(v[2N+1:end,j])*dx
      p_T[j] = sum(v[N+1:2N])*dx
      get_density!(v[:,j], ρ, par_grid,0.0)
      get_current_rel!(v[:,j], S, par_grid)
      Q_T[j] = get_total_charge(ρ,(J, dx))/L # we divide by L because the density is 1, so the total charge is L, this way we compare with 1.
      S_T[j] = sum(S)/N/Q_T[j]
      T[j] = var(v[N+1:2N,j])
      E_mode = abs(rfft(v[2N+1:end,j])[nm+1])
  end
  return Energy_K, Energy_E, EField_T, p_T, Q_T, S_T, T, E_mode
end

function get_averages_threads(v,par_grid,par_evolv, par_f)
  
  (N, L, J, dx, order) = par_grid
  (t_i, t_f, M, M_g, dt) = par_evolv
  (θ, nm, κ) = par_f

  TS = zeros(J, nthreads())
  Tn = zeros(J, nthreads())
  

  par_density = (par_grid, Tn)
  par_current = (par_grid, TS)

  Energy_K = zeros(M_g)
  Energy_E = zeros(M_g)
  EField_T = zeros(M_g)
  p_T = zeros(M_g)
  Q_T = zeros(M_g)
  S_T = zeros(M_g)
  #E_E = 0.0
  E_mode = zeros(M_g)
  T = zeros(M_g)
  #P = zeros(J)
  ρ = zeros(J)
  S = zeros(J)

  for j in 1:M_g
      (Energy_K[j],Energy_E[j]) = get_energy_rel(v[:,j],(L,N,J))
      EField_T[j] = sum(v[2N+1:end,j])
      p_T[j] = sum(v[N+1:2N])*dx
      get_density_threads!(v[:,j], ρ, par_density,0.0)
      get_current_rel_threads!(v[:,j], S, par_current)
      Q_T[j] = get_total_charge(ρ,(J, dx)) / L # we divide by L because the density is 1, so the total charge is L, this way we compare with 1.
      S_T[j] = sum(S)/N/Q_T[j]
      T[j] = var(v[N+1:2N,j])
      E_mode = abs(rfft(v[2N+1:end,j])[nm+1])
  end
  return Energy_K, Energy_E, EField_T, p_T, Q_T, S_T, T, E_mode
end


function get_local_averages_threads(u,par_grid, par_f)
  
  (N, L, J, dx, order) = par_grid
  (θ, nm, κ) = par_f

  TS = zeros(J, nthreads())
  Tn = zeros(J, nthreads())
  
  par_density = (par_grid, Tn)
  par_current = (par_grid, TS)

  ρ = zeros(J)
  S = zeros(J)
  Efield = zeros(J)

  Energy_K, Energy_E = get_energy_rel(u,(L,N,J))
  Efield = u[2N+1:end]
  EField_T = sum(u[2N+1:end])*dx
  p_T = sum(u[N+1:2N])*dx
  get_density_threads!(u[:], ρ, par_density, 0.0)
  get_current_rel_threads!(u[:], S, par_current)
  Q_T = get_total_charge(ρ,(J, dx)) / L # we divide by L because the density is 1, so the total charge is L, this way we compare with 1.
  S_T = sum(S)/N/Q_T
  #T = var(u[N+1:2N])
  T = get_temperature_rel(u,N)
  E_mode = abs(rfft(u[2N+1:end])[nm+1])

  return ρ, S, Efield, Energy_K, Energy_E, EField_T, p_T, Q_T, S_T, T, E_mode
end

function load_averages(file_name, j, par_grid, pars_f)
    ρ, S, Efield, Energy_K, Energy_E, EField_T, p_T, Q_T, S_T, T, E_mode = get_local_averages_threads(u,par_grid, pars_f)
    tiempo = @sprintf("%05d", j)
    jldopen(file_name, "a+") do file
        file["n_$(tiempo)"] = ρ
        file["S_$(tiempo)"] = S
        file["Efield_$(tiempo)"] = Efield
        file["Energy_E_$(tiempo)"] = Energy_E
        file["Energy_K_$(tiempo)"] = Energy_K
        file["EField_T_$(tiempo)"] = EField_T
        file["p_T_$(tiempo)"] = p_T
        file["Q_T_$(tiempo)"] = Q_T
        file["S_T_$(tiempo)"] = S_T
        file["T_$(tiempo)"] = T
        file["E_mode_$(tiempo)"] = E_mode
    end
end

function retrieve_average_data(data, par_grid, par_evolv; M_last=nothing)
  (N, L, J, dx, order) = par_grid
  (t_i, t_f, M, M_g, dt) = par_evolv
  #v = zeros(2N+J,M_g)
  n_t = zeros(J,M_g)
  S_t = zeros(J,M_g)
  Efield_t = zeros(J,M_g)
  EField_T = zeros(M_g)
  Energy_K = zeros(M_g)
  Energy_E = zeros(M_g)
  E_mode = zeros(M_g)
  p_T = zeros(M_g)
  Q_T = zeros(M_g)
  S_T = zeros(M_g)
  T = zeros(M_g)
  if M_last !== nothing # if we gave some values, then use it.
    M_g = M_last
  end
  for j in 1:M_g
      tiempo = @sprintf("%05d", j)
      n_t[:,j] = data["n_$(tiempo)"]
      S_t[:,j] = data["S_$(tiempo)"]
      Efield_t[:,j] = data["Efield_$(tiempo)"]
      Energy_K[j] = data["Energy_K_$(tiempo)"]
      Energy_E[j] = data["Energy_E_$(tiempo)"]
      EField_T[j] = data["EField_T_$(tiempo)"]
      p_T[j] = data["p_T_$(tiempo)"]
      Q_T[j] = data["Q_T_$(tiempo)"]
      S_T[j] = data["S_T_$(tiempo)"]
      T[j] = data["T_$(tiempo)"]
      E_mode[j] = data["E_mode_$(tiempo)"]
  end
  return n_t, S_t, Efield_t, (Energy_E,  Energy_K, EField_T, p_T, Q_T, S_T, T, E_mode)
end

function retrieve_data(data, par_grid, par_evolv)
  (N, L, J, dx, order) = par_grid
  (t_i, t_f, M, M_g, dt) = par_evolv
  v = zeros(2N+J,M_g)
  for j in 1:M_g
      tiempo = @sprintf("%05d", j)
      v[:,j] = data["u/u_$tiempo"]
  end
  return v
end

function retrieve_meta_data(file_name::String)
  data = load(file_name)
  run_name = data["run_name"]
  par_grid = data["par_grid"]
  par_evolv = data["par_evolv"]
  par_f = data["p_Ini"]
  (N, L, J, dx, order) = par_grid
  (t_i, t_f, M, M_g, dt) = par_evolv
  #@show (θ, nm, κ) = par_f
  n0 = N/L
  x = [(i-1)*dx for i in 1:J]
  dT = dt * (M-1) / (M_g-1)
  t_series = [(i-1)*dT for i in 1:M_g]
  return data, run_name, par_grid, par_evolv, par_f, n0, x, t_series
end


########################################################################
# plotting functions
########################################################################

function plot_averages(averages, t_series, N, run_name, save_plots)
  Energy_K, Energy_E, EField_T, p_T, Q_T, S_T = averages
  M_last = length(t_series)
  E1 = Energy_K[1] + Energy_E[1]
  plt = plot(layout=(2,2), size=(800,600))
  plot!(subplot=1, t_series, Energy_K[1:M_last] .- Energy_K[1], label="Energy_K")
  plot!(subplot=1, t_series, Energy_E[1:M_last] .- Energy_E[1], label="Energy_E")
  #plot!(subplot=1, t_series, Energy_K, label="Energy_K")
  #plot!(subplot=1, t_series, Energy_E[1:400], label="Energy_E")
  plot!(subplot=2, t_series, (Energy_K + Energy_E)[1:M_last] ./ E1 .- 1.0, label="Total Energy")
  plot!(subplot=3, t_series, Q_T[1:M_last] .- 1, label="charge")
  plot!(subplot=4, t_series, S_T[1:M_last], label="Total Current", legend=:topleft)
  if save_plots
      png("Images/"  * run_name * "_total_run")
  end
  return plt
end

function plot_energies(Energy_K, Energy_E, t_series, run_name, save_plots)
  M_last = length(t_series)
  E1 = Energy_K[1] + Energy_E[1]
  plt = plot(t_series[2:end], abs.(Energy_K[2:M_last] .- Energy_K[1]), title = "Energy conservation (order = $(order))", label = "Kinetic Energy"
    #, legend = :outertopright
    , legend = :bottomright, ls=:dash)
    plot!(t_series[2:M_last], abs.(Energy_E[2:M_last] .- Energy_E[1]), label = "|Electric Energy|", ls=:dot)
    plot!(t_series[2:M_last], abs.(Energy_K[2:M_last] + Energy_E[2:M_last] .- E1) ./ E1
    , yscale=:log10
    #, xscale=:log10
    , label = "Total Energy / Initial Energy -1 ")
    if save_plots
    png("Images/" * run_name * "_energy_conservation")
    end
  return plt
end

"""
several fits for the energy:
  The first one is a couple of sinusoidal functions. 
  The second one is a square fit.
  The third one is a more ellaborated square fit.
  """
function energy_fit(t_series, Energy_E, model_e, pe, N_i, N_f, run_name, save_plots; yscale=:identity)
  
  fit_energy = curve_fit(model_e, t_series[N_i:N_f], Energy_E[N_i:N_f], pe);

  plt = Plots.scatter(t_series[N_i:N_f],Energy_E[N_i:N_f]
  , markersize=1
  , title = "Electric Energy decay"
  , label= "Electric Energy"
  , yscale=yscale 
  )

  plot!(t_series[N_i:N_f], model_e(t_series[N_i:N_f],fit_energy.param), ls=:dash
  , markersize = 0.2
  #, xlims=(00,100)
  , label="Fit"
  )
  
  if save_plots
      png("Images/" * run_name * "_energy_fit")
  end
  return fit_energy, plt 
end

function temperature_fit(t_series, T, p_tl001, N_i, N_f, run_name, save_plots)
  @. model_tl001(x,p) = p[1] + p[2]*cos(p[3]*x + p[4])*exp(-p[5]*x) + p[6]*cos(p[7]*x + p[8])*exp(-p[9]*x)
  #p_tl001 = [0.001; 3; 2; 0; 0.1; 0; 1.; 0; 0]
  fit_tl001 = curve_fit(model_tl001, t_series[N_i:N_f], T[N_i:N_f], p_tl001)
  #@. model_tl001_s(x,p) = p[1] + p[2]*cos(p[3]*x + p[4])*exp(-p[5]*x)
  #p_tl001_s = [0.001; 3; 2; 0; 0.0]
  #fit_tl001_s = curve_fit(model_tl001_s, t_series[1:end], T[1:end], p_tl001_s)
  plt = Plots.scatter(t_series[N_i:N_f], T[N_i:N_f], ms = 1)
  plot!(t_series[N_i:N_f], model_tl001(t_series[N_i:N_f],fit_tl001.param))
  
  if save_plots
      png("Images/" * run_name * "_temperature_fit")
  end
  return fit_tl001, plt
end

"""
Plots the values of a matrix as a surface plot.
"""
function plot_matrix(A::Matrix{Float64}; fc=:ocean, linealpha=0.3, fillalpha=0.5, camera=(60,40), title = "")
    default(size=(600,600)
#, fc=:thermal
#, fc=:heat
    , fc=fc
    )
    if !(ndims(A) == 2) 
        error("Array must be 2 dimensional and seems to be of dims = $(ndims(A))")
    end
    (n,m) = size(A)
    x, y = 1:n, 1:m
    z = Surface((x,y)->A[x,y], x, y)
    surface(x,y,z, linealpha = linealpha, fillalpha=fillalpha, display_option=Plots.GR.OPTION_SHADED_MESH, camera=camera, title = title)
end



