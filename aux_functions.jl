using Distributed
using DistributedArrays
using DistributedArrays.SPMD
@everywhere using Distributed
@everywhere using DistributedArrays
@everywhere using DistributedArrays.SPMD
using StaticArrays
using Base.Threads
using Printf

"""
The following function evaluates the electric field on a uniform grid from the electric potential.

    // Calculate electric field from potential
"""
function get_E_from_ϕ!(ϕ, E, dx)
      J = length(E)
      for j in 2:J-1
        E[j] = (ϕ[j-1] - ϕ[j+1]) / 2. / dx
      end
      E[1] = (ϕ[J] - ϕ[2]) / 2. / dx;
      E[J] = (ϕ[J-1] - ϕ[1]) / 2. / dx;
end

function get_E_from_ϕ!(ϕ, E, Box_x::Tuple)
    J = size(E[1])
    D = length(J)
    dx = differentials(Box_x,J)
    for i in 1:D
      for k in 1:J[1]
        for j in 1:J[2]
         # E[j,k][i] = (ϕ[j-1] - ϕ[j+1]) / 2. / dx
      end
    end
  end

# UNFINISHED
end

""" The following routine solves Poisson's equation in 1-D to find the instantaneous electric potential on a uniform grid.

// Solves 1-d Poisson equation:
//    d^u / dx^2 = v   for  0 <= x <= L
// Periodic boundary conditions:
//    u(x + L) = u(x),  v(x + L) = v(x)
// Arrays u and v assumed to be of length J.
// Now, jth grid point corresponds to
//    x_j = j dx  for j = 0,J-1
// where dx = L / J. L / (J-1) in Julia
// Also,
//    kappa = 2 pi / L
"""
function get_ϕ!(ϕ, ρ, κ)
  #V = fill(0.0+im*0.,J) 
  #U = fill(0.0+im*0.,J÷2+1) 
  J = length(ρ)
  # Fourier transform source term
  V = rfft(ρ)

  # Calculate Fourier transform of u

  V[1] =  0.
  for j in  2:(J÷2+1)
    V[j] = - V[j] / (j-1)^2 / κ^2
  end

  # Inverse Fourier transform to obtain u
  ϕ[:] = irfft(V,J)
end

function get_ϕ_D!(ϕ, ρ, Box)
  #V = fill(0.0+im*0.,J) 
  #U = fill(0.0+im*0.,J÷2+1) 
  J = size(ρ)
  D = length(Box)÷2
  Box_array = [i for i in Box]
  κ = 2π./(Box_array[2:2:end] - Box_array[1:2:end-1])
  # Fourier transform source term
  V = rfft(ρ)

  # Calculate Fourier transform of u

  V[1,1] =  0.
  for i in  2:(J[1]÷2+1)
    for j in 2:(J[2]÷2+1)
      V[i,j] = - V[i,j] / ((i-1)^2 / κ[1]^2 + (j-1)^2 / κ[2]^2)
    end
  end

  # Inverse Fourier transform to obtain u
  ϕ[:,:] = irfft(V,J)
end

"""
get_E_direct!(E,ρ_m)
computes the Electric field from the constraint equation div E = ρ 
using E_k = E_{k-1} + \rho_{k-1/2}*dx
then substracts a constant field so that E_T = 0
"""
function get_E_direct!(E,ρ_m,par_grid)
  N, L, J, dx, order = par_grid
  E[1] = 0.0
  E_T = 0.0
  for j in 2:J
    E[j] = E[j-1] + ρ_m[j-1]*dx
    E_T = E_T + E[j]
  end
  return E .= E .- E_T/J
end
"""
Takes out the mass of the grid function so that the sum is now null
"""
function filter_constant!(E)
  J = length(E)
  V = rfft(E)
  V[1] = 0.0 #extract the first component
  E[:] = irfft(V,J)
end

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


""" 
v2p(v;m=1)
Given a 3-velocity computes the momentum
"""
#v2p(v;m=1) = m*v/sqrt(1-v^2)
v2p(v;m=1) = m*v/sqrt(1-v'*v) # this works fine
#v2p(v::Array;m=1) = m*v/sqrt(1-v'*v)
"""
p2v(p;m=1) 
Given a 3-momentum computes the 3-velocity
"""
p2v(p;m=1) = p/sqrt(m^2+p'*p)

""" 
The energy function from momentum 
"""
γ_p(p;m=1) = sqrt(m^2 + p'*p)/m 

""" 
The energy function from 3-velocity
"""
γ(v;m=1) = m/sqrt(1 - v'*v)


"""The routine below evaluates the electron number density on an evenly spaced mesh given the 
instantaneous electron coordinates.

Evaluates electron number density n(1:J) from 
array r(1:N) of electron coordinates.

shift can have the values 0 or 1/2 depending whether we want the density at gridpoint or at midpoints.
"""
function get_density!(u, n, par_grid, shift)
  N, L, J, dx, order = par_grid
  n0 = N/L
  r = view(u,1:N)
  fill!(n,0.0)
  # Evaluate number density.
  for i in 1:N
    @inbounds j, y = get_index_and_y(r[i],J,L)
    y += - shift
    for l in (-order):order 
      @inbounds n[mod1(j + l, J)] += Shape(order, -y + l) / dx / n0; # the dx here is from the different definition from the paper
    end
  end
  return n[:] # return rho directly (we need to subtract 1 in cases where we assume positive particles, but this is done elsewhere.)
end

function get_density_2D!(u, n, par_grid, shift)
  N, Box, J, order = par_grid
  #vol = volume(Box)
  fill!(n,0.0)
  #n0 = N/vol # correct expression but not needed
  D = 2
  if D != length(J) 
    error("dimension mismach")
   end
  n0 = N
  j = [1,1]
  y = [0.0,0.0]
  #@show u
  # Evaluate number density.
  for i in 1:N
    s = (i-1)*2D + 1
    u_r = view(u,s:(s+D-1))
    #@inbounds j, y = get_index_and_y!(j,y,u_r,J,Box)
    @inbounds get_index_and_y!(j,y,u_r,J,Box)
    y .= y .- shift # shift must be the same in all directions!
    for l in (-order):order 
      for m in (-order):order
      #@inbounds n[mod1(j + l, J)] += Shape(order, -y + l) / dx / n0; # the dx here is from the different definition from the paper
      @inbounds n[mod1(j[1] + l, J[1]), mod1(j[2] + m, J[2])] += Shape(order, -y[1] + l) * Shape(order, -y[2] + m)/ n0
      end
    end
  end
  return n[:,:] # return rho directly (we need to subtract 1 in cases where we assume positive particles, but this is done elsewhere.)
end

function get_density_threads!_2D!(u, n, par, shift)
  par_grid, Tn = p
  N, Box, J, order = par_grid
  D = 2
  if D != length(J) 
    error("dimension mismach")
   end
  j = fill(Int64(1),nthreads()) 
  y = fill(Float64(1.0),nthreads())
  Tn .= zeros(Float64)
  n0 = N
  
  n0 = N
  # Evaluate number density.
  @threads for i in 1:N
    s = (i-1)*2D + 1
    u_r = view(u,s:(s+D-1))
    #@inbounds j, y = get_index_and_y!(j,y,u_r,J,Box)
    @inbounds get_index_and_y!(j,y,u_r,J,Box)
    @inbounds get_index_and_y!(j[threadid()], y[threadid()], u_r, J, Box)
    y[threadid()] .= y[threadid()] .- shift # shift must be the same in all directions!
    for l in (-order):order 
      for m in (-order):order
      #@inbounds n[mod1(j + l, J)] += Shape(order, -y + l) / dx / n0; # the dx here is from the different definition from the paper
      @inbounds n[mod1(j[1] + l, J[1]), mod1(j[2] + m, J[2])] += Shape(order, -y[1] + l) * Shape(order, -y[2] + m)/ n0
      end
    end
  end
  return n[:,:] # return rho directly (we need to subtract 1 in cases where we assume positive particles, but this is done elsewhere.)
end

function get_density_threads!(u, n, p, shift)
  par_grid, Tn = p
  N, L, J, dx, order = par_grid
  j = fill(Int64(1),nthreads()) 
  y = fill(Float64(1.0),nthreads())
  Tn .= zeros(Float64)
  n0 = N/L
  # Evaluate number density.
  @threads for i in 1:N
    @inbounds j[threadid()], y[threadid()] = get_index_and_y(u[i], J, L)
    y[threadid()] += - shift
    for l in (-order):-j[threadid()]
      @inbounds Tn[J + j[threadid()] + l, threadid()] += Shape(order, -y[threadid()] + l) 
    end
    for l in max(-order,-j[threadid()]+1):min(order,J-j[threadid()])
      @inbounds Tn[j[threadid()] + l, threadid()] += Shape(order, -y[threadid()] + l) 
    end
    for l in J-j[threadid()]+1:order
      @inbounds Tn[j[threadid()] - J + l, threadid()] += Shape(order, -y[threadid()] + l) 
    end
  end
  n .= zeros(Float64)
  @threads for i in 1:J
    for t in 1:nthreads()
      @inbounds n[i] += Tn[i, t]/dx/n0 # the dx here is from the different definition from the paper
    end
  end
  return n[:]
end


function get_total_charge(ρ,par)
  J, dx = par
  Q = 0.0
  for i in 1:J
      Q += ρ[i]
  end
  return Q * dx
end

function get_total_charge(ρ,Vol)
  Q = sum(ρ)*Vol/length(ρ)
end

"""The routine below evaluates the electron current on an evenly spaced mesh given the instantaneous electron coordinates.
// Evaluates electron number density S(0:J-1) from 
// array r(0:N-1) of electron coordinates.
"""

function get_current_rel!(u, S, par_grid)
  N, L, J, dx, order = par_grid
  r = view(u,1:N)
  p = view(u,N+1:2N) # in the relativistic version we compute p instead of v
  fill!(S,0.0)
  n0 = N/L
  for i in 1:N
    @inbounds j, y = get_index_and_y(r[i],J,L)
    @inbounds v = p2v(p[i]) / dx / n0 # the dx here is from the different definition from the paper
    for l in (-order):order 
      @inbounds S[mod1(j + l, J)] += Shape(order, -y + l) * v;
    end
  end
  return S[:] # allready normalized with n0
end


    
function get_current_rel_2D!(u, S, par_grid)
  N, Box, J, order = par_grid 
  D = 2
   if D != length(J) 
    error("dimension mismach")
   end
  #vol = volume(Box)
  fill!(S,[0.0,0.0])
  #n0 = N/vol # correct expression but not needed
  n0 = N
  for i in 1:N
    s = (i-1)*2D + 1
    r = view(u,s:s+D-1)
    p = view(u,s+D:s+2*D-1) # in the relativistic version we compute p instead of v
    j = [1,1]
    y = [0.0,0.0]
    @inbounds j, y = get_index_and_y!(j,y,r,J,Box)
    #@inbounds v = p2v(p) / vol / n0 # correct but can be made simpler
    @inbounds v = p2v(p) / n0 # dividimos aquí para hacerlo más eficiente.
    for l in (-order):order 
      for m in (-order):order
      @inbounds S[mod1(j[1] + l, J[1]), mod1(j[2] + m, J[2])] += Shape(order, -y[1] + l) * Shape(order, -y[2] + m) * v;
      end
    end
  end
  return S[:,:] # allready normalized with n0
end


function get_current_rel_threads!(u, S, p)
  par_grid, TS = p
  N, L, J, dx, order = par_grid
  j = fill(Int64(1),nthreads()) 
  y = fill(Float64(1.0),nthreads())
  TS .= zeros(Float64)
  n0 = N/L
  @threads for i in 1:N
    @inbounds j[threadid()], y[threadid()] = get_index_and_y(u[i], J, L)
    @inbounds v = p2v(u[N+i]) / dx / n0 # the dx here is from the different definition from the paper
    for l in (-order):-j[threadid()]
      @inbounds TS[J + j[threadid()] + l, threadid()] += Shape(order, -y[threadid()] + l) * v
    end
    for l in max(-order,-j[threadid()]+1):min(order,J-j[threadid()])
      @inbounds TS[j[threadid()] + l, threadid()] += Shape(order, -y[threadid()] + l) * v
    end
    for l in J-j[threadid()]+1:order
      @inbounds TS[j[threadid()] - J + l, threadid()] += Shape(order, -y[threadid()] + l) * v
    end
    #= 
    for l in (-order):order
       @inbounds TS[mod1(j + l, J), threadid()] += Shape(order, -y + l) * v
    end
    =#
  end
  S .= zeros(Float64)
  @threads for j in 1:J
    for t in 1:nthreads()
      @inbounds S[j] += TS[j, t]
    end
  end
  S[:]
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

"""
Calculates the RHS of the evolution equation. 
Returns du 
uses several functions which are passed as parameters
p = N, J, L, dx, Density!, Electric!, Poisson1D!
"""


function RHSC_rel(u,t,p_RHSC)
  if nthreads() == 1
    N, J, L, dx, order, n, S, du, get_density!, get_current_rel!, Interpolate = p_RHSC
    par_grid = (N, L, J, dx, order)
    get_current_rel!(u, S, par_grid)
  else
    N, J, L, dx, order, n, S, du, get_density!, get_current_rel_threads!, Interpolate, TS = p_RHSC
    par_grid = (N, L, J, dx, order)
    get_current_rel_threads!(u, S, (par_grid, TS))
  end

    E = view(u,2N+1:2N+J)
    
    for i in 1:N        
      @inbounds du[i] = p2v(u[N+i]) # relativistic factor (u is the momentum)
      @inbounds du[N+i] = - Interpolate(order, E, u[i], J, L)
    end

    for j in 1:J
      @inbounds du[2N+j] =  S[j] # particles have negative sign!
    end
    return du[:]
end

function Coordinate_test(r,L)
    if minimum(r) < 0.0 
        error("negative coordinates")
    end
    if maximum(r) > L
        error("coordintates extend beyond L")
    end
end

function Coordinate_test(r,Box::Tuple)
  D = length(Box)÷2
  for d in 1:D
    if minimum(r[d:2D:2D*N]) < Box[2d-1] || maximum(r[d:2D:2D*N]) > Box[2d]
      error("particle out of Box")
    end
  end
end


function make_periodic!(r,L)
  return mod1.(r,L)
end

function make_periodic!(r,Box::Tuple)
  D = length(Box)÷2
  N = length(r)÷2÷D
  Box_array = [i for i in Box]
  L = Box_array[2:2:end] - Box_array[1:2:end-1]
  for i in 1:N
    for j in 1:D
      r[(i-1)*2*D+j] = mod1(r[(i-1)*2*D+j] - Box_array[2j-1], L[j]) + Box_array[2j-1]
    end
  end
  return r[:]
end



function RK4_Step!(f,y0,t0,h,p)
    k1 = h*f(y0,t0,p)
    k2 = h*f(y0+0.5*k1, t0+0.5*h,p)
    k3 = h*f(y0+0.5*k2, t0+0.5*h,p)
    k4 = h*f(y0+k3, t0+h,p)
    y0 .= y0 + (k1 + 2k2 + 2k3 + k4)/6
end

function get_index_and_distance(s,dx,L)
    #if s < 0
    #    s = s + L
    #end
    #if s > L
    #    s = s - L 
    #end
    s = mod1(s,L)
    j = floor(Int64, s ÷ dx) + 1 #find the grid space where it is.
    if j > J || j < 1
      error("j = $j")
    end
    y = (s % dx)/dx #how far is there 
    return j, y
end

"""
given a number s in between x_j and x_{j+1} computes y = (s - x_j)/dx and return j and y.
get_index_and_y(0.4,2,1)
dx = 0.5, j = 1, y = 0.4/0.5 = 0.8
get_index_and_y(0.7,2,1) = j = 2, 0.2/0.5= 0.4
"""
function get_index_and_y(s,J,L)
  s = (s/L*J + J)%J
  j = floor(Int,s) + 1
  #j = convert(Int64,s) + 1
  y = (s%1)
  return j, y
end

"""
same as before, but for an arbitrary box and dimension.
Jt = (100, 200, 200)
Box = (0.0, 10, 0.0, 20, -20.0, 0.0)
ss = [8.74, 8.74, -11.26]
j=[1,1,1]
y=zeros(3)
get_index_and_y!(j,y,ss,Jt,Box)
([88, 88, 88], [0.4000000000000057, 0.39999999999997726, 0.4000000000000057])
"""
function get_index_and_y!(j::Array{Int64,1}, y::Array{Float64,1}, s, J::Tuple,Box::Tuple)
  D = length(s)
  for i in 1:D
    y[i] =  (s[i]/(Box[2i] - Box[2i-1])*J[i] + J[i])%J[i]
    j[i] = floor(Int,y[i]) + 1
    y[i] = (y[i]%1)
  end
  #return j[:], y[:]
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
  Box, N, J= p
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
  
    get_Fields!(E,B,u)
    energy_E = sum.(E'*E+B'*B)
  
  # return energy_K,  dx * energy_E / 2 * n0
  return energy_K / n0,  dx * energy_E / 2 # normalized version
end
"""
Derivatives of shape functions.
They have support for -(order+1)/2 =< y =< (order+1)/2
"""
@inline function W(order::Int,y::Float64)
  #y = norm(y)
  y = abs(y)
  if order == 0
    return  (y <= 1/2) ? 1 : 0
  elseif order ==1
    return  (y <= 1) ? 1 - y : 0
  elseif order == 2
    return (y <= 1/2) ? 3/4 - y^2  : (((y > 1/2) && (y <= 3/2)) ? (3 - 2*y)^2 / 8 : 0)
  elseif order == 3
    return (y <= 1) ? 2/3 - y^2 + y^3 / 2 : (((y > 1) && (y <= 2)) ? (2 - y)^3 / 6 : 0)
  elseif order == 4
    return (y <= 1/2) ? 115/192 - 5y^2/8 + y^4/4 : (((y > 1/2) && (y <= 3/2)) ? (55 + 20y -120y^2 + 80y^3 - 16y^4)/96 : (((y > 3/2) && (y < 5/2)) ? (5 - 2y)^4/384 : 0))
  elseif order == 5
    #return (y <= 1) ? 11/20 - y^2/2 + y^4/4 - y^5/12 : (((y > 1) && (y <= 2)) ? 17/40 + 5y/8 - 7y^2/4 + 5y^3/4 - 3y^4/8 + y^5/24 : (((y > 2) && (y < 3)) ? (3 - y)^5/120 : 0))
    return (y <= 1) ? 11/20 - y^2/2 + y^4/4 - y^5/12 : (((y <= 2)) ? 17/40 + 5y/8 - 7y^2/4 + 5y^3/4 - 3y^4/8 + y^5/24 : (((y < 3)) ? (3 - y)^5/120 : 0))
  
  else
    error("order = $order not yet implemented ")
  end
end


"""
Alternative definition
"""
function W_alt(order::Int,y::Float64)
  y = abs(y)
  if order == 0
    return  (y > 1/2) ? 0 : 1
  elseif order ==1
    return  (y > 1) ? 0 : 1 - y 
  elseif order == 2
    #return (y <= 1/2) ? 3/4 - y^2  : (((y > 1/2) && (y <= 3/2)) ? (3 - 2*y)^2 / 8 : 0)
    return (y > 3/2) ? 0 : (((y > 1/2) && (y <= 3/2)) ? (3 - 2*y)^2 / 8 : 3/4 - y^2)
  elseif order == 3
    #return (y <= 1) ? 2/3 - y^2 + y^3 / 2 : (((y > 1) && (y <= 2)) ? (2 - y)^3 / 6 : 0)
    return (y > 2) ? 0 : (((y > 1) && (y <= 2)) ? (2 - y)^3 / 6 : 2/3 - y^2 + y^3 / 2)
  elseif order == 4
    #return (y <= 1/2) ? 115/192 - 5y^2/8 + y^4/4 : (((y > 1/2) && (y <= 3/2)) ? (55 + 20y -120y^2 + 80y^3 - 16y^4)/96 : (((y > 3/2) && (y < 5/2)) ? (5 - 2y)^4/384 : 0))
    return (y > 5/2) ? 0 : (((y > 3/2) && (y < 5/2)) ? (5 - 2y)^4/384 : (((y > 1/2) && (y <= 3/2)) ? (55 + 20y -120y^2 + 80y^3 - 16y^4)/96 : 115/192 - 5y^2/8 + y^4/4))
  elseif order == 5
    #return (y <= 1) ? 11/20 - y^2/2 + y^4/4 - y^5/12 : (((y > 1) && (y <= 2)) ? 17/40 + 5y/8 - 7y^2/4 + 5y^3/4 - 3y^4/8 + y^5/24 : (((y > 2) && (y < 3)) ? (3 - y)^5/120 : 0))
    return (y > 3) ? 0 : (((y > 2) && (y < 3)) ? (3 - y)^5/120 : (((y > 1) && (y <= 2)) ? 17/40 + 5y/8 - 7y^2/4 + 5y^3/4 - 3y^4/8 + y^5/24 : 11/20 - y^2/2 + y^4/4 - y^5/12))
  else
    error("order = $order not yet implemented ")
  end
end

"""
Sharp functions are W functions of an order less. 
The dx in the Appendix is added when used so as not to carry dx all over. 
So this definition DIFFERS FROM THE PAPER BY A DX!
They have support on -order/2 =< y =< order/2
The orders goes from 1 to 6 (order zero is not defined)
While the orders of W goes from 0 to 5
"""
@inline function Shape(order::Int,y::Float64) 
  return W(order - 1,y)
end

"""
This are interpolation functions for getting the Electric field correct.
According the SHARP the second is better. Since it keeps momentum conservation.
Modified so as to use the smallest stencils.
"""
@inline function Interpolate_1(order, vector, x, J, L)
  stencil = order÷2
  j, y = get_index_and_y(x,J,L)
  vi = 0.0
    for l in (-stencil):(stencil +1)
      vi += vector[mod1(j+l,J)] * W(order, -y + l)
    end
  return vi
end

@inline function Interpolate_2(order, vector, x, J, L)
  stencil = (order+1)÷2
  j, y = get_index_and_y(x,J,L)
  vi = 0.0
    for l in (-stencil):(stencil) 
      vi += (vector[mod1(j+l,J)] + vector[mod1(j+l+1,J)]) * W(order, -y + 1/2 + l) / 2
    end
  return vi
end
"""
Terminado, pero sin probar.
"""
@inline function Interpolate_per(order, vector, x, J, L)
  j, y = get_index_and_y(x,J,L)
  vi = 0.0
    for l in (-order):-j 
      vi += vector[J+j+l] * W(order, -y + 1/2 + l)
    end
    for l in  max(-order,-j+1):min(order,J-j)
      vi += vector[j+l] * W(order, -y + 1/2 + l)
    end
    for l in J-j+1:order
      vi += vector[j+l] * W(order, -y + 1/2 + l)
    end
    # now we increase j by 1
    j += 1
    for l in (-order):-j 
      vi += vector[J+j+l] * W(order, -y + 1/2 + l)
    end
    for l in  max(-order,-j+1):min(order,J-j)
      vi += vector[j+l] * W(order, -y + 1/2 + l)
    end
    for l in J-j+1:order
      vi += vector[j+l] * W(order, -y + 1/2 + l)
    end
  return vi / 2
end

"""
Structure and functions to work with particles 
This is probably much faster for the position and
velocity are adjacent in memory.

Initially the vector u in 1D was of the form: [x[N],v[N],E[J]]
in more dimensions it would be [x1[N],x2[N],x3[N],v1[N],v2[N],v3[N],E1[Jx*Jy],E2,E3,B1,B2,B3] or [x1,x2,v1,v2,E1,E2,B3]
So the sizes are: 1D 2N+J, 2D 4N+Jx*Jy*3, 3D 6N+Jx*Jy*Jz*6
uor (ordered u) has the form [x1[1],x2[1],x3[1],v1[1],v2[1],v3[1], etc.]
"""

mutable struct Particles
  r ::Float64
  v ::Float64
end

#pars = Vector{Particles}(undef, N)
#par1 = par(0.0,0.0)

#pp = fill(par1,N)
"""
Makes particles out of a vector with positions first and then
velocities
"""

function make_particles!(par_dis, pars)
  N = length(par_dis)÷2
  for i in 1:N
      pars[i] = Particles(par_dis[i], par_dis[N+i])
  end
end

"""
reorder particle vector so that position and velocity are contiguous
"""
function reorder_particles!(u,uro,dimension)
  N = length(u)÷2÷dimension # particles number
  for i in 1:N
    for d in dimension
      uro[2i - 1 + d] = u[dimension*i+d]
      uro[2i*d] = u[N+i]
    end
  end
end

"""
This function assumes that the particle vector is ordered one particle phase space after the other.
    D = 2
    x = zeros(D)
    par_dis = [1, 2, 3, 4, 5, 6, 7, 8]
    i = 2
    get_positions!(x,i,par_dis)
    > [5.0, 6.0]
"""
function get_positions!(x,i,par_dis)
  D = length(x)
  x[:] = par_dis[(i-1)*2*D+1:(i-1)*2*D+D]
end
"""
This function assumes that the particle vector is ordered one particle phase space after the other.
  D = 2
  p = zeros(D)
  par_dis = [1, 2, 3, 4, 5, 6, 7, 8]
  i = 2
  get_positions!(p,i,par_dis)
  > [7.0, 8.0]
"""
function get_momenta!(p,i,par_dis)
  D = length(p)
  p[:] = par_dis[(i-1)*2*D+1+D:i*2*D]
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


############## GRID FUNCTIONS #################
"""
Grid volume
"""
function volume(Box) 
  if length(Box) == 2
    return abs(Box[2]-Box[1])
  elseif length(Box) == 4
    return abs((Box[2]-Box[1])*(Box[4]-Box[3]))
  elseif length(Box) == 6
    return abs((Box[2]-Box[1])*(Box[4]-Box[3])*(Box[6]-Box[5]))
  else 
    error("wrong length for Box")
  end
end

"""
Grid Differential 
"""
function differentials(Box,J,periodic=true)
  D = length(J)
  dd = zeros(D)
  if periodic
    for i in 1:D
      dd[i] = (Box[2i] - Box[2i-1])/J[i]
    end
  else
    for i in 1:D
      dd[i] = (Box[2i] - Box[2i-1])/(J[i]-1)
    end
  end
  return dd[:]
end
