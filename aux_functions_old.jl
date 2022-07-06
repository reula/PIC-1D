using Distributed
using DistributedArrays
using DistributedArrays.SPMD
#@everywhere using SharedArrays
@everywhere using Distributed
@everywhere using DistributedArrays
@everywhere using DistributedArrays.SPMD

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

"""
Takes out the mass of the grid function so that the sum is now null
"""
function filter_constant!(E)
  J = length(E)
  V = rfft(E)
  V[1] = 0.0 #extract the first component
  E[:] = irfft(V,J)
end

"""The routine below evaluates the electron number density on an evenly spaced mesh given the instantaneous electron coordinates.

// Evaluates electron number density n(0:J-1) from 
// array r(0:N-1) of electron coordinates.
"""
function get_density_old_old!(r, n, dx)
  # Initialize 
  N = length(r)
  J = length(n)
  fill!(n,0.0)
  # Evaluate number density.
  for i in 1:N
    j, y = get_index_and_distance(r[i],dx,L)
    if j < 1
        error("j = $j")
    end 
    n[j] += (1. - y) / dx;
    if (j == J) 
        n[1] += y / dx
    else 
        n[j+1] += y / dx
    end
  end
end

function get_density_old!(u, n, par_grid)
  N, L, J, dx, order = par_grid
  r = view(u,1:N)
  fill!(n,0.0)
  # Evaluate number density.
  for i in 1:N
    j, y = get_index_and_y(r[i],J,L)
    #if j < 1
    #    error("j = $j")
    #end 
    if order == 1
      n[j] += (1. - y) / dx;
      if (j == J) 
          n[1] += y / dx
      else 
        n[j+1] += y / dx
      end
    end
    if order == 2
      if y <= 1/2
        n[j] += (3/4 - y^2) / dx
        if j == 1
          n[J] += (1-2y)^2/8 / dx
          n[2] += (1+2y)^2/8 / dx
        elseif j==J
          n[1] += (1+2y)^2/8 / dx
          n[J-1] += (1-2y)^2/8 / dx
        else
          n[j+1] += (1+2y)^2/8 / dx
          n[j-1] += (1-2y)^2/8 / dx
        end
      elseif y <= 3/2
        n[j] += (3 - 2y)^2/8 / dx
        if j == J-1
          n[J] += (3/4 - (1-y)^2) / dx
          n[1] += (1-2y)^2/8 / dx
        elseif j==J
          n[1] += (3/4 - (1-y)^2) / dx
          n[2] += (1-2y)^2/8 / dx
        else
          n[j+1] += (3/4 - (1-y)^2) / dx
          n[j+2] += (1-2y)^2/8 / dx
        end
      end
    end
  end
    #n .= n/n0 - 1.0 # return rho directly
end

function get_density!(u, n, par_grid)
  N, L, J, dx, order = par_grid
  n0 = N/J
  r = view(u,1:N)
  fill!(n,0.0)
  # Evaluate number density.
  for i in 1:N
    @inbounds j, y = get_index_and_y(r[i],J,L)
    for l in (-order):order 
      @inbounds n[mod1(j + l, J)] += W(order, -y + l) / dx;
    end
  end
    n .= n/n0 # return rho directly (we need to subtract 1 in cases where we assume positive particles, but this is done elsewhere.)
end

function get_total_charge(ρ,par)
  J, dx = par
  Q = 0.0
  for i in 1:J
      Q += ρ[i]
  end
  return Q * dx
end

function get_density_ro!(u, n, par_grid)
  N, L, J, dx, order = par_grid
  r = view(u,1:2:2N-1)
  fill!(n,0.0)
  n0 = N/J
  # Evaluate number density.
  for i in 1:N
    j, y = get_index_and_y(r[i],J,L)
    for l in (-order):order 
      n[mod1(j + l, J)] += W(order, -y + l) / dx;
    end
  end
    n .= n/n0 #- 1.0 # return rho directly
end

"""The routine below evaluates the electron current on an evenly spaced mesh given the instantaneous electron coordinates.

// Evaluates electron number density S(0:J-1) from 
// array r(0:N-1) of electron coordinates.
"""
function get_current_old!(u, S, par_grid)
  N, L, J, dx, order = par_grid
  r = view(u,1:N)
  v = view(u,N+1:2N)
  fill!(S,0.0)
  # Evaluate number density.
  for i in 1:N
    #j, y = get_index_and_distance(r[i],dx,L)
    j, y = get_index_and_y(r[i],J,L)
    #if j < 1 || j > J
    #    error("j = $j")
    #end 
    if order == 1
      S[j] += (1. - y)*v[i] / dx;
      if (j == J) 
        S[1] += y*v[i] / dx
      else 
        S[j+1] += y*v[i] / dx
      end
    end
    if order == 2
      if y <= 1/2
        S[j] += (3/4 - y^2) *v[i] / dx
        if j == 1
          S[J] += (1-2y)^2/8 *v[i] / dx
          S[2] += (1+2y)^2/8 *v[i] / dx
        elseif j==J
          S[1] += (1+2y)^2/8 *v[i] / dx
          S[J-1] += (1-2y)^2/8 *v[i] / dx
        else
          S[j+1] += (1+2y)^2/8 *v[i] / dx
          S[j-1] += (1-2y)^2/8 *v[i] / dx
        end
      elseif y <= 3/2
        S[j] += (3 - 2y)^2/8 *v[i] / dx
        if j == J-1
          S[J] += (3/4 - (1-y)^2) *v[i] / dx
          S[1] += (1-2y)^2/8 *v[i] / dx
        elseif j==J
          S[1] += (3/4 - (1-y)^2) *v[i] / dx
          S[2] += (1-2y)^2/8 *v[i] / dx
        else
          S[j+1] += (3/4 - (1-y)^2) *v[i] / dx
          S[j+2] += (1-2y)^2/8 *v[i] / dx
        end
      end
    end
  end
end

function get_current!(u, S, par_grid)
  N, L, J, dx, order = par_grid
  r = view(u,1:N)
  v = view(u,N+1:2N)
  fill!(S,0.0)

  for i in 1:N
    @inbounds j, y = get_index_and_y(r[i],J,L)
    for l in (-order):order 
      @inbounds S[mod1(j + l, J)] += W(order, -y + l) * v[i]
      #S[mod1(j + l, J)] += W_alt(order, -y + l) * v[i]
    end
  end
  S / dx
end

function get_current_rel!(u, S, par_grid)
  N, L, J, dx, order = par_grid
  r = view(u,1:N)
  p = view(u,N+1:2N) # in the relativistic version we compute p instead of v
  fill!(S,0.0)
  n0 = N/J
  for i in 1:N
    @inbounds j, y = get_index_and_y(r[i],J,L)
    @inbounds v = p[i]/sqrt(1+p[i]^2)/dx/n0
    for l in (-order):order 
      @inbounds S[mod1(j + l, J)] += W(order, -y + l) * v;
      #S[mod1(j + l, J)] += W_alt(order, -y + l) * v / dx;
    end
  end
  return S[:] # allready normalized with n0
end

function get_current_ro!(u, S, par_grid)
  N, L, J, dx, order = par_grid
  r = view(u,1:2:2N-1)
  v = view(u,2:2:2N)
  fill!(S,0.0)
  n0 = N/J
  for i in 1:N
    @inbounds j, y = get_index_and_y(r[i],J,L)
    for l in (-order):order 
      @inbounds S[mod1(j + l, J)] += W(order, -y + l) * v[i] / dx/ n0;
    end
  end
end

function get_current(u, S, par_grid)
  N, L, J, dx, order = par_grid
  r = view(u,1:N)
  v = view(u,N+1:2N)
  fill!(S,0.0)
  n0 = N/J

  for i in 1:N
    @inbounds j, y = get_index_and_y(r[i],J,L)
    for l in (-order):order 
      @inbounds  S[mod1(j + l, J)] += W(order, -y + l) * v[i] / dx / n0;
    end
  end
  return S[:]
end

function get_mini_current(pp, S, par_grid)
  N, L, J, dx, order = par_grid
  fill!(S,0.0)
  n0 = N/J
    j, y = get_index_and_y(pp.r,J,L)
    for l in (-order):order 
      S[mod1(j + l, J)] = W(order, -y + l) * pp.v / dx / n0;
    end
  return S[:]
end


function get_current_threads!(u, S, par_grid)
  N, L, J, dx, order = par_grid
  j = fill(Int64(1),nthreads()) 
  y = fill(Float64(1.0),nthreads())
  TS .= zeros(Float64)
  n0 = N/J
  @threads for i in 1:N
    @inbounds j[threadid()], y[threadid()] = get_index_and_y(u[i], J, L)
    #@inbounds p = u[N+i] # in the relativistic version we carry p, so we need to transform to velocities
    #@inbounds p = p/sqrt(1+p^2) / dx # we divide also by dx so as to save computations
    for l in (-order):-j[threadid()]
      @inbounds TS[J + j[threadid()] + l, threadid()] += W(order, -y[threadid()] + l) * u[N+i]/sqrt(1+u[N+i]^2);
    end
    for l in max(-order,-j[threadid()]+1):min(order,J-j[threadid()])
      @inbounds TS[j[threadid()] + l, threadid()] += W(order, -y[threadid()] + l) * u[N+i]/sqrt(1+u[N+i]^2);
    end
    for l in J-j[threadid()]+1:order
      @inbounds TS[j[threadid()] - J + l, threadid()] += W(order, -y[threadid()] + l) * u[N+i]/sqrt(1+u[N+i]^2);
    end
    # for l in (-order):order
    #   @inbounds TS[mod1(j + l, J), threadid()] += W(order, -y + l) * v[i] / dx;
    # end
  end

  S .= zeros(Float64)
  @threads for i in 1:J
    for t in 1:nthreads()
      @inbounds S[i] += TS[i, t]/dx/n0
    end
  end
  return S[:]
end

function get_current_rel_threads!(u, S, p)
  par_grid, TS = p
  N, L, J, dx, order = par_grid
  j = fill(Int64(1),nthreads()) 
  y = fill(Float64(1.0),nthreads())
  TS .= zeros(Float64)
  n0 = N/J
  @threads for i in 1:N
    @inbounds j[threadid()], y[threadid()] = get_index_and_y(u[i], J, L)
    for l in (-order):-j[threadid()]
      @inbounds TS[J + j[threadid()] + l, threadid()] += W(order, -y[threadid()] + l) * u[N+i]/sqrt(1+u[N+i]^2)
    end
    for l in max(-order,-j[threadid()]+1):min(order,J-j[threadid()])
      @inbounds TS[j[threadid()] + l, threadid()] += W(order, -y[threadid()] + l) * u[N+i]/sqrt(1+u[N+i]^2)
    end
    for l in J-j[threadid()]+1:order
      @inbounds TS[j[threadid()] - J + l, threadid()] += W(order, -y[threadid()] + l) * u[N+i]/sqrt(1+u[N+i]^2)
    end
    #= 
    for l in (-order):order
       @inbounds TS[mod1(j + l, J), threadid()] += W(order, -y + l) * u[N+i]/sqrt(1+u[N+i]^2)
    end
    =#
  end
  S .= zeros(Float64)
  @threads for i in 1:J
    for t in 1:nthreads()
      @inbounds S[i] += TS[i, t]/dx/n0
    end
  end
  S[:]
end


function get_current_ro_par(Dpars::DArray, DS::DArray, par_grid)
  N, L, J, dx, order = par_grid
  ran_pars = first(localindices(Dpars))
  #assuming ran_pars[1] is odd
  odd_range = first(ran_pars):2:last(ran_pars)
  LS = DS[:L]
  n0 = N/J
  fill!(LS,0.0)
  for i in odd_range
    j, y = get_index_and_y(Dpars[i],J,L)
    for l in (-order):order 
      LS[mod1(j + l, J),:] .+= W(order, -y + l) * Dpars[i+1] / dx / n0;
    end
  end
end

"""
Calculates the RHS of the evolution equation. 
Returns du 
uses several functions which are passed as parameters
p = N, J, L, dx, Density!, Electric!, Poisson1D!
"""
function RHS(u,t,p_RHC)
    N, J, L, dx, n, E, ϕ, du, Density!, Electric!, Poisson1D! = p_RHC
    p = L, N, J, κ, dx
    get_density!(u[1:N], n, p)
    n0 = N/L
    get_ϕ!(ϕ, n/n0 .-1.0,κ)
    get_E_from_ϕ!(ϕ,E,dx)

    for i in 1:N
        #j, y = get_index_and_distance(u[i],dx,L)
        j, y = get_index_and_y(u[i],J,L)
        if (j == J)
            Efield = E[j] * (1. - y) + E[1] * y;
        else
            Efield = E[j] * (1. - y) + E[j+1] * y;
        end

        du[i] = u[N+i]/sqrt(1 + u[N+i]^2) # relativistic factor
        du[N+i] = - Efield
    end
    return du[:]
end


"""
Calculates the RHS of the evolution equation. 
Returns du 
uses several functions which are passed as parameters
p = N, J, L, dx, n, S, du, get_current! 
"""
function RHSC_old(u,t,p_RHSC)
  N, J, L, dx, order, n, S, du, get_density!, get_current! = p_RHSC
    p = L, N, J, κ, dx, order
    #get_density!(u, n, p)
    get_current!(u, S, p)
    E = view(u,2N+1:2N+J)
    n0 = N/L

    for i in 1:N
        #j, y = get_index_and_distance(u[i],dx,L)
        j, y = get_index_and_y(u[i],J,L)
        
        if order == 1
          if (j == J)
              Efield = E[j] * (1. - y) + E[1] * y;
          else
              Efield = E[j] * (1. - y) + E[j+1] * y;
          end
        elseif order == 2
          if y <= 1/2
            if j == 1
              Efield = E[j] * (3/4 - y^2) + E[2] * (1+2y)^2/8 + E[J] * (1-2y)^2/8
            elseif j==J
              Efield = E[j] * (3/4 - y^2) + E[1] * (1+2y)^2/8 + E[J-1] * (1-2y)^2/8
            else
              Efield = E[j] * (3/4 - y^2) + E[j+1] * (1+2y)^2/8 + E[j-1] * (1-2y)^2/8
            end
          elseif y <= 3/2
            if j == J-1
              Efield = E[j] * (3 - 2y)^2/8 + E[J] * (3/4 - (1-y)^2) + E[1] * (1-2y)^2/8
            elseif j == J
              Efield = E[j] * (3 - 2y)^2/8 + E[1] * (3/4 - (1-y)^2) + E[2] * (1-2y)^2/8
            else
              Efield = E[j] * (3 - 2y)^2/8 + E[j+1] * (3/4 - (1-y)^2) + E[j+2] * (1-2y)^2/8
            end
          end
        end

        du[i] = u[N+i]
        du[N+i] = - Efield
    end

      for j in 1:J
        du[2N+j] =  S[j]/n0 # particles have negative sign!
      end

    return du[:]
end

function RHSC(u,t,p_RHSC)
  if nthreads() == 1
    N, J, L, dx, order, n, S, du, get_density!, get_current!, Interpolate = p_RHSC
    p = (N, L, J, κ, dx, order)
    get_current!(u, S, p)
  else
    N, J, L, dx, order, n, S, du, get_density!, get_current!, Interpolate, TS = p_RHSC
    p = (N, L, J, κ, dx, order, TS)
    get_current_threads!(u, S, p)
  end
    
    E = view(u,2N+1:2N+J)
    
    #@threads for i in 1:N
    for i in 1:N
        #j, y = get_index_and_y(u[i],J,L)
        #Efield = 0.0
        #for l in (-order):order 
        #  Efield += E[mod1(j+l,J)] * W(order, -y + l)
        #end

        @inbounds du[i] = u[N+i]/sqrt(1 + u[N+i]^2) # relativistic factor
        @inbounds du[N+i] = - Interpolate(order, E, u[i], J, L)
    end

    #@threads for j in 1:J
    for j in 1:J
        @inbounds du[2N+j] =  S[j] # particles have negative sign!
    end
    return du[:]
end

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
    #n0 = N/J
    
    #@threads for i in 1:N
    for i in 1:N
        #j, y = get_index_and_y(u[i],J,L)
        #Efield = 0.0
        #for l in (-order):order 
        #  Efield += E[mod1(j+l,J)] * W(order, -y + l)
        #end

        @inbounds du[i] = u[N+i]/sqrt(1 + u[N+i]^2) # relativistic factor (u is the momentum)
        @inbounds du[N+i] = - Interpolate(order, E, u[i], J, L)
    end

    #@threads for j in 1:J
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

function make_periodic_old!(r,L)
    N = length(r)
    for i in 1:N
      if (r[i] < 0.) 
        r[i] = r[i] + L;
      end
      if (r[i] > L) 
        r[i] = r[i] - L;
      end
    end
    return r[:]
end

function make_periodic!(r,L)
  return mod1.(r,L)
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

function get_energy(u,p)
  L, N, J = p
  dx = L/J
  n0 = N/L
  energy_K = 0.0
  energy_E = 0.0
  for i in 1:N
    energy_K = energy_K + u[N+i]^2
  end
  for j in 1:J
    energy_E = energy_E + u[2N+j]^2
  end
  
  return energy_K/2,  dx*energy_E /2 * n0
end

function get_energy_rel(u,p)
  L, N, J = p
  dx = L/J
  n0 = N/L
  energy_K = 0.0
  energy_E = 0.0
  for i in 1:N
    energy_K = energy_K + (sqrt(1+u[N+i]^2) - 1)
  end
  for j in 1:J
    energy_E = energy_E + u[2N+j]^2
  end
  
  return energy_K,  dx*energy_E /2 * n0
end

"""
Derivatives of shape functions
"""
function W(order::Int,y::Float64)
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
    return (y <= 1) ? 11/20 - y^2/2 + y^4/4 - y^5/12 : (((y > 1) && (y <= 2)) ? 17/40 + 5y/8 - 7y^2/4 + 5y^3/4 - 3y^4/8 + y^5/24 : (((y > 2) && (y < 3)) ? (3 - y)^5/120 : 0))
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
This are interpolation functions for getting the Electric field correct.
According the SHARP the second is better. Since it keeps momentum conservation.
"""
function Interpolate_1(order, vector, x, J, L)
  j, y = get_index_and_y(x,J,L)
  vi = 0.0
    for l in (-order+1):order 
      vi += vector[mod1(j+l,J)] * W(order, -y + l)
    end
  return vi
end

function Interpolate_2(order, vector, x, J, L)
  j, y = get_index_and_y(x,J,L)
  vi = 0.0
    for l in (-order):order 
      vi += (vector[mod1(j+l,J)] + vector[mod1(j+l+1,J)])/ 2 * W(order, -y + 1/2 + l)
    end
  return vi
end
"""
Terminado, pero sin probar.
"""
function Interpolate_per(order, vector, x, J, L)
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
function reorder_particles!(u,uro)
  N = length(u)÷2
  for i in 1:N
    uro[2i - 1] = u[i]
    uro[2i] = u[N+i]
  end
end

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
  #E_K = zeros(J)
  #P = zeros(J)
  ρ = zeros(J)
  S = zeros(J)

  for j in 1:M_g
      (Energy_K[j],Energy_E[j]) = get_energy_rel(v[:,j],(L,N,J))
      EField_T[j] = sum(v[2N+1:end,j])*dx
      p_T[j] = sum(v[N+1:2N])*dx
      get_density!(v[:,j], ρ, par_grid)
      get_current_rel!(v[:,j], S, par_grid)
      Q_T[j] = get_total_charge(ρ,(J, dx))
      S_T[j] = sum(S)/N/Q_T[j]
  end
  return Energy_K, Energy_E, EField_T, p_T, Q_T, S_T
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

