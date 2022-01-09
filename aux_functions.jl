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


"""The routine below evaluates the electron number density on an evenly spaced mesh given the instantaneous electron coordinates.

// Evaluates electron number density n(0:J-1) from 
// array r(0:N-1) of electron coordinates.
"""
function get_density!(r, n, dx)
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

function get_density!(u, n, p)
  L, N, J, κ, dx = p
  r = view(u,1:N)
  fill!(n,0.0)
  # Evaluate number density.
  for i in 1:N
    j, y = get_index_and_distance(r[i],dx,L)
    #if j < 1
    #    error("j = $j")
    #end 
    n[j] += (1. - y) / dx;
    if (j == J) 
        n[1] += y / dx
    else 
        n[j+1] += y / dx
    end
  end
  #n .= n/n0 - 1.0 # return rho directly
end

"""The routine below evaluates the electron current on an evenly spaced mesh given the instantaneous electron coordinates.

// Evaluates electron number density S(0:J-1) from 
// array r(0:N-1) of electron coordinates.
"""
function get_current!(u, S, p)
  L, N, J, κ, dx = p
  r = view(u,1:N)
  v = view(u,N+1:2N)
  fill!(S,0.0)
  # Evaluate number density.
  for i in 1:N
    j, y = get_index_and_distance(r[i],dx,L)
    #if j < 1 || j > J
    #    error("j = $j")
    #end 
    S[j] += (1. - y)*v[i] / dx;
    if (j == J) 
        S[1] += y*v[i] / dx
    else 
        S[j+1] += y*v[i] / dx
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
        j, y = get_index_and_distance(u[i],dx,L)
        
        if (j == J)
            Efield = E[j] * (1. - y) + E[1] * y;
        else
            Efield = E[j] * (1. - y) + E[j+1] * y;
        end

        du[i] = u[N+i]
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
function RHSC(u,t,p_RHSC)
  N, J, L, dx, n, S, du, get_density!, get_current! = p_RHSC
    p = L, N, J, κ, dx
    #get_density!(u, n, p)
    get_current!(u, S, p)
    E = view(u,2N+1:2N+J)
    n0 = N/L

    for i in 1:N
        j, y = get_index_and_distance(u[i],dx,L)
        
        if (j == J)
            Efield = E[j] * (1. - y) + E[1] * y;
        else
            Efield = E[j] * (1. - y) + E[j+1] * y;
        end

        du[i] = u[N+i]
        du[N+i] = - Efield
    end

    if false
      du[2N+1] =  (S[J] + S[1])/n0/2
      for j in 2:J
        du[2N+j] =  (S[j-1] + S[j])/n0/2 # particles have negative sign!
      end
    end
    if true
      for j in 1:J
        du[2N+j] =  S[j]/n0 # particles have negative sign!
      end
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

function make_periodic!(r,L)
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


function RK4_Step!(f,y0,t0,h,p)
    k1 = h*f(y0,t0,p)
    k2 = h*f(y0+0.5*k1, t0+0.5*h,p)
    k3 = h*f(y0+0.5*k2, t0+0.5*h,p)
    k4 = h*f(y0+k3, t0+h,p)
    y0 .= y0 + (k1 + 2k2 + 2k3 + k4)/6
end

function get_index_and_distance(s,dx,L)
    if s < 0
        s = s + L
    end
    if s > L
        s = s - L 
    end
    j = convert(Int64, s ÷ dx) + 1 #find the grid space where it is.
    if j > J || j < 1
      error("j = $j")
    end
    y = (s % dx)/dx #how far is there 
    return j, y
end
