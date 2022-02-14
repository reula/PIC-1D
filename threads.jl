#Código de Charly

using FileIO
using Base.Threads
using BenchmarkTools

const order = 1
const L = 5
const N = 80000
#const N = 20000
const J = 50
const κ = 2π/L # for Fourier Transform
const dx = L/J

include("aux_functions.jl")

println("threads = $(nthreads())")

function get_current_threads_soa(r, v, p)
  L, N, J, κ, dx, order = p

  TS = zeros(J, nthreads())
  @threads for i in 1:N
    @inbounds j, y = get_index_and_y(r[i], J, L)
    for l in (-order):-j
      @inbounds TS[J + j + l, threadid()] += W(order, -y + l) * v[i] / dx;
    end
    for l in max(-order,-j+1):min(order,J-j)
      @inbounds TS[j + l, threadid()] += W(order, -y + l) * v[i] / dx;
    end
    for l in J-j+1:order
      @inbounds TS[j - J + l, threadid()] += W(order, -y + l) * v[i] / dx;
    end
    # for l in (-order):order
    #   @inbounds TS[mod1(j + l, J), threadid()] += W(order, -y + l) * v[i] / dx;
    # end
  end

  S = zeros(J)
  @threads for i in 1:J
    for t in 1:nthreads()
      @inbounds S[i] += TS[i, t]
    end
  end
  S
end

function get_current_threads(u, p)
  L, N, J, κ, dx, order = p

  TS = zeros(J, nthreads())
  @threads for i in 1:N
    @inbounds j, y = get_index_and_y(u[i], J, L)
    for l in (-order):-j
      @inbounds TS[J + j + l, threadid()] += W(order, -y + l) * u[N+i] / dx;
    end
    for l in max(-order,-j+1):min(order,J-j)
      @inbounds TS[j + l, threadid()] += W(order, -y + l) * u[N+i] / dx;
    end
    for l in J-j+1:order
      @inbounds TS[j - J + l, threadid()] += W(order, -y + l) * u[N+i] / dx;
    end
    # for l in (-order):order
    #   @inbounds TS[mod1(j + l, J), threadid()] += W(order, -y + l) * v[i] / dx;
    # end
  end

  S = zeros(J)
  @threads for i in 1:J
    for t in 1:nthreads()
      @inbounds S[i] += TS[i, t]
    end
  end
  S
end

function get_current_threads_aos(rv, p)
  L, N, J, κ, dx, order = p

  TS = zeros(J, nthreads())
  @threads for i in 1:N
    @inbounds j, y = get_index_and_y(rv[1,i], J, L)
    for l in (-order):-j
      @inbounds TS[J + j + l, threadid()] += W(order, -y + l) * rv[2,i] / dx;
    end
    for l in max(-order,-j+1):min(order,J-j)
      @inbounds TS[j + l, threadid()] += W(order, -y + l) * rv[2,i] / dx;
    end
    for l in J-j+1:order
      @inbounds TS[j - J + l, threadid()] += W(order, -y + l) * rv[2,i] / dx;
    end
    # for l in (-order):order
    #   @inbounds S[mod1(j + l, J), threadid()] += W(order, -y + l) * rv[2,i] / dx;
    # end
  end

  S = zeros(J)
  @threads for i in 1:J
    for t in 1:nthreads()
      @inbounds S[i] += TS[i, t]
    end
  end
  S
end


"""
reorder particle vector so that position and velocity are contiguous
"""
function split_particles_soa!(u, r, v)
  r[:] = u[1:N]
  v[:] = u[N+1:2N]
end

function split_particles_aos!(u, rv)
  for i in 1:N
    rv[1, i] = u[i]
    rv[2, i] = u[N+i]
  end
end

function load_soa()
  pars = zeros(2N)
  #pars = dones((N,), workers(), [nworkers()] );
  par_dis = load("Initial_Distributions/par_dis_L5_N8_5_theta001.jld2", "par_dis");

  r = zeros(N)
  v = zeros(N)
  split_particles_soa!(par_dis, r, v)
  (r, v)
end

function load_aos()
  pars = zeros(2N)
  #pars = dones((N,), workers(), [nworkers()] );
  par_dis = load("Initial_Distributions/par_dis_L5_N8_5_theta001.jld2", "par_dis");

  rv = zeros(2,N)
  split_particles_aos!(par_dis, rv)
  rv
end

x = [dx*(i-1) for i in 1:J]
p = (L, N, J, κ, dx, order)

#rv = load_aos()
#@benchmark get_current_threads_aos(rv, p)

(r,v) = load_soa()
#@benchmark get_current_threads_soa(r, v, p)

u = [r;v]
#@benchmark get_current_threads(u, p)



