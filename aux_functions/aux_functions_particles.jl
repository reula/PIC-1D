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
  x[:] = par_dis[range_x(i, D)]
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
  p[:] = par_dis[range_p(i,D)]
end

@inline range_x(i, D) =  (i-1)*2*D+1:(i-1)*2*D+D

@inline range_p(i, D) = (i-1)*2*D+1+D:i*2*D


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
      error("particle out of Box min=$(minimum(r[d:2D:2D*N])), max=$(maximum(r[d:2D:2D*N]))")
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

function get_all_positions(u,D,N)
  x = zeros(D)
  pos = Array{Float64,1}(undef,2D*N)
  for i in N
    get_positions!(x,i,u[1:2D*N])
    for l in 1:D
      pos[i+(l-1)] = x[l]
    end
  end
  return pos[:]
end
    
