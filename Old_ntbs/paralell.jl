
using Transducers
using FileIO

nth = Threads.nthreads()
println("nth = $(nth)")

order = 1
const L = 5
#N = 80000
const N = 20000
const J = 50
const κ = 2π/L # for Fourier Transform
dx = L/J
x = [dx*(i-1) for i in 1:J] ;
p = (L, N, J, κ, dx, order)

function get_current_par!(u, S, p, th)
    L, N, J, κ, dx, order, nth = p
    #th = threadid()
    #nth = nthreads()
    NT = N ÷ nth
    r = view(u,(th-1)*NT+1:NT*th)
    v = view(u,N+(th-1)*NT+1:N+NT*th)
    fill!(S[:,th],0.0)
  
    for i in 1:N
      j, y = get_index_and_y(r[i],J,L)
      for l in (-order):order 
        S[mod1(j + l, J),th] += W(order, -y + l) * v[i] / dx;
      end
    end
    return S[:]
  end