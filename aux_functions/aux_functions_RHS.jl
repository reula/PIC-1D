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

function RHS_D(u,t,p_RHSC)
    if nthreads() == 1
      N, J, Box, order, n, S, du, get_density!, get_current_rel_2D!, Interpolate, Dx, Dy = p_RHSC
      par_grid = (N, J, Box, order)
      get_current_rel_2D!(u, S, par_grid)
    else
      N, J, Box, order, n, S, du, get_density!, get_current_threads_2D!, Interpolate, TS, Dx, Dy  = p_RHSC
      par_grid = (N, J, Box, order)
      get_current_threads_2D!(u, S, (par_grid, TS))
    end
  
      #E = view(u,2N+1:2N+J)
      F = reshape(u[N+1:end],(3,J...))
      E = F[1:2,:,:]
      B = F[3,:,:]
      
      for i in 1:N        
        v = p2v(u[range_p(i, D)])
        @inbounds du[range_x(i, D)] = v # relativistic factor (u is the momentum)
        @inbounds du[range_p(i, D)] = - Interpolate(order, E - [v[2],-v[1]]*B, u[range_x(i, D)], J, Box)
      end
      
      dF = reshape(du[N+1:end],(3,J...))

      @threads for i in 1:J[1]
        mul!(view(dF,1,i,:),Dy, view(F,3,i,:),one(eltype(F)))
        mul!(view(dF,3,i,:),Dy, view(F,1,i,:),one(eltype(F)))
        end
        @threads for j in 1:J[2]
        mul!(view(dF,2,:,j),Dx, view(F,3,:,j),-one(eltype(F)))
        mul!(view(dF,3,:,j),Dx, view(F,2,:,j),-one(eltype(F)),one(eltype(F)))
        end

      @threads for j in 1:J[2]
        for i in 1:J[1]
            for l in 1:2
        @inbounds du[l,i,j] +=  S[l,i,j] # particles have negative sign!
            end
        end
      end
      return du[:]
  end



function RK4_Step!(f,y0,t0,h,p)
    k1 = h*f(y0,t0,p)
    k2 = h*f(y0+0.5*k1, t0+0.5*h,p)
    k3 = h*f(y0+0.5*k2, t0+0.5*h,p)
    k4 = h*f(y0+k3, t0+h,p)
    y0 .= y0 + (k1 + 2k2 + 2k3 + k4)/6
end


function RK4_Step!!(f,y0,t0,h,p)
    k1 = h*f(y0,t0,p)
    k2 = h*f(y0+0.5*k1, t0+0.5*h,p)
    k3 = h*f(y0+0.5*k2, t0+0.5*h,p)
    k4 = h*f(y0+k3, t0+h,p)
    y0 .= y0 + (k1 + 2k2 + 2k3 + k4)/6
end
