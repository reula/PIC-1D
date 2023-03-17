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
      N, J, Box, order, n, S, du, get_density!, get_current, Interpolate,  Dx, Δx, σx, Dy, Δy, σy = p_RHSC
      par_grid = (N, J, Box, order)
      get_current(u, S, par_grid)
    else
      N, J, Box, order, n, S, du, get_density!, get_current_threads, Interpolate, TS,  Dx, Δx, σx, Dy, Δy, σy  = p_RHSC
      par_grid = (N, J, Box, order)
      get_current_threads(u, S, (par_grid, TS))
    end
  
      #E = view(u,2N+1:2N+J)
      Fu = view(u,4N+1:4N+3*prod(J))
      F = reshape(Fu,(3,J...))
      E = F[1:2,:,:]
      B = F[3,:,:]
      
      dFu = view(du,4N+1:4N+3*prod(J))
      dF = reshape(dFu,(3,J...))

      @threads for i in 1:J[1]
        mul!(view(dF,1,i,:), Dy, view(F,3,i,:),one(eltype(F)))
        mul!(view(dF,3,i,:), Dy , view(F,1,i,:),one(eltype(F)))
        end
        @threads for j in 1:J[2]
        mul!(view(dF,2,:,j), Dx, view(F,3,:,j),-one(eltype(F)))
        mul!(view(dF,3,:,j), Dx, view(F,2,:,j),-one(eltype(F)),one(eltype(F)))
        end

      @threads for i in 1:J[1]
        mul!(view(dF,1,i,:), Δy, view(F,1,i,:), σy, one(eltype(F)))
        mul!(view(dF,2,i,:), Δy, view(F,2,i,:), σy, one(eltype(F)))
        mul!(view(dF,3,i,:), Δy, view(F,3,i,:), σy, one(eltype(F)))
        end
        @threads for j in 1:J[2]
        mul!(view(dF,1,:,j), Δx, view(F,1,:,j), σx, one(eltype(F)))
        mul!(view(dF,2,:,j), Δx, view(F,2,:,j), σx, one(eltype(F)))
        mul!(view(dF,3,:,j), Δx, view(F,3,:,j), σx, one(eltype(F)))
        end

      @threads for j in 1:J[2]
        for i in 1:J[1]
            for l in 1:2
         dF[l,i,j] +=  S[l,i,j] # particles have negative sign!
            end
        end
      end

      for i in 1:N        
        v = p2v(u[range_p(i, D)])
         du[range_x(i, D)] = v # relativistic factor (u is the momentum)
         du[range_p(i, D)] = - Interpolate(order, E, B, v, u[range_x(i, D)], J, Box)
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



