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
      N, J, Box, order, n, S, du, get_density!, get_current, Interpolate,  Dx, Δx, σx, Dy, Δy, σy, no_maxwell, dissipation = p_RHSC
      par_grid = (N, J, Box, order)
      get_current(u, S, par_grid)
    else
      N, J, Box, order, n, S, du, get_density!, get_current_threads, Interpolate,  Dx, Δx, σx, Dy, Δy, σy, no_maxwell, dissipation  = p_RHSC
      par_grid = (N, J, Box, order)
      S = get_current_threads(Val(order), Box_x, u)
    end
    #@show norm(S)
    make_periodic!(u,Box_x,N)
    Fu = view(u,4N+1:4N+3*prod(J))
    F = reshape(Fu,(3,J...))
    E = F[1:2,:,:]
    B = F[3,:,:]
      
    du .= 0.0
    dFu = view(du,4N+1:4N+3*prod(J))
    dF = reshape(dFu,(3,J...))

    if no_maxwell == false  #take away waves if false
      @threads for i in 1:J[1]
        mul!(view(dF,1,i,:), Dy, view(F,3,i,:),one(eltype(F)))
        mul!(view(dF,3,i,:), Dy , view(F,1,i,:),one(eltype(F)))
        end
      @threads for j in 1:J[2]
        mul!(view(dF,2,:,j), Dx, view(F,3,:,j),-one(eltype(F)))
        mul!(view(dF,3,:,j), Dx, view(F,2,:,j),-one(eltype(F)),one(eltype(F)))
      end
        if dissipation # take away dissipation
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
        end
    else
        du[4N+1:4N+3*prod(J)] .= 0.0
    end

      @threads for j in 1:J[2]
        for i in 1:J[1]
            for l in 1:2
         #dF[l,i,j] +=  S[l,i,j] # particles have negative sign!
         dF[l,i,j] +=  S[i,j,l]
            end
        end
      end

      @threads for i in 1:N        
        @views v = p2v(u[i*2D-D+1:i*2D])
        # v = p2v(u[range_p(i, D)])(i-1)*2*D+1+D:i*2*D
         du[range_x(i, D)] = v # relativistic factor (u is the momentum)
         du[i*2D-D+1:i*2D] = - Interpolate(order, E, B, v, u[range_x(i, D)], J, Box)
      end
      #@show norm(du[4N+1:4N+3*prod(J)])
    return du[:]
end

function RHS_D_slim!(u,t,p_RHSC) #version to optimize
  Order ,N, J, Box, _, n, S, du, get_density!, get_current, Interpolate,  Dx, Δx, σx, Dy, Δy, σy, no_maxwell, dissipation  = p_RHSC
  (order, N, J, Box_x, order, n, S, du, get_density_2D!, get_current_slim, Interpolate_All_EBv_2_slim, Dx, Δx, σx, Dy, Δy, σy, no_maxwell, dissipation) ;
  par_grid = (N, J, Box, Val{Order})
  L = [(Box[2d] - Box[2d-1]) for d = 1:D]
  make_periodic!(u,Box_x,N)
  #r = [u[(i-1)*2D+d] for i in 1:N, d in 1:D] # no se como hacerlo funcionar con threads
  r = zeros(Float64,N,D)
  @threads for i in 1:N
              for d in 1:D
              @inbounds   r[i,d] = u[(i-1)*2D+d]
              end
            end
  local_results = zeros(Float64, J[1], J[2], 2, Threads.nthreads())
  idx = ones(Int64, N, 2)
  y = zeros(Float64, N, 2)
  v = zeros(Float64, N, 2)
  n0 = N/prod(J) # dividimos también por el número total de grillas para obtener una densidad independiente del grillado.
  get_indices_and_y_trans!(idx, y, r, J, L)
  v_trans!(Val(D), v, N, u)
  S = get_current(Val(Order), Box_x, J, local_results, idx, y, v)
     
  #@show norm(S)
  Fu = view(u,4N+1:4N+3*prod(J))
  F = reshape(Fu,(3,J...))
  E = F[1:2,:,:]
  B = F[3,:,:]
      
  du .= 0.0
  dFu = view(du,4N+1:4N+3*prod(J))
  dF = reshape(dFu,(3,J...))

  if no_maxwell == false  #take away waves if true
      @threads for i in 1:J[1]
        mul!(view(dF,1,i,:), Dy, view(F,3,i,:),one(eltype(F)))
        mul!(view(dF,3,i,:), Dy , view(F,1,i,:),one(eltype(F)))
        end
        @threads for j in 1:J[2]
        mul!(view(dF,2,:,j), Dx, view(F,3,:,j),-one(eltype(F)))
        mul!(view(dF,3,:,j), Dx, view(F,2,:,j),-one(eltype(F)),one(eltype(F)))
        end
        if dissipation # add dissipation
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
        end
  else
        du[4N+1:4N+3*prod(J)] .= 0.0
  end

  if false #put to false to run maxwell in vacuum
    @threads for j in 1:J[2]
      for i in 1:J[1]
          for l in 1:2
              #dF[l,i,j] +=  S[i,j,l] # particles have negative sign!
              @inbounds dF[l,i,j] +=  S[i,j,l]
          end
      end
    end
  end

      interp = Interpolate(Val(Order), E, B, v, idx, y, J, Box)
      @threads for i in 1:N
        #@inbounds @views v = p2v(u[i*2D-D+1:i*2D])
        # v = p2v(u[range_p(i, D)])
        for d in 1:D
          @inbounds du[(i-1)*2D+d] = v[i,d] # relativistic factor (u is the momentum)
          @inbounds du[i*2D-D+d] = -interp[i, d]
        end
      end
      return du[:]
  end

  function RHS_D_opt(u,t,p_RHSC) #version to optimize
    if nthreads() == 1
      N, J, Box, order, n, S, du, get_density!, get_current, Interpolate,  Dx, Δx, σx, Dy, Δy, σy, maxwell, dissipation = p_RHSC
      par_grid = (N, J, Box, order)
      get_current(u, S, par_grid)
    else
      N, J, Box, order, n, S, du, get_density!, get_current_threads, Interpolate,  Dx, Δx, σx, Dy, Δy, σy, maxwell, dissipation  = p_RHSC
      par_grid = (N, J, Box, order)
      #S = get_current_threads(Val(order), Box_x, u)
    end 
    make_periodic!(u,Box_x,N)
    Fu = view(u,4N+1:4N+3*prod(J))
    F = reshape(Fu,(3,J...))
    E = F[1:2,:,:]
    B = F[3,:,:]
      
    du .= 0.0
    dFu = view(du,4N+1:4N+3*prod(J))
    dF = reshape(dFu,(3,J...))

    if maxwell  #take away waves if false
      @threads for i in 1:J[1]
        mul!(view(dF,1,i,:), Dy, view(F,3,i,:),one(eltype(F)))
        mul!(view(dF,3,i,:), Dy , view(F,1,i,:),one(eltype(F)))
        end
        @threads for j in 1:J[2]
        mul!(view(dF,2,:,j), Dx, view(F,3,:,j),-one(eltype(F)))
        mul!(view(dF,3,:,j), Dx, view(F,2,:,j),-one(eltype(F)),one(eltype(F)))
        end
        if dissipation # take away dissipation
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
        end
      else
        du[4N+1:4N+3*prod(J)] .= 0.0
      end

      @threads for j in 1:J[2]
        for i in 1:J[1]
            for l in 1:2
         #dF[l,i,j] +=  S[l,i,j] # particles have negative sign!
         @inbounds dF[l,i,j] +=  0.0 # S[i,j,l]
            end
        end
      end

      @threads for i in 1:N        
        @inbounds @views v = p2v(u[i*2D-D+1:i*2D])
        # v = p2v(u[range_p(i, D)])
        @inbounds du[range_x(i, D)] = v # relativistic factor (u is the momentum)
        @inbounds du[i*2D-D+1:i*2D] = - Interpolate(order, E, B, v, u[range_x(i, D)], J, Box)
      end
      return du[:]
  end

function RHS_Flux(u,t,par_RHS_Flux) #version to optimize
  if nthreads() == 1
    N, J, Box, dx, order, n, S, du, get_density!, get_current, Interpolate, no_maxwell, par_WENOZ = par_RHS_Flux
    #par_grid = (N, L, J, dx, order)
  else
    Order, N, J, Box, dx, _, n, S, du, get_density!, get_current, Interpolate, no_maxwell, par_WENOZ = par_RHS_Flux
    #par_grid = (N, L, J, dx, order)
  end
  par_grid = (N, J, Box, Val{Order})
  L = [(Box[2d] - Box[2d-1]) for d = 1:D]
  make_periodic!(u,Box_x,N)
  #r = [u[(i-1)*2D+d] for i in 1:N, d in 1:D] # no se como hacerlo funcionar con threads
  r = zeros(Float64,N,D)
  @threads for i in 1:N
              for d in 1:D
              @inbounds   r[i,d] = u[(i-1)*2D+d]
              end
            end
  local_results = zeros(Float64, J[1], J[2], 2, Threads.nthreads())
  idx = ones(Int64, N, 2)
  y = zeros(Float64, N, 2)
  v = zeros(Float64, N, 2)
  n0 = N/prod(J) # dividimos también por el número total de grillas para obtener una densidad independiente del grillado.
  get_indices_and_y_trans!(idx, y, r, J, L)
  v_trans!(Val(D), v, N, u)
  S = get_current(Val(Order), Box_x, J, local_results, idx, y, v)
     
    #@show norm(S)
    Fu = view(u,4N+1:4N+3*prod(J))
    F = reshape(Fu,(3,J...))
    E = F[1:2,:,:]
    B = F[3,:,:]
      
    du .= 0.0
    dFu = view(du,4N+1:4N+3*prod(J)) # I use them for the sources
    dF = reshape(dFu,(3,J...))

    if no_maxwell == false  #take away waves if false
        wenoz!(dFu, F, par_WENOZ, t)
    end

      @threads for j in 1:J[2]
        for i in 1:J[1]
            for l in 1:2
         #dF[l,i,j] +=  S[i,j,l] # particles have negative sign!
         @inbounds dF[l,i,j] +=  S[i,j,l]
            end
        end
      end

      interp = Interpolate(Val(Order), E, B, v, idx, y, J, Box)
      @threads for i in 1:N
        #@inbounds @views v = p2v(u[i*2D-D+1:i*2D])
        # v = p2v(u[range_p(i, D)])
        for d in 1:D
          @inbounds du[(i-1)*2D+d] = v[i,d] # relativistic factor (u is the momentum)
          @inbounds du[i*2D-D+d] = -interp[i, d]
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


function TVD3_Step!(f,y0,t0,h,pf)
  y1 = y0 + h*f(y0,t0,pf) 
  y2 = (3*y0 + y1 + h*f(y1,t0,pf))/4
  y0 .= (y0 + 2*y2 + 2*h*f(y2,t0,pf))/3
end
