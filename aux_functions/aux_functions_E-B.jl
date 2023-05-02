@inline function range_E(N, J)
    D = length(J)
    if D==1
      return  2N+1:2N+J[1]
    elseif D==2
      return 4N+1:4N+D*J[1]*J[2]
    elseif D==3
      return 6N+1:6N+D*J[1]*J[2]*J[3]
    else
      error("Is is not defined for D=$D")
    end
  end
  
  @inline function range_B(N, J)
    D = length(J)
    if D==1
      error("Not defined for D=1") # no magnetic field in this case.
    elseif D==2
      return 4N+D*J[1]*J[2]+1:4N+(D+1)*J[1]*J[2] #here B is one component
    elseif D==3
      return 6N+D*J[1]*J[2]*J[3]+1:6N+2D*J[1]*J[2]*J[3] #here B is 3 components
    else
      error("Is is not defined for D=$D")
    end
  end
  
  @inline function get_E(u,N,J)
    D = length(J)
    return reshape(u[range_E(N,J)],(D,J...))
  end

  
  @inline function get_B(u,N,J)
    D = length(J)
    if D==2
      return reshape(u[range_B(N,J)],J)
    elseif D ==3
      return reshape(u[range_B(N,J)],(D,J...))
    else
      error("Function not defined for D=$D")
    end
  end
  
  @inline function get_Fields!(E,B,u,N,J)
    E = get_E(u,N,J)
    B = get_B(u,N,J)
  end

  

  @inline function ϕ_test(x,x0,Box,r0,p)
    D = length(x)
    desp = zeros(D)
    dx = zeros(D)
    L = [Box[2d] - Box[2d-1] for d in 1:D]
    @assert 0.0 < abs(r0) && abs(r0) < maximum(L)
    for d in 1:D
      if x0[d] - r0 < Box[2d-1]
        desp[d] = r0
      elseif x0[d] + r0 > Box[2d]
        desp[d] = -r0
      else
        desp[d] = 0.0
      end
        dx[d] = (x[d] - x0[d] + desp[d])%L[d] - desp[d]
    end
    r2 = dx'*dx
    r02 = r0^2
    if r2 < r02
      return (r2-r02)^p/r02^p
    else
      return 0
    end
  end

  @inline function ∇ϕ_test(x,x0,Box,r0,p)
    local D = length(x)
    desp = zeros(D)
    L = [Box[2d] - Box[2d-1] for d in 1:D]
    @assert 0.0 < abs(r0) && abs(r0) < maximum(L)
    for d in 1:D
      if x0[d] - r0 < Box[2d-1]
        desp[d] = r0
      elseif x0[d] + r0 > Box[2d]
        desp[d] = -r0
      else
        desp[d] = 0.0
      end
        dx[d] = (x[d] - x0[d] + desp[d])%L[d] - desp[d]
    end
    r2 = dx'*dx
    r02 = r0^2
    if r2 < r0^2
      return 2p*(r2-r02)^(p-1)/r02^p*dx
    else
      return [0,0]
    end
  end

"""
This function is not working properly, it outputs zero...
"""
function div_E!(Div,E,Dx,Dy,J)
    for i in 1:J[1]
        mul!(Div[i,:],Dy,E[2,i,:],one(eltype(Div)))
    end
    for j in 1:J[2]
        mul!(Div[:,j],Dx,E[1,:,j],one(eltype(Div)),one(eltype(Div)))
    end
end

function constraint_test(E,ρ,J,Box, ϕ, ∇ϕ, pars)
  D = length(J)
  dx = differentials(Box,J)
  x_p = [dx[1]*(i-1) + Box[1] for i in 1:J[1]] ;
  y_p = [dx[2]*(i-1) + Box[3] for i in 1:J[2]] ;
  x0, r0, p = pars
  M = zeros(J)
  for i in 1:J[1]
      for j in 1:J[2]
        M[i,j] = ϕ([x_p[i],y_p[j]],x0,Box,r0,p)
      end
  end

  DM = zeros(J...,D)
  for i in 1:J[1]
      for j in 1:J[2]
          for d in 1:D
            DM[i,j,d] = ∇ϕ([x_p[i],y_p[j]],x0,Box,r0,p)[d]
          end
      end
  end

  DivE_∇ϕ_test = 0
  for i in 1:J[1]
      for j in 1:J[2]
          for d in 1:D
              DivE_∇ϕ_test += E[d,i,j]*DM[i,j,d]
          end
      end
  end

  ρf_ϕ_test = 0
  for i in 1:J[1]
      for j in 1:J[2]
          ρf_ϕ_test += ρ[i,j]*M[i,j]
      end
  end

  return DivE_∇ϕ_test, ρf_ϕ_test, abs((DivE_∇ϕ_test - ρf_ϕ_test)/ρf_ϕ_test)
end

  