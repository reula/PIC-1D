@inline function range_E(N, J)
    D = length(J)
    if D==1
      return  N+1:N+J[1]
    elseif D==2
      return N+1:N+D*J[1]*J[2]
    elseif D==3
      return N+1:N+D*J[1]*J[2]*J[3]
    else
      error("Is is not defined for D=$D")
    end
  end
  
  @inline function range_B(N, J)
    D = length(J)
    if D==1
      error("Not defined for D=1") # no magnetic field in this case.
    elseif D==2
      return N+D*J[1]*J[2]+1:N+(D+1)*J[1]*J[2] #here B is one component
    elseif D==3
      return N+D*J[1]*J[2]*J[3]+1:N+2D*J[1]*J[2]*J[3] #here B is 3 components
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


  
  