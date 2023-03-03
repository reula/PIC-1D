function int_mid_point_f(f, par, n::Array, Box_p::Tuple)
    #Vol = volume(Box) 
    D = length(n)
    if D == 1
        dp = (Box_p[2]-Box_p[1])/(n[1]-1)
        pp = [Box_p[1] + dp*(i-1) for i in 1:n[1]]
        F = [f(p,par) for p in pp]
        return (sum(F) - 1//2*(F[1] + F[end]))*dp*2 # the 2 is for the negative side (ONLY FOR SYMMETRIC DISTRIBUTIONS)
    elseif D == 2
        dp1 = (Box_p[2]-Box_p[1])/(n[1]-1)
        dp2 = (Box_p[4]-Box_p[3])/(n[2]-1)
        pp1 = [Box_p[1] + dp1*(i-1) for i in 1:n[1]]
        pp2 = [Box_p[3] + dp2*(i-1) for i in 1:n[2]]
        F = [f([p1,p2],par) for p1 in pp1, p2 in pp2]
        return (sum(F) - 1//2*(sum(F[1,:] + F[end,:]) + sum(F[:,1] + F[:,end])) + 1//4*(F[1,1]+F[1,end]+F[end,1]+F[end,end]))*dp1*dp2
    elseif D == 3
      dp1 = (Box_p[2]-Box_p[1])/(n[1]-1)
      dp2 = (Box_p[4]-Box_p[3])/(n[2]-1)
      dp3 = (Box_p[6]-Box_p[5])/(n[3]-1)
      pp1 = [Box_p[1] + dp1*(i-1) for i in 1:n[1]]
      pp2 = [Box_p[3] + dp2*(i-1) for i in 1:n[2]]
      pp3 = [Box_p[5] + dp3*(i-1) for i in 1:n[3]]
      F = [f([p1,p2,p3],par) for p1 in pp1, p2 in pp2, p3 in pp3]
      Int = sum(F) 
      Int += - 1//2*(sum(F[1,:,:] + F[end,:,:]) + sum(F[:,1,:] + F[:,end,:]) + sum(F[:,:,1] + F[:,:,end])) 
      Int += + 1//4*sum(F[1,1,:]+F[1,end,:]+F[end,1,:]+F[end,end,:])
      Int += + 1//4*sum(F[1,:,1]+F[1,:,end]+F[end,:,1]+F[end,:,end]) 
      Int += + 1//4*sum(F[:,1,1]+F[:,1,end]+F[:,end,1]+F[:,end,end])
      Int += - 1//8*(F[1,1,1]+F[1,1,end]+F[1,end,1]+F[1,end,end]+F[end,1,1]+F[end,1,end]+F[end,end,1]+F[end,end,end])
      Int = Int*dp1*dp2*dp3
      return Int
    else 
      error("more than 3 dimensions is not implemented")
    end
  end