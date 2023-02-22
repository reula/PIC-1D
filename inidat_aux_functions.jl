include("aux_functions.jl")


"""The following routine returns a random sample according to a given distribution function.
The algorithm used to achieve this is called the rejection method.
f is the distribution function and p are its parameters. We need also
f_max which is a function of p, and the interval of the sampling region
int = (x_min,x_max)
"""
function random_sampling_from_distribution(f,f_max,par_f,interval)
 (x_max,x_min) = interval
 fmax = f_max(par_f)
 x = x_min + (x_max - x_min) * (rand());
 #  Accept/reject value
 f_v = f(x,par_f)
 x_t = fmax * rand()
 if (x_t > f_v) return random_sampling_from_distribution(f,f_max,par_f,interval)
 else return x
 end
end

function random_sampling_from_distribution(f,f_max,par_f,Box::Tuple)
  fmax = f_max(par_f)
  D = length(Box) ÷ 2
  x = zeros(D)
  for i in 1:D
    x[i] = Box[i*2-1] + (Box[i*2] - Box[i*2-1]) * (rand());
  end
  #  Accept/reject value
  f_v = f(x,par_f)
  x_t = fmax * rand()
  if (x_t > f_v) return random_sampling_from_distribution(f,f_max,par_f,Box)
  else return x[:]
  end
 end
"""
build_initial_data(data_name::String, par_grid, f_r, f_r_max, f_p, f_p_max, par_f_r, par_f_p)
builds complete initial data for distributions which are product of two distributions.
It uses random_sampling_from_distribution()
It stores the data in Initial_Distributions/
data_name = string with the name of the initial_data_set
par_grid = grid parameters, 
f_r = Distribution function for the position, 
f_r_max = Maximum of Distribution function for the position, 
f_p = Distribution funtions for momentum, 
f_p_max = Maximum of Distribution function for momentum, 
par_f_r = parameters for f_r, 
par_f_p = parameters for f_p
"""
function build_initial_data(data_name::String, pars, f_x, f_x_max, par_f_x, interval_x, f_p, f_p_max, par_f_p, interval_p)
    N, = pars
    vp = zeros(N÷2)
    r = zeros(N)
    
    
    for i in 1:N
    r[i] = random_sampling_from_distribution(f_x,f_x_max,par_f_x,interval_x)
    end
    

    for i in 1:N÷2
        vp[i] = random_sampling_from_distribution(f_p,f_p_max,par_f_p,interval_p)
    end

    par_dis = [r;vp;-vp];


    file_name = "Initial_Distributions/" * data_name * ".jld2"
    
    field_name = "par_dis"
    run_pars = Dict("data_name" => data_name, "pars" => pars, "par_f_x" => par_f_x,"par_f_p" => par_f_p)
    save(file_name, run_pars)
    
    #save("Initial_Distributions/" * data_name * ".jld2", "field_name", par_dis)
    jldopen(file_name, "a+") do file
        file[field_name] = par_dis;
    end
  end

function build_initial_data_D(data_name::String, pars, f_x, f_x_max, par_f_x, Box_x, f_p, f_p_max, par_f_p, Box_p)
      N, = pars
      D = length(Box_x)÷2
      par_dis = zeros(D*N*2)

      @show D
      
      # first the space distribution 

      for i in 1:N
        par_dis[(i-1)*2*D+1:(i-1)*2*D+D] = random_sampling_from_distribution(f_x,f_x_max,par_f_x,Box_x)
      end
  
      # We set part of the distribution in antisymmetric form so that the total momentum vanishes.
      for i in 1:2D:D*N
          par_dis[i+D:(i-1)+2D] = random_sampling_from_distribution(f_p,f_p_max,par_f_p,Box_p)
          par_dis[D*N + i+D:D*N + i-1+2D] = - par_dis[i+D:(i-1)+2D]
      end
  
      file_name = "Initial_Distributions/" * data_name * ".jld2"
      
      field_name = "par_dis"
      run_pars = Dict("data_name" => data_name, "pars" => pars, "par_f_x" => par_f_x, "Box_x" => Box_x, "par_f_p" => par_f_p, "Box_p" => Box_p)
      save(file_name, run_pars)
      
      #save("Initial_Distributions/" * data_name * ".jld2", "field_name", par_dis)
      jldopen(file_name, "a+") do file
          file[field_name] = par_dis;
      end
end

function retrieve_initial_data_D(file_name::String)
  data = load(file_name)
  run_name = data["data_name"]
  pars = data["pars"]
  par_f_x = data["par_f_x"]
  Box_x = data["Box_x"]
  par_f_p = data["par_f_p"]
  Box_p = data["Box_p"]

  return data["par_dis"], run_name, pars, par_f_x, Box_x, par_f_p, Box_p
end

function retrieve_initial_data(file_name::String)
  data = load(file_name)
  run_name = data["data_name"]
  pars = data["pars"]
  par_f_x = data["par_f_x"]
  par_f_p = data["par_f_p"]

  return data["par_dis"], run_name, pars, par_f_x, par_f_p
end

# momentum distributions

""" to normalize the momentum distribution first run with norm=1, 
and then calculate the norm, and use it again to get the correct normalization. 
"""
function norm_f_p_rel(f_p_rel, par_f_p, n, p_max)
  dp = p_max/(n-1) 
  p = [dp*(i-1) for i in 1:n]
  F(p) = f_p_rel(p,par_f_p)
  return (sum(F.(p)) - 0.5*(F(0) + F(p_max)))*dp*2 # the 2 is for the negative side (ONLY FOR SYMMETRIC DISTRIBUTIONS)
end

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

f_p_rel(p,(θ,norm)) = exp((1 - sqrt(1+p^2))/θ) / sqrt(θ*π*2) / norm 
f_p_rel(p::Array,(θ,norm,D)) = exp((1 - sqrt(1+p'*p))/θ) / sqrt(θ*π*2)^(length(p)) / norm 
f_p_rel_max((θ,norm,D)) = 1 / sqrt(θ*π*2)^D / norm 

f_p(p,(θ,)) = exp(- p^2/θ/2) / sqrt(θ*π*2) 
f_p(p::Array,(θ,D)) = exp(- (p'*p)/θ/2) / sqrt(θ*π*2)^D
f_p_max((θ,D)) = 1 / sqrt(θ*π*2)^D




"""The following routine returns a random velocity distributed on a double Maxwellian distribution function 
corresponding to two counter-streaming beams. The algorithm used to achieve this is called the rejection method, 
and will be discussed later in this course.

  // Function to distribute electron velocities randomly so as 
  // to generate two counter propagating warm beams of thermal
  // velocities unity and mean velocities +/- vb.
  // Uses rejection method.
"""
function distribution_stream(vb)
  #Generate random v value
  fmax = 0.5 * (1. + exp(-2.0 * vb * vb));
  vmin = - 5.0 * vb;
  vmax = + 5.0 * vb;
  v = vmin + (vmax - vmin) * (rand());

  #  Accept/reject value
  f = 0.5 * (exp(-(v - vb) * (v - vb) / 2.0) + exp(-(v + vb) * (v + vb) / 2.0))
  x = fmax * rand()
  if (x > f) return distribution_stream(vb)
  else return v
  end
end

# space distributions 

function f_x(x::Float64,par_f_x) 
    α, mn, L = par_f_x
    k = 2π*mn/L
    return (1 + α *cos(k'*x))/L
end

function f_x_max(par_f_x)
    if typeof(par_f_x[3]) == typeof(1.0)
        α, mn, L = par_f_x
        return (1+abs(α))/L
    else 
      α, mn, Box = par_f_x
        return (1+abs(α))/volume(Box)
    end
end

"""
Example:
Box = (0.0,3.0,-1.0,1.0)
m = [1,3]
α = [0.1,0.2]
par = (α, m, Box)
f_x([3,4],par)
> 0.75
"""
function f_x(x::Array,par_f_x)
  α, m, Box = par_f_x
  k = zeros(length(m))
  Box_array = [i for i in Box]
  L = Box_array[2:2:end] - Box_array[1:2:end-1]
  k = 2π*m./L
  return (1 + α*cos(k'*x))/volume(Box)
end



    


