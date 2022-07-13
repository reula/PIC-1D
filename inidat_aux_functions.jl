
"""The following routine returns a random sample according to a given distribution function.
The algorithm used to achieve this is called the rejection method.
f is the distribution function and p are its parameters. We need also
f_max which is a function of p, and the interval of the sampling region
int = (x_min,x_max)
"""
function random_sampling_from_distribution(f,f_max,p,interval)
 (x_max,x_min) = interval
 fmax = f_max(p)
 x = x_min + (x_max - x_min) * (rand());
 #  Accept/reject value
 f_v = f(x,p)
 x_t = fmax * rand()
 if (x_t > f_v) return random_sampling_from_distribution(f,f_max,p,interval)
 else return x
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
function build_initial_data(data_name::String, pars, f_x, f_x_max, par_f_x, f_p, f_p_max, par_f_p)
    (N, L) = pars
    vp = zeros(N÷2)
    r = zeros(N)
    
    interval_x = (0.0, L)
    for i in 1:N
    r[i] = random_sampling_from_distribution(f_x,f_x_max,par_f_x,interval_x)
    end
    
    p_max = 1 + 10*θ #p_max = sqrt(20*θ)
    interval_p = (0.0, p_max)

    for i in 1:N÷2
        vp[i] = random_sampling_from_distribution(f_p,f_p_max,par_f_p,interval_p)
    end

    par_dis = [r;vp;-vp];


    file_name = "Initial_Distributions/" * data_name * ".jld2"
    
    field_name = "par_dis"
    run_pars = Dict("data_name" => data_name, "pars" => (N,L), "par_f_x" => par_f_x,"par_f_p" => par_f_p)
    save(file_name, run_pars)
    
    #save("Initial_Distributions/" * data_name * ".jld2", "field_name", par_dis)
    jldopen(file_name, "a+") do file
        file[field_name] = par_dis;
    end


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

f_p_rel(p,(θ,)) = exp((1 - sqrt(1+p^2))/θ) 
f_p_rel_max(θ) = 1


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

function f_x(x,par_f_x) 
    α, k = par_f_x
    return 1 + α *cos(k*x)
end

function f_x_max(par_f_x)
    α, k = par_f_x
    return 1+α
end



