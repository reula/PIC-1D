using FileIO
using JLD2
using Base.Threads
using Distributions
using LaTeXStrings
using Printf
using ArraysOfArrays
using SummationByPartsOperators
using LinearAlgebra
using FFTW
FFTW.set_provider!("mkl")
using BenchmarkTools
using Profile
using ProfileView
println("nthreads = $(nthreads())")

@show pwd()
#cd("./Tests")
include("../aux_functions/aux_functions.jl")
include("../aux_functions/inidat_aux_functions.jl")


## initial data configurations

trys = false
thermal = false
weibel = false
weibel_norel = false
#trys = true
thermal = true
#weibel = true
#weibel_norel = true

const order = 5
const D = 2


# Particle numbers
const N_exp = 5 #7#6#5 #6
N = 10^(N_exp)

if  thermal
    #run_name = "thermal_norel_"
    run_name = "thermal_rel_"
    #data_name = "par_dis_norel_thermal_"
    data_name = "par_dis_rel_thermal_"
    J = (50,50)
    Box_x = (0.0,1.0,0.0,1.0) 
    nm = [1,1]
    Box_array = [i for i in Box_x]
    k = 2π*nm./(Box_array[2:2:end]-Box_array[1:2:end-1]) # this here is k, not \hat{k}
    alpha_exp = 1 # 2 3 8
    α = 10.0^(-alpha_exp) # 0.1 0.001
    par_f_x = (α, nm, Box_x)
    exp_Theta = 3
    θ = 10.0^(-exp_Theta)
    p_max = sqrt((1+10*θ)^2-1)
    Box_p = (-p_max,p_max,-p_max,p_max)
    par_f_p = (θ,D)
    #par_f_p_1 = (θ,1,D) #to compute the norm of f_p
    #norm = int_mid_point_f(f_p_rel, par_f_p_1, [20_000,20_000], Box_p)
    #par_f_p_rel = (θ,norm,D)
    par_init = (α, nm)
    data_name = data_name * "nm_[$(nm[1]),$(nm[2])]_"
    data_name = data_name * "alp$(alpha_exp)_N$(N_exp)_Th$(exp_Theta)"
    run_name = run_name * "Box_$(Box_x[2]-Box_x[1])x$(Box_x[4]-Box_x[3])_J_$(J[1])x$(J[2])_nm_[$(nm[1]),$(nm[2])]__Th$(exp_Theta)_alp$(alpha_exp)"
elseif weibel_norel
    run_name = "weibel_norel_"
    data_name = "par_dis_norel_weibel_"
    J = (100,100)
    Box_x = (0.0,1.0,0.0,1.0) 
    Ax = 1/2 # velocity anisotropy in the x direction
    alpha_exp = 8 # 8 so that is extremely small i.e. zero
    α = 10.0^(-alpha_exp) # 0.1 0.001
    par_f_x = (α, [0.0,0.0], Box_x)
    exp_Theta = 3
    θ = 10.0^(-exp_Theta)
    p_max = sqrt((1+10*θ)^2-1)
    Box_p = (-p_max,p_max,-p_max,p_max)
    par_f_p = (θ,D,Ax)
    par_init = (α, Ax)
    data_name = data_name * "Ax_(1d2)_"
    data_name = data_name * "alp$(alpha_exp)_N$(N_exp)_Th$(exp_Theta)"
    run_name = run_name * "Box_$(Box_x[2]-Box_x[1])x$(Box_x[4]-Box_x[3])_J_$(J[1])x$(J[2])_Ax_(1d2)_Th$(exp_Theta)_alp$(alpha_exp)"
elseif weibel 
    run_name = "weibel_"
    data_name = "par_dis_rel_weibel_"
    J = (50,50)::NTuple{2,Int64} #no parece funcionar el dar el tipo NTuple, al menos no así.
    Box_x = (0.0,1.0,0.0,1.0) 
    Ax = 25 # velocity anisotropy in the x direction
    alpha_exp = 8 # 8 so that is extremely small i.e. zero
    α = 10.0^(-alpha_exp) # 0.1 0.001
    par_f_x = (α, [0.0,0.0], Box_x)
    exp_Theta = 3
    θ1 = 10.0^(-exp_Theta)
    θ2 = Ax * θ1
    p_max = 10*θ2
    Box_p = (-p_max,p_max,-p_max,p_max)
    par_init = (α, Ax)
    @show pars_f = (θ1,θ2,Ax)
    data_name = data_name * "Ax_$(Ax)_"
    data_name = data_name * "alp$(alpha_exp)_N$(N_exp)_Th$(exp_Theta)"
    run_name = run_name * "Box_$(Box_x[2]-Box_x[1])x$(Box_x[4]-Box_x[3])_J_$(J[1])x$(J[2])_Ax_$(Ax)_Th$(exp_Theta)_alp$(alpha_exp)"
end



@show data_name 

# Evolution parameters
exp_t = 1 #0# 1 #2
t = 0.0
t_i = 0.0
t_f = 10.0^(exp_t)
M = 5001 #1001 #11 #16001# 4001 # 81 # 2001 # time steps 
M_g = 101 #8001 #400 + 1 #number of outputs, starting from the initial data
dx = differentials(Box_x,J)
@show dx_min = minimum(dx)
@show dt = (t_f-t_i)/(M-1)
@show CFL = dt/dx_min


# for plotting and other things
x_p = [dx[1]*(i-1) for i in 1:J[1]] ;
y_p = [dx[2]*(i-1) for i in 1:J[2]] ;

Dx = periodic_derivative_operator(derivative_order=1, accuracy_order=6, xmin=Box_x[1], xmax=Box_x[2], N=J[1])
Dy = periodic_derivative_operator(derivative_order=1, accuracy_order=6, xmin=Box_x[3], xmax=Box_x[4], N=J[2])
Δx = dissipation_operator(Dx;
                     #mode=D.coefficients.mode
                     #,mode=ThreadedMode()
                     )
Δy = dissipation_operator(Dy;
                     #mode=D.coefficients.mode
                     #,mode=ThreadedMode()
                     )
const σx = 0.0 #1.0 #dissipation strength
const σy = 0.0 #1.0 #dissipation strength
dissipation = false
maxwell = true

@show par_evolv = (t_i, t_f, M, M_g, dt)
@show par_grid = (N, J, Box_x, order)
#@show pars_f = (θ1,θ2,)



println("t_f = $(t_f), M = $M, dt = $(dt), exp_Theta = $(exp_Theta)")



run_name = run_name * "tf_$(convert(Int,10*t_f))_N$(N_exp)_M$(M)_o$(order)"
println(run_name)
println(data_name)

par_dis, data_name, pars, par_f_x, Box_x, par_f_p, Box_p = retrieve_initial_data_D("../Initial_Distributions/" * data_name * ".jld2")

get_density_2D_trans = Density2DTrans(N, J)
n = get_density_2D_trans(Val(order), Box_x, par_dis)
get_current_2D_trans = Current2DTrans(N, J)
S = get_current_2D_trans(Val(order), Box_x, par_dis)

B0 = 2.0 #initial magnetic field

B = [B0 for i in 1:J[1], j in 1:J[2]]

run_name = run_name * "_B0_$(convert(Int,B0))"

E = Array{Float64,3}(undef,(2,J...))
#ρ = n.-1.0/prod(J)
ρ = n.-1.0
@show sum(ρ)
get_E!(E,ρ,Box_x);

average_outputs = false # detailed output (for runs which are too long)
full_outputs = false
animation = false # to run animations
phase_space_show = false #show phase space diagrams
#phase_space_show = true
remote_server = false # if run in remote_server as a script avoid all plots 
save_plots = false # put true is you want so save your plots.
full_outputs = true
animation = false

u = Vector{Float64}(undef, 4N + 3*prod(J)); # contains r, v and E and B
du = Vector{Float64}(undef, 4N + 3*prod(J)); # contains r, v and E
du_ref = Vector{Float64}(undef, 4N + 3*prod(J)); # contains r, v and E
u[1:4N] = par_dis 
Fu = view(u,4N+1:4N+3*prod(J))
F = reshape(Fu,3,J...)
F[1:2,:,:] = E
F[3,:,:] = B;

Energy_Ks, Energy_Es = get_energy_rel(u,(Box_x, N, J))
println("Energies= EK = $(Energy_Ks), E = $(Energy_Es)")

Coordinate_test(u,Box_x,N)

if nthreads() > 1
    #TS = zeros(Float64, (2,J...,nthreads()))
    #p_RHS_D = (N, J, Box_x, order, n, S, du, get_density_2D!, get_current_threads_2D!, Interpolate_EBv_1, TS, Dx, Δx, σx, Dy, Δy, σy) ;
    p_RHS_D = (N, J, Box_x, order, n, S, du, get_density_2D!, get_current_slim, Interpolate_All_EBv_1_slim, Dx, Δx, σx, Dy, Δy, σy, maxwell, dissipation) ;
    p_RHS_D_ref = (N, J, Box_x, order, n, S, du_ref, get_density_2D!, get_current_2D_trans, Interpolate_EBv_1, Dx, Δx, σx, Dy, Δy, σy, maxwell, dissipation) ;
else
    p_RHS_D = (N, J, Box_x, order, n, S, du, get_density_2D!, get_current_rel_2D!, Interpolate_EBv_1, Dx, Δx, σx, Dy, Δy, σy, maxwell, dissipation) ;
    
end

val_order = Val(order)

RHS_D(u,0.0,p_RHS_D_ref); # here we save the output before changes
#RHS_D_opt(du,u,p_RHS_D);
RHS_D_slim!(val_order,u,0.0,p_RHS_D);
#@profile RHS_D_slim(u,0.0,p_RHS_D);
#VSCodeServer.@profview RHS_D_slim!(val_order,u,0.0,p_RHS_D);

@show norm(u)
@show norm(du_ref), norm(du)
@show norm(du_ref - du)
#@btime RHS_D_slim!($val_order,$u,0.0,$p_RHS_D);


# N = 10^5
# J = (50,50)
# threads = 2, RHS_D($du,$u,$p_RHS_D), 1.959 s (46123146 allocations: 1.91 GiB)
# threads = 2, RHS_D(du,u,p_RHS_D), 1.917 s (46123146 allocations: 1.91 GiB)

#jupyter_CCAD
# with no dissipation
# N = 10^6
# J = (50,50)
# threads = 2,  12.275 s (461023153 allocations: 19.04 GiB)
# threads = 20,  9.567 s (461024065 allocations: 19.05 GiB)
# threads = 40, 14.994 s (461025064 allocations: 19.05 GiB)

# with no dissipation nor maxwell
# N = 10^6
# J = (50,50)
# threads = 2,  12.345 s (461022726 allocations: 19.04 GiB)
# threads = 20,  9.128 s (461023429 allocations: 19.05 GiB)
# threads = 40, 17.194 s (461024189 allocations: 19.05 GiB)
# 
# with no dissipation nor maxwell nor S, only interpolations
# N = 10^6
# J = (50,50)
# threads = 2, 11.652 s (459017648 allocations: 18.88 GiB)
# threads = 20, 9.225 s (459017866 allocations: 18.88 GiB)
# threads = 40, 13.880 s (459018110 allocations: 18.88 GiB)

# slim version with no dissipation nor maxwell only S, and interpolations
# N = 10^6
# J = (50,50)
# threads = 2, 12.836 s (460000118 allocations: 19.06 GiB)
# threads = 20, 8.439 s (460000940 allocations: 19.06 GiB)
# threads = 20, 9.012 s (460021069 allocations: 19.06 GiB)
# threads = 40, 13.963 s (460001826 allocations: 19.06 GiB)



# B0 = 1   (norm(du_ref), norm(du)) = (20.08351535598122, 2.4834756109661114)
# B0 = 0   (norm(du_ref), norm(du)) = (14.30248902971554, 2.4576354850116586)
# B0 = 2   (norm(du_ref), norm(du)) = (31.613679521935694, 2.5580533146069464)