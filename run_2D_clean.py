# Cell 0
#using DifferentialEquations
using Plots
using Statistics
using FFTW
#FFTW.set_provider!("mkl")
#import Pkg; Pkg.add("FileIO")
using FileIO
using JLD2
using Base.Threads
using Distributions
#Pkg; Pkg.add("DistributedArrays")
println("nthreads = $(nthreads())")
using Printf
#import Pkg; Pkg.add("IJuliaBell")
#using IJuliaBell
using ArraysOfArrays
using SummationByPartsOperators
using LinearAlgebra
using Revise

# Cell 2
includet("aux_functions/aux_functions.jl")
includet("aux_functions/inidat_aux_functions.jl")
includet("aux_functions/choques_utils.jl") 
# todas las funciones necesarias


# Cell 5
## initial data configurations

trys = false
thermal = false
weibel = false
weibel_norel = false
damped = false
#trys = true
thermal = true
#weibel = true
#weibel_norel = true
#damped = true

## simulation configurations
no_maxwell = false #adds the Maxwell equations to the system
sbp = false
holbo = false
lanczos = false
weno = false

## derivatives approximation
#holbo = true
sbp = true
#lanczos = true
#weno = true
#no_maxwell = true

# Interpolator for E and B
interpolation = 1
#interpolation = 2

# Cell 6
const order = 5 #particles smoothness
const D = 2 # Dimension
const factor = 2 #factor between dx for grid and dx for shape. A factor bigger than 1 makes the shapes extend over a bigger region of the grid.


# Particle numbers
const N_exp = 5 #5 #6 #7 #6 #7 #7#6#5 #6
const N = 10^(N_exp)

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
    run_name = "weibel_I2"
    #data_name = "par_dis_rel_weibel_"
    data_name = "par_dis_rel_weibel_MN_"
    J = (128,128)::NTuple{2,Int64} #(64,64)::NTuple{2,Int64} #no parece funcionar el dar el tipo NTuple, al menos no así.
    Box_x = (0.0,15.0,0.0,15.0) 
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
elseif damped
    run_name = "damped_I2_"
    data_name = "par_dis_rel_damped_"
    J = (930,100)::NTuple{2,Int64}
    Box_x = (0.0,7.455,0.0,1.0) # corner box coordinates
    Box_array = [i for i in Box_x]
    nm = [15,0] # wave numbers
    k = 2π*nm./(Box_array[2:2:end]-Box_array[1:2:end-1]) # this here is k, not \hat{k}
    alpha_exp = 2 # 3 8
    α = 10.0^(-alpha_exp) # 0.1, 0.01 0.001
    par_f_x = (α, nm, Box_x)
    exp_Theta = 3
    θ = 10.0^(-exp_Theta)
    p_max = sqrt((1+10*θ)^2-1) #maximum of the momentum box
    Box_p = (-p_max,p_max,-p_max,p_max)
    par_f_p = (θ,D)
    #par_f_p_1 = (θ,1,D) #to compute the norm of f_p
    #norm = int_mid_point_f(f_p_rel, par_f_p_1, [20_000,20_000], Box_p)
    #par_f_p_rel = (θ,norm,D)
    data_name = data_name * "nm_[$(nm[1]),$(nm[2])]_"
    data_name = data_name * "alp$(alpha_exp)_N$(N_exp)_Th$(exp_Theta)"
    run_name = run_name * "Box_$(Box_x[2]-Box_x[1])x$(Box_x[4]-Box_x[3])_J_$(J[1])x$(J[2])_Th$(exp_Theta)_alp$(alpha_exp)"
end



@show data_name 



# EVOLUTION PARAMETERS

exp_t = 1 #0# 1 #2
t = 0.0
t_i = 0.0
t_f = 1.0*10.0^(exp_t)
M = 10001 #10001 #10001 #1001 #11 #16001# 4001 # 81 # 2001 # time steps 
M_g = 1000 #8001 #400 + 1 #number of outputs, starting from the initial data
dx = differentials(Box_x,J)
@show dx_min = minimum(dx)
@show dt = (t_f-t_i)/(M-1)
@show CFL = dt/dx_min


# for plotting and other things
x_p = [dx[1]*(i-1) for i in 1:J[1]] ;
y_p = [dx[2]*(i-1) for i in 1:J[2]] ;

if sbp
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
    dissipation = false #adds KO dissipation to the system
    run_name = run_name * "tf_$(convert(Int,10*t_f))_N$(N_exp)_M$(M)_o$(order)_sbp"

elseif holbo
    stencil_width=11
    Dx = periodic_derivative_operator(Holoborodko2008(); derivative_order=1, accuracy_order=4, stencil_width=stencil_width, xmin=Box_x[1], xmax=Box_x[2], N=J[1])
    Dy = periodic_derivative_operator(Holoborodko2008(); derivative_order=1, accuracy_order=4, stencil_width=stencil_width, xmin=Box_x[3], xmax=Box_x[4], N=J[2])
    dissipation = false
    Δx = 0.0
    Δy = 0.0
    σx = 0.0
    σy = 0.0
    run_name = run_name * "tf_$(convert(Int,10*t_f))_N$(N_exp)_M$(M)_o$(order)_holbo_stencil$(stencil_width)"

elseif lanczos
    stencil_width=11
    Dx = periodic_derivative_operator(LanczosLowNoise(); derivative_order=1, accuracy_order=4, stencil_width=stencil_width, xmin=Box_x[1], xmax=Box_x[2], N=J[1])
    Dy = periodic_derivative_operator(LanczosLowNoise(); derivative_order=1, accuracy_order=4, stencil_width=stencil_width, xmin=Box_x[3], xmax=Box_x[4], N=J[2])
    dissipation = false
    Δx = 0.0
    Δy = 0.0
    σx = 0.0
    σy = 0.0
    run_name = run_name * "tf_$(convert(Int,10*t_f))_N$(N_exp)_M$(M)_o$(order)_lanczos_stencil$(stencil_width)"

elseif weno
    function Flux_x!(F,u,F_pars)
        Ea = u[1:2]
        B = u[3]
        F[1] = 0.0
        F[2] = B 
        F[3] = Ea[2]
    end
    
    function Flux_y!(F,u,F_pars)
        Ea = u[1:2]
        B = u[3]
        F[1] = -B
        F[2] = 0.0
        F[3] = -Ea[1]
    end
    
    SpeedMax(U,eq_pars) = 1.0
    
    N_FIELDS = 3
    eqpars = (1.0,)
    one_over_dx = 1.0./dx
    (t_i, t_f) = (0.0, 4.0)
    tspan = (t_i, t_f)
    auxvectorsWENOZ = createWENOZvectors(N_FIELDS)
    parWENOZ = (eqpars, one_over_dx, J, N_FIELDS, (Flux_x!,Flux_y!), SpeedMax, auxvectorsWENOZ)
    run_name = run_name * "_weno"
    
elseif no_maxwell
    Dx = 0
    Dy = 0
    dissipation = false
    Δx = 0.0
    Δy = 0.0
    σx = 0.0
    σy = 0.0
    run_name = run_name * "_no_maxwell"
else 
    println("No SBP operator selected")
end




@show par_evolv = (t_i, t_f, M, M_g, dt)
@show par_grid = (N, J, Box_x, order)
#@show pars_f = (θ1,θ2,)



println("t_f = $(t_f), M = $M, dt = $(dt), exp_Theta = $(exp_Theta)")

#interpolation defined in above cell

if interpolation == 1 # we use the _S because it preserves energy
    if nthreads()>1
        interp = Interpolate_All_EBv_1_slim_S
        run_name = run_name * "_int1_threads"
    else
        interp = Interpolate_EBv_1_S
        run_name = run_name * "_int1"
    end
elseif interpolation == 2 #do not use, they do not preserve energy
    if nthreads()>1
        interp = Interpolate_All_EBv_2_slim_W
        run_name = run_name * "_int2_threads"
    else
        interp = Interpolate_EBv_2_W
        run_name = run_name * "_int2"
    end
else
    error("No interpolation selected")
end

println(run_name)
println(data_name)


# Cell 7
par_dis, data_name, pars, par_f_x, Box_x, par_f_p, Box_p = retrieve_initial_data_D("Initial_Distributions/" * data_name * ".jld2")
#@show data_name_from_inidat, pars, par_f_x, Box_x, par_f_p, Box_p


# Cell 9
plot(histogram2d(par_dis[1:2D:end],par_dis[2:2D:end], title="space distribution", bins=50, c=:thermal))

# Cell 10
plot(histogram2d(par_dis[3:4:end],par_dis[4:4:end] )
,aspectratio=1
,title = "momentum histogram"
)

# Cell 12
T = get_temperature_rel_D(par_dis,N,D)

# Cell 14
#includet("aux_functions/aux_functions.jl")
#par_grid = (N, Box_x, J, 5)
#n = zeros(J)
#get_density_2D!(par_dis, n, par_grid)

get_density_2D_trans = Density2DTrans(N, J)
n = get_density_2D_trans(Val(order), Box_x, par_dis, factor=factor);


# Cell 15
#plot_matrix(n)
surface(y_p,x_p, n[:,:])

# Cell 16
plot(n[:,2])

# Cell 18
sum(n)/prod(J)

# Cell 20
par_grid = (N, J, Box_x, 5)
#S = [0.0::Float64 for l in 1:2, i in 1:J[1], j in 1:J[2]]
#S_new  = [0.0::Float64 for i in 1:J[1], j in 1:J[2], l in 1:2]
#TS = zeros(Float64,(2,J...,nthreads()))

#par_current_threads_S = (par_grid, TS)

#get_current_threads_2D!(par_dis, S, par_current_threads_S)
#S_new = get_current_2D_trans(Val(order), N, J, Box_x, par_dis)


get_current_2D_trans = Current2DTrans(N, J)
S = get_current_2D_trans(Val(order), Box_x, par_dis, factor=factor)
@show sum(S[:,:,1])

# Cell 21
surface(y_p,x_p, S[:,:,1])

# Cell 23
# We prescrive a constant magnetic field and the E is just a solution from Poisson's equations. 

B0 = 0.0 #initial magnetic field

B = [B0 for i in 1:J[1], j in 1:J[2]]

run_name = run_name * "_B0_$(convert(Int,trunc(B0)))"

# Cell 24
E = Array{Float64,3}(undef,(2,J...))
#ρ = n.-1.0/prod(J)
ρ = n.-1.0
@show sum(ρ)
get_E!(E,ρ,Box_x);

if false #testing evolution in vacuum with a smooth electric field
    #includet("aux_functions/aux_functions_E-B.jl")
    r0 = 30*minimum(abs.(dx))
    x0 = [2.0,2.0]
    p = 6
    #x_p = [dx[1]*(i-1) + Box_x[1] for i in 1:J[1]] ;
    #y_p = [dx[2]*(i-1) + Box_x[3] for i in 1:J[2]] ;
    for i in 1:J[1]
        for j in 1:J[2]
            E[1,i,j] = ∇ϕ_test([x_p[i],y_p[j]],x0,Box_x,r0,p)[2]
            E[2,i,j] = -∇ϕ_test([x_p[i],y_p[j]],x0,Box_x,r0,p)[1]
        end
    end
end

# Cell 25
surface(y_p,x_p, E[2,:,:])

# Cell 26
plot(x_p, E[1,:,5])

# Cell 27
Err = 0.0
    Max = 0.0
    r0 = 0.8
    M_test = 1000
    x0=zeros(M_test,D)
    err = zeros(M_test)
    rho = zeros(M_test)
    erab = zeros(M_test)
    diV = zeros(M_test)
    nn = zeros(J...)
for i in 1:M_test
    x0[i,:]=[Box_x[1] + (Box_x[2]-Box_x[1])*rand(), Box_x[3] + (Box_x[4]-Box_x[3])*rand()]
    #diV[i], rho[i], err[i], erab[i] = constraint_test(E, n .- 1.0, J, Box_x, ϕ_test, ∇ϕ_test, (x0[i,:], r0, 6))
    diV[i], rho[i], err[i], erab[i] = constraint_test(E, nn, J, Box_x, ϕ_test, ∇ϕ_test, (x0[i,:], r0, 6))
    global Err = Err + abs(err[i])
    if err[i] > Max
        @show Max = err[i]
        @show x0[i,:]
    end
end
println("error = $(Err/M_test)")
println("Max = $(Max)")

# Cell 28

p1 = Plots.scatter(x0[:,1],x0[:,2],marker_z=-err[:], markersize=10, alpha=0.1, label="-err", aspectratio=1)
p2 = Plots.scatter(x0[:,1],x0[:,2],marker_z=-1.0./abs.(rho[:])./100, markersize=3, alpha=1, label="rho", aspectratio=1)
p3 = Plots.scatter(x0[:,1],x0[:,2],marker_z=-erab[:], markersize=10, alpha=0.1, label="-erab", aspectratio=1)
plot(p1,p2,p3, layout = (2, 2))

# Cell 29
histogram(erab)

# Cell 30
average_outputs = false # detailed output (for runs which are too long)
full_outputs = false
animation = false # to run animations
phase_space_show = false #show phase space diagrams
#phase_space_show = true
remote_server = false # if run in remote_server as a script avoid all plots 
save_plots = false # put true is you want so save your plots.
#full_outputs = true
average_outputs = true
#animation = true

# Cell 32
u = Vector{Float64}(undef, 4N + 3*prod(J)); # contains r, v and E and B

du = Vector{Float64}(undef, 4N + 3*prod(J)); # contains r, v and E

u[1:4N] = par_dis 

Fu = view(u,4N+1:4N+3*prod(J))

F = reshape(Fu,3,J...)

#F = view(u[4N+1:end],3,J...)
F[1:2,:,:] = E
F[3,:,:] = B;

# the total electric energy is:

(norm(E)^2 + norm(B)^2)*prod(dx)
#u
#F[3,:,:]



# Cell 33
#surface(y_p,x_p,F[1,:,:])

# Cell 34
Energy_Ks, Energy_Es = get_energy_rel(u,(Box_x, N, J))

# Cell 35
#surface(y_p,x_p,F[3,:,:])

# Cell 38

Coordinate_test(u,Box_x,N)

#println("n_total = $(sum(n .- 1.0))")
#println("v_total = $(sum(u[N+1:2N]))")
#println("E_total = $(sum(E_i))")

#println("S_total = $(sum(S))")

# Cell 40
#includet("aux_functions/aux_functions_RHS.jl")
if nthreads() > 1
    if sbp || holbo || lanczos || no_maxwell
        #TS = zeros(Float64, (2,J...,nthreads()))
        #p_RHS_D = (N, J, Box_x, order, n, S, du, get_density_2D!, get_current_threads_2D!, Interpolate_EBv_1, TS, Dx, Δx, σx, Dy, Δy, σy) ;
        p_RHS_D = (order, N, J, Box_x, order, n, S, du, get_density_2D!, get_current_slim, interp, Dx, Δx, σx, Dy, Δy, σy, no_maxwell, dissipation, factor) ;
        RHS = RHS_D_slim!
        println(RHS)
    elseif weno
        p_RHS_D = (order, N, J, Box_x, dx, order, n, S, du, get_density!, get_current_slim, interp, no_maxwell, parWENOZ, factor)
        RHS = RHS_Flux
    end
else
    if sbp || holbo || lanczos || no_maxwell
        p_RHS_D = (N, J, Box_x, order, n, S, du, get_density_2D!, get_current_rel_2D!, interp, Dx, Δx, σx, Dy, Δy, σy, no_maxwell, dissipation, factor) ;
        RHS = RHS_D
    elseif weno
        p_RHS_D = (order, N, J, Box_x, dx, order, n, S, du, get_density_2D!,get_current_rel_2D!, interp, no_maxwell , parWENOZ, factor)
        RHS = RHS_Flux
    end
end


# Cell 41
RHS

# Cell 42
t = 0.0
j = 1

run_pars = Dict("run_name" => run_name, "par_grid" => par_grid, "par_evolv" => par_evolv, "p_Ini" => (par_f_x, par_f_p))

if full_outputs
    run_name = run_name * "_full"
elseif average_outputs
    run_name = run_name * "_ave"
end

file_name = "Results/"* run_name * ".jld2"
#rm(file_name)




save(file_name, run_pars)

if false # solo para testear
file = jldopen(file_name, "r+")
close(file)
end

if full_outputs
    field_name = "u"
    tiempo = @sprintf("%05d", j)
    jldopen(file_name, "a+") do file
        file[field_name * "/u_$(tiempo)"] = u;
    end
end

if average_outputs
    load_averages_D(file_name, j, par_grid)#, pars_f)
end

if animation
    par = Array{Float64,2}(undef,M_g,length(u));
    par[1,:] = u[:]
end

# Cell 43
@show file_name
@show full_outputs
@show average_outputs


# Cell 44
#includet("aux_functions/aux_functions.jl")
#includet("aux_functions/aux_functions_RHS.jl")
for k in 2:M
    if sbp || holbo || lanczos || no_maxwell
      RK4_Step!(RHS,u,t,dt,p_RHS_D)
    elseif weno
      TVD3_Step!(RHS,u,t,dt,p_RHS_D)
    end
  #global u = [make_periodic!(u[1:4N],Box_x); u[4N+1:end]]
  # make_periodic!(u,Box_x,N) # done inside the RHS_D function at every RK step
  #filter_constant!(u[2N+1:end])
  global t = t + dt
  if (k-1) % (M÷(M_g-1)) == 0
    local j = (k-1)÷(M÷(M_g-1))+1
    make_periodic!(u,Box_x,N) # just to save the correct data
    global Energy_Ks, Energy_Es = get_energy_rel(u,(Box_x, N, J))

    if average_outputs
      load_averages_D(file_name, j, par_grid)#, pars_f)
    end

    if full_outputs
        local tiempo = @sprintf("%05d", j)
        jldopen(file_name, "a+") do file
            file[field_name * "/u_$(tiempo)"] = u
      end
    end

    println("j = $j , t = $t, k = $k, nthreads = $(nthreads()), Total_Energy = $(Energy_Ks + Energy_Es), E_Energy = $(Energy_Es)")

    if animation
      par[j,:] = u[:]
    end
  end
end

# Cell 45
j = 11

Plots.scatter(par[j,1:2D:4N], par[j,2:2D:4N], ts=0.2
, thickness_scaling = 0.3
, markersize = 0.3
, title = "space", legend =:false)

#png("weibel_space_t=10")

# Cell 46
#include("aux_functions/aux_functions.jl")


#make_periodic!(par[M_g,:],Box_x,N)
Coordinate_test(par[M_g,:],Box_x,N)

# Cell 47
j = 11

histogram2d(par[j,1:2D:4N], par[j,2:2D:4N], ts=0.2
, thickness_scaling = 0.3
, markersize = 0.3
, title = "space", legend =:false)

# Cell 48
j = 11
Plots.scatter(par[j,3:2D:4N], par[j,4:2D:4N], ts=0.2
, thickness_scaling = 0.3
, markersize = 0.3
, title = "momentum", legend =:false)

#png("weibel_momentum_t=10")

# Cell 49
j = M_g
Fuj = view(par[j,:],4N+1:4N+3*prod(J))

Fj = reshape(Fuj,3,J...)

#F = view(u[4N+1:end],3,J...)
Ej = Fj[1:2,:,:]
Bj = Fj[3,:,:];

#surface(y_p,x_p,Ej[2,:,:])
#surface(y_p,x_p,Bj[:,:])

#norm(Bj)

# Cell 50
plot(Ej[1,:,50])
plot!(Bj[:,50])

# Cell 51
l = 1
P = 5
plot()
for j in 1:M_g
    Fuj = view(par[j,:],4N+1:4N+3*prod(J))
    Fj = reshape(Fuj,3,J...)
    Ej = Fj[1:2,:,:]
    Bj = Fj[3,:,:];
    plot!(Ej[l,:,P])
    plot!(Bj[:,P])
end
plot!(legend=false)

# Cell 52
norm(par[1,4N+1:4N+3*prod(J)] - par[11,4N+1:4N+3*prod(J)])/norm(par[1,4N+1:4N+3*prod(J)])

# Cell 54


# Cell 56
j = M_g
par_grid = (N, Box_x, J, 5)

#nf = get_density_2D_trans(Val(order), Box_x, par[j,1:4N])

#Ef = Array{Float64,3}(undef,(2,J...))
#ρf = nf.-1.0/prod(J)
#ρf = nf .-1.0
#@show sum(ρf)
#get_E!(Ef,ρf,Box_x);

Fuj = view(par[j,:],4N+1:4N+3*prod(J))
Fj = reshape(Fuj,3,J...)
#F = view(u[4N+1:end],3,J...)
Ej = Fj[1:2,:,:]
@show(norm(Ej))

# Cell 58
nf = get_density_2D_trans(Val(order), Box_x, par[M_g,1:4N])
ρf = nf .-1.0;

# Cell 59
#include("aux_functions/aux_functions_E-B.jl")
x0=[Box_x[1] + (Box_x[2]-Box_x[1])*rand(), Box_x[3] + (Box_x[4]-Box_x[3])*rand()]
@show x0
pars = (x0, 0.2, 6)
constraint_test(Ej, ρf, J, Box_x, ϕ_test, ∇ϕ_test, pars)

# Cell 60
Err = 0.0
Max = 0.0
M = 100
for i in 1:M
    x0=[Box_x[1] + (Box_x[2]-Box_x[1])*rand(), Box_x[3] + (Box_x[4]-Box_x[3])*rand()]
    pars = (x0, 0.2, 6)
    div, rho, err = constraint_test(Ej, ρf, J, Box_x, ϕ_test, ∇ϕ_test, pars)
    Err = Err + abs(err)
    if err > Max
        @show Max = err
        @show x0
    end
    
end
@show Err/M
@show Max;

