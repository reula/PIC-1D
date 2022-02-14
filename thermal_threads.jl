using Statistics
using FFTW
FFTW.set_provider!("mkl")
import Pkg; Pkg.add("FileIO")
using FileIO
using Base.Threads

include("aux_functions.jl")

run_name = "long_"
order = 1
const L = 5
#N = 80000
const N = 20000
const J = 50
exp_Theta = 1
exp_t = 1 # 5
θ = 10.0^(-exp_Theta)
t = 0.0
t_f = 10.0^(exp_t)
M = 1_000_001
M_g = 1000 + 1 #number of outputs, starting from the initial data
dt = t_f / (M-1)
t_i = 0.0
#M = convert(Int64,t_f/dt)
#M=1
const κ = 2π/L # for Fourier Transform
dx = L/J
x = [dx*(i-1) for i in 1:J] ;
p = (L, N, J, κ, dx, order)

animation = false

run_name = run_name * "t$(convert(Int,t_f))_L$(L)_N2_5_J$(J)_M$(M)_o$(order)_T$(exp_Theta)"

E = zeros(J)
ϕ = zeros(J)
n = zeros(J) #charge density
S = zeros(J) #carge current
du = zeros(2*N+J); # contains r, v and E

par_dis = load("Initial_Distributions/par_dis_L5_N2_5_theta001.jld2", "par_dis");

@assert length(par_dis) ÷ 2 == N

@time get_density!(par_dis, n, p)
n0 = N/L
#get_ϕ!(ϕ, -n/n0 .+ 1., κ) # chenge the sign here to make it consistent with charge conservation and the time derivative of E
@time get_ϕ!(ϕ, n/n0 .- 1., κ)
@time get_E_from_ϕ!(ϕ,E,dx)
global u = [par_dis;E];
length(u)

@time get_current!(u, S, p)

@time Coordinate_test(u[1:N],L)

println("n_total = $(sum(n .- n0))")
println("v_total = $(sum(u[N+1:2N]))")
println("E_total = $(sum(E))")

println("S_total = $(sum(S))")

if nthreads() > 1
    TS = zeros(J, nthreads())
    p_RHSC = (N, J, L, dx, order, n, S, du, get_density!, get_current_threads!, Interpolate_2, TS) ;
else
    p_RHSC = (N, J, L, dx, order, n, S, du, get_density!, get_current!, Interpolate_2) ;
end
    
Energy_K = zeros(M_g)
Energy_E = zeros(M_g)
E_T = zeros(M_g)
v_T = zeros(M_g)
D_T = zeros(M_g)
S_T = zeros(M_g)
T = zeros(M_g)
if animation
    par = zeros(M_g,2N)
end
Energy_K[1], Energy_E[1]  = get_energy(u,(L, N, J))
E_T[1] = sum(u[2N+1:end])
v_T[1] = sum(u[N+1:2N])
get_density!(u, n, p)
get_current!(u, S, p)
D_T[1] = sum(n)/n0/J - 1
S_T[1] = sum(S)/n0/J
T[1] = var(u[N+1:2N])
if animation
    par[1,:] = u[1:2N]
end
    
    
    
global t = 0.0
global j = 1
    

    for k in 2:(M+1)
        RK4_Step!(RHSC,u,t,dt,p_RHSC)
        global u = [make_periodic!(u[1:N],L); u[N+1:end]]
        #filter_constant!(u[2N+1:end])
        global t = t + dt
        if (k-1) % (M÷(M_g-1)) == 0
          global j = j+1
          #scatter(plt, u[1:N], u[N+1:2*N])
          Energy_K[j], Energy_E[j] = get_energy(u,(L, N, J))
          E_T[j] = sum(u[2N+1:end])
          v_T[j] = sum(u[N+1:2N])
          get_density!(u, n, p)
          get_current!(u, S, p)
          D_T[j] = sum(n)/n0/J - 1
          S_T[j] = sum(S)/n0/J
          T[j] = var(u[N+1:2N])
          println("j = $j , t = $t, k = $k")
          if animation
          par[j,:] = u[1:2N]
          end
        end
    end

n_F = zeros(J)
S_F = zeros(J)
get_density!(u, n_F, p)
get_current!(u, S_F, p)
save(run_name * "results.jld2", Dict("p" => p, "Energy_E" => Energy_E, "Energy_K" => Energy_K, "E_f" => u[2N+1:end], "n_F" => n_F, "S_F" => S_F, "E_T"=> E_T, "v_T" => v_T, "S_T" => S_T, "D_T" => D_T, "T" => T))

#plot(layout=(2,2))
#plot!(subplot=1,E_T,title="Total Electric Field")
#plot!(subplot=2,v_T./N, title="Total velocity")
#plot!(subplot=3,D_T,title = "Total density")
#plot!(subplot=4,S_T,title = "Total Current")
#png("$run_name * totals")

E_F = zeros(J)
ϕ_F = zeros(J)
n_F = zeros(J)
get_density!(u, n_F, p)
n0 = N/L
get_ϕ!(ϕ_F, n_F/n0 .+ 1, κ)
#get_ϕ!(ϕ_F, n_F .- n0, κ)
get_E_from_ϕ!(ϕ_F,E_F,dx)
#plot(x,u[2N+1:end], label = "dynamical", title="Electric Field with order = $(order)")
#plot!(x,E_F,label="from constraint", ls=:dash)
#t_f = 40
#png("$run_name * E")

println("averaged total E field = $(sum(E_F))")
println("E_total = $(sum(u[2N+1:end])/J)")
println("Total velocity = $(sum(u[N+1:2N])/N)")
println("Total Charge = $(sum(n_F .- n0))")
println("Final Energy = $(get_energy(u,(L, N, J)))")
get_current!(u, S, p)
println("Total_current = $(sum(S)/J)")

#plot(x,u[2N+1:end]-E_F, label = "difference", title="Electric Field")
#png("$run_name * electric_diff")

#scatter(u[1:N],u[N+1:2N]
#, markersize = 0.3
#, title = "phase-space", legend =:false)
#png("run_name * ps")

#load(run_name * "results.jld2", "S_F")

