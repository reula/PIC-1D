#using DifferentialEquations
make_plots = false
if make_plots
using Plots
end
using Statistics
using FFTW
FFTW.set_provider!("mkl")
#import Pkg; Pkg.add("FileIO")
using FileIO
using Base.Threads
using Distributions
#Pkg; Pkg.add("DistributedArrays")
println("nthreads = $(nthreads())")

include("aux_functions.jl")

undamped = false
damped = false
#undamped = true
damped = true

if undamped
    run_name = "landau_undamped_"
    const L = 39.738 
    const J = 3522
elseif damped
    run_name = "landau_damped_"
    const L = 7.455
    const J = 930
end
order = 5
#const N = 80_000
#const N = 800_000
const N = 8_000_000
exp_Theta = 3
exp_t = 1
θ = 10.0^(-exp_Theta)
t = 0.0
t_f = 1.5*10.0^(exp_t)
M = 1501
M_g = 1500 + 1 #number of outputs, starting from the initial data
dt = t_f / (M-1)
t_i = 0.0
#M = convert(Int64,t_f/dt)
#M=1
const κ = 2π/L # for Fourier Transform
dx = L/J
x = [dx*(i-1) for i in 1:J] ;
p = (L, N, J, κ, dx, order)

println("t_f = $(t_f), M = $M, dt = $(dt), exp_Theta = $(exp_Theta)")

animation = false
phase_space_show = false

run_name = run_name * "t$(convert(Int,t_f))_L$(L)_N8_6_J$(J)_M$(M)_o$(order)_T$(exp_Theta)"
println(run_name)

        
#test_parameters(M, M_g, dt, 0.0, t_f)

E = zeros(J)
ϕ = zeros(J)
n = zeros(J) #charge density
S = zeros(J) #carge current
du = zeros(2*N+J); # contains r, v and E


if undamped
    #par_dis = load("Initial_Distributions/par_dis_8_5_abig_undamping.jld2", "par_dis");
    #par_dis = load("Initial_Distributions/par_dis_8_5_undamping.jld2", "par_dis");
elseif damped
    #par_dis = load("Initial_Distributions/par_dis_8_5_abig_damping.jld2", "par_dis");
    #par_dis = load("Initial_Distributions/par_dis_8_5_damping.jld2", "par_dis");
    #par_dis = load("Initial_Distributions/par_dis_8_4_damping.jld2", "par_dis");
    par_dis = load("Initial_Distributions/par_dis_8_6_damping.jld2", "par_dis");
end

#length(par_dis)
@assert length(par_dis) ÷ 2 == N

@time get_density!(par_dis, n, p)
n0 = N/L
#get_ϕ!(ϕ, -n/n0 .+ 1., κ) # chenge the sign here to make it consistent with charge conservation and the time derivative of E
@time get_ϕ!(ϕ, n/n0 .- 1., κ)
@time get_E_from_ϕ!(ϕ,E,dx)
u = [par_dis;E];
length(u)

@time get_current!(u, S, p)

println(maximum(u[1:N]))
println(minimum(u[1:N]))

@time Coordinate_test(u[1:N],L)

println("n_total = $(sum(n .- n0))")
println("v_total = $(sum(u[N+1:2N]))")
println("E_total = $(sum(E))")

println("S_total = $(sum(S))")

if make_plots
    plot(layout=(2,2))
    plot!(subplot=1,x,n/n0, title = "density", legend = :false)
    plot!(subplot=2,x,ϕ, title = "potential", legend = :false)
    plot!(subplot=3,x,E, title = "Electric Field", legend = :false)
    plot!(subplot=4,x,S, title = "Current", legend = :false)

    plot(layout=(2,2))
    histogram!(subplot=1,u[1:N], title = "density", legend = :false, bins = 300)
    histogram!(subplot=2,u[N+1:2*N], title = "velocity", legend = :false)

    histogram!(subplot=3,S, title = "current", legend = :false)

    if phase_space_show
        Plots.scatter(u[1:N],u[N+1:2N]
        , thickness_scaling = 0.3
        , markersize = 0.3
        , title = "phase-space", legend =:false)
    end
end

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



t = 0.0
j = 1



for k in 2:M
  RK4_Step!(RHSC,u,t,dt,p_RHSC)
  global u = [make_periodic!(u[1:N],L); u[N+1:end]]
  #filter_constant!(u[2N+1:end])
  global t = t + dt
  if (k-1) % (M÷(M_g-1)) == 0
    j = (k-1)÷(M÷(M_g-1))+1
    #scatter(plt, u[1:N], u[N+1:2*N])
    Energy_K[j], Energy_E[j] = get_energy(u,(L, N, J))
    E_T[j] = sum(u[2N+1:end])
    v_T[j] = sum(u[N+1:2N])
    get_density!(u, n, p)
    get_current!(u, S, p)
    D_T[j] = sum(n)/n0/J - 1
    S_T[j] = sum(S)/n0/J
    T[j] = var(u[N+1:2N])
    println("j = $j , t = $t, k = $k, nthreads = $(nthreads())")
    if animation
    par[j,:] = u[1:2N]
    end
  end
end

n_F = zeros(J)
S_F = zeros(J)
get_density!(u, n_F, p)
get_current!(u, S_F, p)
run = Dict("run_name" => run_name, "p" => p, "Energy_E" => Energy_E, "Energy_K" => Energy_K, "E_f" => u[2N+1:end], "n_F" => n_F, "S_F" => S_F, "E_T"=> E_T, "v_T" => v_T, "S_T" => S_T, "D_T" => D_T, "T" => T)
save(run_name * "th$(nthreads())_results.jld2", run)

if make_plots
    plot(abs.(Energy_K[2:end] .- Energy_K[1]), title = "Energy conservation", label = "Kinetic Energy"
    #, legend = :outertopright
    , legend = :bottomright)
    plot!(abs.(Energy_E[2:end] .- Energy_E[1]), label = "|Electric Energy|")
    plot!(abs.(Energy_K[2:end]  + Energy_E[2:end] .- (Energy_K[1]+Energy_E[1])) ./ (Energy_K[1]+Energy_E[1]) , yscale=:log10
    #, xscale=:log10
    , label = "Total Energy / Initial Energy -1 ")
    #png(run_name * "energy_conservation")

    abs.(Energy_K[end]  + Energy_E[end] .- (Energy_K[1]+Energy_E[1])) ./ (Energy_K[1]+Energy_E[1])

    #ω = 1
    #tv = [dt*(i-1) for i in 1:M]
    plot(T, label= "θ", title = "Temperature = var", legend = false)
    #plot!(sin.(ω*tv))
    #png(run_name * "temperature")
end

Th = rfft(T)
plot(real.(Th))
argmin(real.(Th))*t_f/length(T) * sqrt(θ)
#length(T)

T[end]

if make_plots
    plot(layout=(2,2))
    plot!(subplot=1,E_T,title="Total Electric Field")
    plot!(subplot=2,v_T./N, title="Total velocity")
    plot!(subplot=3,D_T,title = "Total density")
    plot!(subplot=4,S_T,title = "Total Current")
    #png(run_name * "totals")
end

E_F = zeros(J)
ϕ_F = zeros(J)
n_F = zeros(J)
get_density!(u, n_F, p)
n0 = N/L
get_ϕ!(ϕ_F, n_F/n0 .+ 1, κ)
#get_ϕ!(ϕ_F, n_F .- n0, κ)
get_E_from_ϕ!(ϕ_F,E_F,dx)

if make_plots
    plot(x,u[2N+1:end], label = "dynamical", title="Electric Field with order = $(order)")
    plot!(x,E_F,label="from constraint", ls=:dash)
    #t_f = 40
    #png(run_name * "E")
end

println("averaged total E field = $(sum(E_F))")
println("E_total = $(sum(u[2N+1:end])/J)")
println("Total velocity = $(sum(u[N+1:2N])/N)")
println("Total Charge = $(sum(n_F .- n0))")
println("Final Energy = $(get_energy(u,(L, N, J)))")
get_current!(u, S, p)
println("Total_current = $(sum(S)/J)")

if make_plots
plot(x,u[2N+1:end]-E_F, label = "difference", title="Electric Field")
#png(run_name * "electric_diff")
if make_plots

if phase_space_show
Plots.scatter(u[1:N],u[N+1:2N]
#, thickness_scaling = 0.3
, markersize = 0.3
, title = "phase-space", legend =:false)
#png(run_name * "ps")
end

if animation
    anim = @animate for i = 1:M_g
        Plots.scatter(par[i,1:N], par[i,N+1:2N]
        , markersize = 0.3
        , title = "phase-space"
        , legend=false
        , ylim = (-0.3,0.3)
        , xlim = (0,5)
        )
    end
 
    gif(anim, run_name * "ps_fps5.gif", fps = 5)
end

#load("Results/" * run_name * "results.jld2", "S_F")

#(p, Energy_E, Energy_K, E_f, n_F, S_F, E_T, v_T, S_T, D_T, T) = load("Results/" * run_name * "results.jld2", "p", "Energy_E", "Energy_K", "E_f", "n_F",  "S_F", "E_T", "v_T", "S_T", "D_T", "T");


#using LsqFit


