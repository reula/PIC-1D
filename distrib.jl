using Plots
#using Pkg; Pkg.add("DistributedArrays")
#using Transducers
using FileIO
using Distributed
using DistributedArrays
using DistributedArrays.SPMD
#@everywhere using SharedArrays
@everywhere using Distributed
@everywhere using DistributedArrays
@everywhere using DistributedArrays.SPMD
using BenchmarkTools

#rmprocs(2:37)
d_closeall()

N_workers = 1
addprocs(N_workers)

println(workers())

@everywhere include_string(Main, $(read("aux_functions.jl", String)), "aux_functions.jl")

order = 1
const L = 5
const N = 80000
#const N = 20000
const J = 50
const κ = 2π/L # for Fourier Transform
dx = L/J
x = [dx*(i-1) for i in 1:J] ;
p = (L, N, J, κ, dx, order)

DS = dzeros((J,nworkers()), workers(), [1,nworkers()]);

pars = zeros(2N)
#pars = dones((N,), workers(), [nworkers()] );
par_dis = load("Initial_Distributions/par_dis_L5_N8_5_theta001.jld2", "par_dis");
reorder_particles!(par_dis,pars)

Dpars = distribute(pars);

if false
    @benchmark begin
        Sf = zeros(J)
    spmd(get_current_ro_par, Dpars, DS, p)
    for j in 1:J
        for i in 1:nworkers()
            Sf[j] += DS[j,i]
        end
    end
    end
end

if true
    @benchmark begin
    S = zeros(J)
    get_current_ro!(pars, S, p)
    end
end


#(S - Sf)'*(S-Sf)

