using BenchmarkTools
const N = 2000000::Int64
const L = 1.0::Float64
const J = 100::Int64
dx = L/J
const order = 5::Int64
const par_grid =(N, L, J, dx, order)

include("../aux_functions.jl")

global u_r = Float64[]

for i in 1:N
           u1 = L*rand()
           u2 = (1.0 - 2.0*rand())
           global u_r = append!(u_r,[u1,u2])
           #@show u_r
end

u = deepcopy(u_r)

S = Float64[0.0 for i in 1:J]


TS = zeros(Float64,J,nthreads())

par_threads = (par_grid, TS)

@btime get_current_rel_threads!($u, $S, $par_threads);