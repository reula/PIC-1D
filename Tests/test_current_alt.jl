using BenchmarkTools
const N = 2_000_000::Int64
const J = (100::Int64,200::Int64)
const Box = (0.0::Float64,1.0::Float64,0.0::Float64,1.0::Float64)
const order = 5::Int64
const par_grid =(N, Box, J, order)
include("../aux_functions.jl")

global u_r = Float64[]

for i in 1:N
           u1 = Box[2]*rand()
           u2 = Box[4]*rand()
           u3 = (1.0 - 2.0*rand())
           u4 = (1.0 - 2.0*rand())
           global u_r = append!(u_r,[u1,u2,u3,u4])
           #@show u_r
end

u = deepcopy(u_r)
S1 = [0.0::Float64 for l in 1:2, i in 1:J[1], j in 1:J[2]]
TS = zeros(Float64,(nthreads(),2,J...))

par_threads = (par_grid, TS)

@btime get_current_threads_2D_alt!($u, $S1, $par_threads)


