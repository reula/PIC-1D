using BenchmarkTools
const N = 2_000::Int64
const J = (100::Int64,200::Int64)
const Box = (0.0::Float64,1.0::Float64,0.0::Float64,1.0::Float64)
const order = 5::Int64
const par_grid =(N, Box, J, order)
include("../aux_functions.jl")

function smvectors(x::Matrix{T}, ::Val{N}) where {T,N}
    size(x,1) == N || error("sizes mismatch")
    isbitstype(T) || error("use for bitstypes only")
    reinterpret(MVector{N,T}, vec(x))
end

println("starting to compute")

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

us = MVector{4N}(u)

S1 = [[0.0::Float64,0.0::Float64] for i in 1:J[1], j in 1:J[2]]
Ss = MVector{2*J[1]*J[2]}
TS = zeros(Float64,(2,J...,nthreads()))
TSs = MVector{2*J[1]*J[2]*nthreads()}

par_threads = (par_grid, TSs)

@btime get_current_threads_2D!($us, $Ss, $par_threads)


