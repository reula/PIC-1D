using BenchmarkTools

const N = 200_000::Int64
const J = (100::Int64,200::Int64)
const Box = (0.0::Float64,1.0::Float64,0.0::Float64,1.0::Float64)
const corder = Val(5)
include("../aux_functions/aux_functions.jl")

u_r = Float64[]
for i in 1:N
    u1 = Box[2] * rand()
    u2 = Box[4] * rand()
    u3 = (1.0 - 2.0 * rand())
    u4 = (1.0 - 2.0 * rand())
    append!(u_r, [u1, u2, u3, u4])
    #@show u_r
end

u = deepcopy(u_r)
@btime get_current_2D_trans(corder, $N, $J, $Box, $u)