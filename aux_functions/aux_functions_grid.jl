############## GRID FUNCTIONS #################
"""
Grid volume
"""
function volume_old(Box) 
  if length(Box) == 2
    return abs(Box[2]-Box[1])
  elseif length(Box) == 4
    return abs((Box[2]-Box[1])*(Box[4]-Box[3]))
  elseif length(Box) == 6
    return abs((Box[2]-Box[1])*(Box[4]-Box[3])*(Box[6]-Box[5]))
  else 
    error("wrong length for Box")
  end
end

function volume(Box) 
  Box_array = [i for i in Box]
  L = Box_array[2:2:end] - Box_array[1:2:end-1]
  return prod(L)
end

"""
Grid Differential 
"""
function differentials(Box,J,periodic=true)
  D = length(J)
  dd = zeros(D)
  if periodic
    for i in 1:D
      dd[i] = (Box[2i] - Box[2i-1])/J[i]
    end
  else
    for i in 1:D
      dd[i] = (Box[2i] - Box[2i-1])/(J[i]-1)
    end
  end
  return dd[:]
end

"""
Volume differential
"""
function vol_diff(Box,J,periodic=true)
  prod(differentials(Box, J))
end


function get_index_and_distance(s,dx,L)
    #if s < 0
    #    s = s + L
    #end
    #if s > L
    #    s = s - L 
    #end
    s = mod1(s,L)
    j = floor(Int64, s ÷ dx) + 1 #find the grid space where it is.
    if j > J || j < 1
      error("j = $j")
    end
    y = (s % dx)/dx #how far is there 
    return j, y
end

"""
given a number s in between x_j and x_{j+1} computes y = (s - x_j)/dx and return j and y.
get_index_and_y(0.4,2,1)
dx = 0.5, j = 1, y = 0.4/0.5 = 0.8
get_index_and_y(0.7,2,1) = j = 2, 0.2/0.5= 0.4
"""
@inline function get_index_and_y(s, J, L, yshift=0.0)
  y, j = modf((s / L * J + J) % J)
  floor(Int64, j) + 1, y - yshift
end

"""
same as before, but for an arbitrary box and dimension.
Jt = (100, 200, 200)
Box = (0.0, 10, 0.0, 20, -20.0, 0.0)
ss = [8.74, 8.74, -11.26]
j=[1,1,1]
y=zeros(3)
get_index_and_y!(j,y,ss,Jt,Box)
([88, 88, 88], [0.4000000000000057, 0.39999999999997726, 0.4000000000000057])
"""
@inline function get_index_and_y(r, J::Tuple, Box::Tuple, yshift=0.0)
  j = zeros(Int64, length(J))
  y = zeros(length(J))
  for i in 1:length(J)
    idx, y = get_index_and_y(r[i], J[i], Box[2i] - Box[2i-1], yshift)

    # tmp = (r[i]/(Box[2i] - Box[2i-1])*J[i] + J[i])%J[i]
    # int = floor(Int,tmp)
    # j[i], y[i] = int + 1, tmp - int

    # y[i], int = modf((r[i]/(Box[2i] - Box[2i-1])*J[i] + J[i])%J[i])
    # j[i] = convert(Int, int) + 1
  end
  j, y
end

@inline function get_indices_and_y_trans!(idx, y, r, sz, L, yshift=0.0)
  N = size(r)[1]
  D = size(r)[2]
  for d = 1:D
    @threads for i = 1:N
      @inbounds tmp = (r[i, d] / L[d] * sz[d] + sz[d]) % sz[d]
      idx[i, d] = floor(Int64, tmp) + 1
      y[i, d] = (tmp % 1) - yshift
    end
  end
end

@inline function get_indices_and_y_trans(r, sz, L, yshift=0.0)
  N = size(r)[1]
  D = size(r)[2]
  idx = Matrix{Int64}(undef, N, D)
  y = Matrix{Float64}(undef, N, D)
  get_indices_and_y_trans!(idx, y, r, sz, L, yshift)
  idx, y
end

"""
same as before, but for an arbitrary box and dimension.
Jt = (100, 200, 200)
Box = (0.0, 10, 0.0, 20, -20.0, 0.0)
ss = [8.74, 8.74, -11.26]
j=[1,1,1]
y=zeros(3)
get_index_and_y!(j,y,ss,Jt,Box)
([88, 88, 88], [0.4000000000000057, 0.39999999999997726, 0.4000000000000057])
"""
@inline function get_index_and_y_old!(j::AbstractArray{Int64,1}, y::AbstractArray{Float64,1}, r, J::Tuple,Box::Tuple) 
  for i in 1:length(J)
    # j[i], y[i] = get_index_and_y(r[i], J[i],Box[2i] - Box[2i-1])
    y[i] =  (r[i]/(Box[2i] - Box[2i-1])*J[i] + J[i])%J[i]
    j[i] = floor(Int,y[i]) + 1
    y[i] = (y[i]%1)
  end
end

@inline function get_index_and_y!(j::AbstractArray{Int64,1}, y::AbstractArray{Float64,1}, r, J::Tuple, Box::Tuple, yshift)
  for i in eachindex(J)
    j[i], y[i] = get_index_and_y(r[i], J[i], Box[2i] - Box[2i-1], yshift)
    # y[i] =  (r[i]/(Box[2i] - Box[2i-1])*J[i] + J[i])%J[i]
    # j[i] = floor(Int,y[i]) + 1
    # y[i] = (y[i]%1)
  end
end

function get_index_and_y_alt!(j::Array{Int64,1}, y::Array{Float64,1}, r, J::Tuple,Box::Tuple)
  J_ar = [i for i in J]
  Box_ar = [i for i in Box] 
    y .=  (r./(Box_ar[2:2:end] - Box_ar[1:2:end-1]).*J_ar .+ J_ar).%J_ar
    j .= floor.(Int,y) .+ 1
    y .= mod1.(y,1)
  #return j[:], y[:]
end