
"""
This are interpolation functions for getting the Electric field correct.
According the SHARP the second is better. Since it keeps momentum conservation.
Modified so as to use the smallest stencils.
"""
@inline function Interpolate_1_S(order::Int64, vector::Array{Float64,1}, x, J::Int64, L::Float64)
  #stencil = order÷2 
  stencil = Int64(ceil((order+1)/2))
  #stencil = order
  j, y = get_index_and_y(x,J,L)
  vi = 0.0
    for l in (-stencil):(stencil +1)
      vi += vector[mod1(j+l,J)] * Shape(order, -y + l)
    end
  return vi
end

@inline function Interpolate_1_W(order::Int64, vector::Array{Float64,1}, x, J::Int64, L::Float64)
  #stencil = order÷2 
  stencil = Int64(ceil((order+1)/2))
  #stencil = order
  j, y = get_index_and_y(x,J,L)
  vi = 0.0
    for l in (-stencil):(stencil +1)
      vi += vector[mod1(j+l,J)] * W(order, -y + l)
    end
  return vi
end

@inline function Interpolate_2_W(order::Int64, vector, x, J::Int64, L::Float64)
  #stencil = (order+1)÷2
  stencil = Int64(ceil((order+2)/2))
  j, y = get_index_and_y(x,J,L)
  vi = 0.0
    for l in (-stencil):(stencil) 
      vi += (vector[mod1(j+l,J)] + vector[mod1(j+l+1,J)]) * W(order, -y + 1/2 + l) / 2
    end
  return vi
end
"""
Tested, OK
"""
@inline function Interpolate_per_W(order::Int64, vector, x, J::Int64, L::Float64)
  stencil = Int64(ceil((order+1)/2))
  j, y = get_index_and_y(x,J,L)
  vi = 0.0
    for l in (-stencil):-j 
      vi += vector[J+j+l] * W(order, -y + 1/2 + l)
    end
    for l in  max(-stencil,-j+1):min(stencil,J-j)
      vi += vector[j+l] * W(order, -y + 1/2 + l)
    end
    for l in J-j+1:stencil
      vi += vector[j-J+l] * W(order, -y + 1/2 + l)
    end
    # now we increase j by 1
    if j<J
      j += 1
      for l in (-stencil):-j 
        vi += vector[J+j+l] * W(order, -y + 1/2 + l)
      end
      for l in  max(-stencil,-j+1):min(stencil,J-j)
        vi += vector[j+l] * W(order, -y + 1/2 + l)
      end
      for l in J-j+1:stencil
        vi += vector[j-J+l] * W(order, -y + 1/2 + l)
      end
    else  
      j = 1
      for l in (-stencil):-j 
        vi += vector[J+j+l] * W(order, -y + 1/2 + l)
      end
      for l in  max(-stencil,-j+1):min(stencil,J-j)
        vi += vector[j+l] * W(order, -y + 1/2 + l)
      end
      for l in J-j+1:stencil
        vi += vector[j-J+l] * W(order, -y + 1/2 + l)
      end
    end
  return vi / 2
end

# D versions 

"""
Multidimensional version 1-2D
1D checked against the original version
"""
@inline function Interpolate_1_W(order::Int64, vector, x, J::Tuple, Box::Tuple)
  #stencil = order÷2
  stencil = Int64(ceil((order+1)/2))
  D = length(J)
  if D==1
    j, y = get_index_and_y(x,J[1],Box[2])
    vi = 0.0
    for l in (-stencil):(stencil +1)
      vi += vector[mod1(j+l,J[1])] * W(order, -y + l)
    end
    return vi
  elseif D==2
    j = [1,1]
    y = [0.0,0.0]
    @inbounds j, y = get_index_and_y!(j,y,x,J,Box)
    #j, y = get_index_and_y!(x,J,Box)
    vi = similar(vector[1,1])
    vi .= 0.0
    for l in (-stencil):(stencil +1)
      for m in (-stencil):(stencil +1)
        vi += vector[mod1(j[1]+l,J[1]),mod1(j[2]+m,J[2])] * W(order, -y[1] + l) * W(order, -y[2] + m)
      end
    end
    return vi
  else
    error("Not yet implemented for D=$D")
  end
end
"""
Interpolate function for the whole of E + v x B
"""
@inline function Interpolate_EBv_1_S(order::Int64, E::Array{Float64,3}, B::Array{Float64,2}, v::Array{Float64,1}, x, J::NTuple, Box::NTuple)
  #stencil = order÷2
  stencil = Int64(ceil((order+1)/2))
  D = length(J)
  EBv = Array{Float64,1}(undef,2)
  val_order = Val(order)
  if D==1
    j, y = get_index_and_y(x,J[1],Box[2]-Box[1])
    vi = 0.0
    for l in (-stencil):(stencil +1)
      vi += vector[mod1(j+l,J[1])] * Shape(val_order, -y + l)
    end
    return vi
  elseif D==2
    j = [1,1]
    y = [0.0,0.0]

    get_index_and_y!(j,y,x,J,Box)
    #j, y = get_index_and_y!(j,y,x,J,Box)
    #j, y = get_index_and_y!(x,J,Box)
    vi = similar(E[:,1,1])
    vi .= 0.0
    @fastmath for l in (-stencil):(stencil +1)
      s1 = Shape(val_order, -y[1] + l)
      for m in (-stencil):(stencil +1)
        @inbounds EBv[1] = E[1,mod1(j[1]+l,J[1]),mod1(j[2]+m,J[2])] - v[2]*B[mod1(j[1]+l,J[1]),mod1(j[2]+m,J[2])]
        @inbounds EBv[2] = E[2,mod1(j[1]+l,J[1]),mod1(j[2]+m,J[2])] + v[1]*B[mod1(j[1]+l,J[1]),mod1(j[2]+m,J[2])]
        @inbounds vi[:] += EBv * s1 * Shape(val_order, -y[2] + m)
      end
    end
    return vi[:]
  else
    error("Not yet implemented for D=$D")
  end
end

@inline function Interpolate_EBv_1_slim_W(::Val{Order}, E::Array{Float64,3}, B::Array{Float64,2}, v::Array{Float64,1}, j, y , J::NTuple, Box::NTuple) where {Order}
  #stencil = order÷2
  stencil = Int64(ceil((Order+1)/2))
  D = length(J)
  EBv = Array{Float64,1}(undef,2)
  val_order = Val(Order)
  if D==1
    #j, y = get_index_and_y(x,J[1],Box[2]-Box[1])
    vi = 0.0
    for l in (-stencil):(stencil +1)
      @inbounds   vi += vector[mod1(j+l,J[1])] * W(val_order, -y + l)
    end
    return vi
  elseif D==2
    #j = [1,1]
    #y = [0.0,0.0]

    #get_index_and_y!(j,y,x,J,Box)
    #j, y = get_index_and_y!(j,y,x,J,Box)
    #j, y = get_index_and_y!(x,J,Box)
    vi = zeros(2)
    @inbounds for m in (-stencil):(stencil +1)
      w2 = W(val_order, -y[2] + m)
      j2 = mod1(j[2]+m,J[2])
      @inbounds for l in (-stencil):(stencil +1)
        w1 = W(val_order, -y[1] + l)
        j1 = mod1(j[1]+l,J[1])
        ws = w2 * w1
        EBv[1] = E[1,j1,j2] - v[2]*B[j1,j2]
        EBv[2] = E[2,j1,j2] + v[1]*B[j1,j2]
        vi[:] += EBv * ws
      end
    end
    return vi[:]
  else
    error("Not yet implemented for D=$D")
  end
end

function Interpolate_All_EBv_1_slim_W(::Val{Order}, E::Array{Float64,3}, B::Matrix{Float64}, v::Matrix{Float64}, idx, y , J::NTuple, Box::NTuple) where {Order}
  stencil = Int64(ceil((Order+1)/2))
  D = length(J)
  val_order = Val(Order)


  if D==2
    Fi = zeros(Float64, N, 2)
    @threads for i in 1:N
      @inbounds for m in (-stencil):(stencil +1)
                  w2 = W(val_order, -y[i,2] + m)
                  idx2 = mod1(idx[i,2]+m,J[2])
        @inbounds for l in (-stencil):(stencil +1)
                    w1 = W(val_order, -y[i,1] + l)
                    ws = w2 * w1
                    idx1 = mod1(idx[i,1]+l,J[1])

                    EBv1 = E[1,idx1,idx2] - v[i,2]*B[idx1,idx2]
                    EBv2 = E[2,idx1,idx2] + v[i,1]*B[idx1,idx2]
                    Fi[i, 1] += EBv1 * ws
                    Fi[i, 2] += EBv2 * ws
        end
      end
    end
    return Fi
  else
    error("Not yet implemented for D=$D")
  end
end

function Interpolate_All_EBv_1_slim_S(::Val{Order}, E::Array{Float64,3}, B::Matrix{Float64}, v::Matrix{Float64}, idx, y , J::NTuple, Box::NTuple) where {Order}
  stencil = Int64(ceil((Order+1)/2))
  D = length(J)
  val_order = Val(Order)

  if D==2
    Fi = zeros(Float64, N, 2)
    @threads for i in 1:N
      @inbounds for m in (-stencil):(stencil +1)
                  s2 = Shape(val_order, -y[i,2] + m)
                  idx2 = mod1(idx[i,2]+m,J[2])
        @inbounds for l in (-stencil):(stencil +1)
                    s1 = Shape(val_order, -y[i,1] + l)
                    ss = s2 * s1
                    idx1 = mod1(idx[i,1]+l,J[1])

                    EBv1 = E[1,idx1,idx2] - v[i,2]*B[idx1,idx2]
                    EBv2 = E[2,idx1,idx2] + v[i,1]*B[idx1,idx2]
                    Fi[i, 1] += EBv1 * ss
                    Fi[i, 2] += EBv2 * ss
        end
      end
    end
    return Fi
  else
    error("Not yet implemented for D=$D")
  end
end

@inline function Interpolate_All_EBv_2_slim(::Val{Order}, E::Array{Float64,3}, B::Matrix{Float64}, v::Matrix{Float64}, idx, y , J::NTuple, Box::NTuple) where {Order}
  stencil = Int64(ceil((Order+1)/2))
  D = length(J)
  val_order = Val(Order)
  #dV = prod(differentials(Box,J))

  if D==2
    Fi = zeros(Float64, N, 2)
    @threads for i in 1:N
      @inbounds for m in (-stencil):(stencil +1)
        w2 = W(val_order, -y[i,2] + 1/2 + m)/4
        idx2 = mod1(idx[i,2]+m,J[2])
        idx2p1 = mod1(idx[i,2]+m+1,J[2])
        @inbounds for l in (-stencil):(stencil +1)
          w1 = W(val_order, -y[i,1] + 1/2 + l)
          ws = w2 * w1
          idx1 = mod1(idx[i,1]+l,J[1])
          idx1p1 = mod1(idx[i,1]+l+1,J[1])
          EBv1 = E[1,idx1,idx2] + E[1,idx1p1,idx2] + E[1,idx1,idx2p1] + E[1,idx1p1,idx2p1]
          EBv2 = E[2,idx1,idx2] + E[2,idx1p1,idx2] + E[2,idx1,idx2p1] + E[2,idx1p1,idx2p1]
          Bz = B[idx1,idx2] + B[idx1p1,idx2] + B[idx1,idx2p1] + B[idx1p1,idx2p1]
          EBv1 = EBv1 - v[i,2]*Bz
          EBv2 = EBv2 + v[i,1]*Bz
          Fi[i, 1] += EBv1 * ws
          Fi[i, 2] += EBv2 * ws
        end
      end
    end
    return Fi
  else
    error("Not yet implemented for D=$D")
  end
end