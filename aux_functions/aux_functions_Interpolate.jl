
"""
This are interpolation functions for getting the Electric field correct.
According the SHARP the second is better. Since it keeps momentum conservation.
Modified so as to use the smallest stencils.
"""
@inline function Interpolate_1(order::Int64, vector::Array{Float64,1}, x, J::Int64, L::Float64)
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

@inline function Interpolate_2(order::Int64, vector, x, J::Int64, L::Float64)
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
@inline function Interpolate_per(order::Int64, vector, x, J::Int64, L::Float64)
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
@inline function Interpolate_1(order::Int64, vector, x, J::Tuple, Box::Tuple)
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
@inline function Interpolate_EBv_1(order::Int64, E::Array{Float64,3}, B::Array{Float64,2}, v::Array{Float64,1}, x, J::NTuple, Box::NTuple)
  #stencil = order÷2
  stencil = Int64(ceil((order+1)/2))
  D = length(J)
  EBv = Array{Float64,1}(undef,2)
  if D==1
    j, y = get_index_and_y(x,J[1],Box[2]-Box[1])
    vi = 0.0
    for l in (-stencil):(stencil +1)
      vi += vector[mod1(j+l,J[1])] * W(order, -y + l)
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
    for l in (-stencil):(stencil +1)
      for m in (-stencil):(stencil +1)
        EBv[1] = E[1,mod1(j[1]+l,J[1]),mod1(j[2]+m,J[2])] - v[2]*B[mod1(j[1]+l,J[1]),mod1(j[2]+m,J[2])]
        EBv[2] = E[2,mod1(j[1]+l,J[1]),mod1(j[2]+m,J[2])] + v[1]*B[mod1(j[1]+l,J[1]),mod1(j[2]+m,J[2])]
        vi += EBv * W(order, -y[1] + l) * W(order, -y[2] + m)
      end
    end
    return vi[:]
  else
    error("Not yet implemented for D=$D")
  end
end

