
"""
This are interpolation functions for getting the Electric field correct.
According the SHARP the second is better. Since it keeps momentum conservation.
Modified so as to use the smallest stencils.
"""
@inline function Interpolate_1(order, vector, x, J, L)
  
  #stencil = order÷2 
  stencil = Int64(ceil((order+1)/2))
  j, y = get_index_and_y(x,J,L)
  vi = 0.0
    for l in (-stencil):(stencil +1)
      vi += vector[mod1(j+l,J)] * W(order, -y + l)
    end
  return vi
end

@inline function Interpolate_2(order, vector, x, J, L)
  #stencil = (order+1)÷2
  stencil = Int64(ceil((order+1)/2))
  j, y = get_index_and_y(x,J,L)
  vi = 0.0
    for l in (-stencil):(stencil) 
      vi += (vector[mod1(j+l,J)] + vector[mod1(j+l+1,J)]) * W(order, -y + 1/2 + l) / 2
    end
  return vi
end
"""
Terminado, pero sin probar.
"""
@inline function Interpolate_per(order, vector, x, J, L)
  j, y = get_index_and_y(x,J,L)
  vi = 0.0
    for l in (-order):-j 
      vi += vector[J+j+l] * W(order, -y + 1/2 + l)
    end
    for l in  max(-order,-j+1):min(order,J-j)
      vi += vector[j+l] * W(order, -y + 1/2 + l)
    end
    for l in J-j+1:order
      vi += vector[j+l] * W(order, -y + 1/2 + l)
    end
    # now we increase j by 1
    j += 1
    for l in (-order):-j 
      vi += vector[J+j+l] * W(order, -y + 1/2 + l)
    end
    for l in  max(-order,-j+1):min(order,J-j)
      vi += vector[j+l] * W(order, -y + 1/2 + l)
    end
    for l in J-j+1:order
      vi += vector[j+l] * W(order, -y + 1/2 + l)
    end
  return vi / 2
end

# 2D versions
@inline function Interpolate_1(order::Int64, vector::Array{Float64}, x, J::Tuple, Box::Tuple)
  #stencil = order÷2
  stencil = Int64(ceil((order+1)/2))
  j, y = get_index_and_y!(x,J,Box)
  vi = 0.0
    for l in (-stencil):(stencil +1)
      vi += vector[mod1(j+l,J)] * W(order, -y + l)
    end
  return vi
end

