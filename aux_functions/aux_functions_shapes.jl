"""
Derivatives of shape functions.
They have support for -(order+1)/2 =< y =< (order+1)/2
"""
@inline function W(order::Int,y::Float64)
  #y = norm(y)
  y = abs(y)
  if order == 0
    return  (y <= 1/2) ? 1 : 0
  elseif order ==1
    return  (y <= 1) ? 1 - y : 0
  elseif order == 2
    #return (y <= 1/2) ? 3/4 - y^2  : (((y > 1/2) && (y <= 3/2)) ? (3 - 2*y)^2 / 8 : 0)
    return (y <= 1/2) ? 3/4 - y^2  : (((y <= 3/2)) ? (3 - 2*y)^2 / 8 : 0)
  elseif order == 3
    #return (y <= 1) ? 2/3 - y^2 + y^3 / 2 : (((y > 1) && (y <= 2)) ? (2 - y)^3 / 6 : 0)
    return (y <= 1) ? 2/3 - y^2 + y^3 / 2 : ((y <= 2) ? (2 - y)^3 / 6 : 0)
  elseif order == 4
    #return (y <= 1/2) ? 115/192 - 5y^2/8 + y^4/4 : (((y > 1/2) && (y <= 3/2)) ? (55 + 20y -120y^2 + 80y^3 - 16y^4)/96 : (((y > 3/2) && (y < 5/2)) ? (5 - 2y)^4/384 : 0))
    return (y <= 1/2) ? 115/192 - 5y^2/8 + y^4/4 : (((y <= 3/2)) ? (55 + 20y -120y^2 + 80y^3 - 16y^4)/96 : (((y < 5/2)) ? (5 - 2y)^4/384 : 0))
  elseif order == 5
    #return (y <= 1) ? 11/20 - y^2/2 + y^4/4 - y^5/12 : (((y > 1) && (y <= 2)) ? 17/40 + 5y/8 - 7y^2/4 + 5y^3/4 - 3y^4/8 + y^5/24 : (((y > 2) && (y < 3)) ? (3 - y)^5/120 : 0))
    return (y <= 1) ? 11/20 - y^2/2 + y^4/4 - y^5/12 : (((y <= 2)) ? 17/40 + 5y/8 - 7y^2/4 + 5y^3/4 - 3y^4/8 + y^5/24 : (((y < 3)) ? (3 - y)^5/120 : 0))
  
  else
    error("order = $order not yet implemented ")
  end
end


"""
Alternative definition
"""
function W_alt(order::Int,y::Float64)
  y = abs(y)
  if order == 0
    return  (y > 1/2) ? 0 : 1
  elseif order ==1
    return  (y > 1) ? 0 : 1 - y 
  elseif order == 2
    #return (y <= 1/2) ? 3/4 - y^2  : (((y > 1/2) && (y <= 3/2)) ? (3 - 2*y)^2 / 8 : 0)
    return (y > 3/2) ? 0 : (((y > 1/2) && (y <= 3/2)) ? (3 - 2*y)^2 / 8 : 3/4 - y^2)
  elseif order == 3
    #return (y <= 1) ? 2/3 - y^2 + y^3 / 2 : (((y > 1) && (y <= 2)) ? (2 - y)^3 / 6 : 0)
    return (y > 2) ? 0 : (((y > 1) && (y <= 2)) ? (2 - y)^3 / 6 : 2/3 - y^2 + y^3 / 2)
  elseif order == 4
    #return (y <= 1/2) ? 115/192 - 5y^2/8 + y^4/4 : (((y > 1/2) && (y <= 3/2)) ? (55 + 20y -120y^2 + 80y^3 - 16y^4)/96 : (((y > 3/2) && (y < 5/2)) ? (5 - 2y)^4/384 : 0))
    return (y > 5/2) ? 0 : (((y > 3/2) && (y < 5/2)) ? (5 - 2y)^4/384 : (((y > 1/2) && (y <= 3/2)) ? (55 + 20y -120y^2 + 80y^3 - 16y^4)/96 : 115/192 - 5y^2/8 + y^4/4))
  elseif order == 5
    #return (y <= 1) ? 11/20 - y^2/2 + y^4/4 - y^5/12 : (((y > 1) && (y <= 2)) ? 17/40 + 5y/8 - 7y^2/4 + 5y^3/4 - 3y^4/8 + y^5/24 : (((y > 2) && (y < 3)) ? (3 - y)^5/120 : 0))
    return (y > 3) ? 0 : (((y > 2) && (y < 3)) ? (3 - y)^5/120 : (((y > 1) && (y <= 2)) ? 17/40 + 5y/8 - 7y^2/4 + 5y^3/4 - 3y^4/8 + y^5/24 : 11/20 - y^2/2 + y^4/4 - y^5/12))
  else
    error("order = $order not yet implemented ")
  end
end

"""
Shape functions are W functions of an order less. 
The dx in the Appendix is added when used so as not to carry dx all over. 
So this definition DIFFERS FROM THE PAPER BY A DX!
They have support on -order/2 =< y =< order/2
The orders goes from 1 to 6 (order zero is not defined)
While the orders of W goes from 0 to 5
"""
@inline function Shape(order::Int,y::Float64) 
  return W(order - 1,y)
end
