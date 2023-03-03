"""
The following function evaluates the electric field on a uniform grid from the electric potential.

    // Calculate electric field from potential
"""
function get_E_from_ϕ!(ϕ, E, dx)
      J = length(E)
      for j in 2:J-1
        E[j] = (ϕ[j-1] - ϕ[j+1]) / 2. / dx
      end
      E[1] = (ϕ[J] - ϕ[2]) / 2. / dx;
      E[J] = (ϕ[J-1] - ϕ[1]) / 2. / dx;
end

function get_E_from_ϕ!(ϕ, E, Box_x::Tuple)
    J = size(E[1])
    D = length(J)
    dx = differentials(Box_x,J)
    for i in 1:D
      for k in 1:J[1]
        for j in 1:J[2]
         # E[j,k][i] = (ϕ[j-1] - ϕ[j+1]) / 2. / dx
      end
    end
  end

# UNFINISHED
end

""" The following routine solves Poisson's equation in 1-D to find the instantaneous electric potential on a uniform grid.

// Solves 1-d Poisson equation:
//    d^u / dx^2 = v   for  0 <= x <= L
// Periodic boundary conditions:
//    u(x + L) = u(x),  v(x + L) = v(x)
// Arrays u and v assumed to be of length J.
// Now, jth grid point corresponds to
//    x_j = j dx  for j = 0,J-1
// where dx = L / J. L / (J-1) in Julia
// Also,
//    kappa = 2 pi / L
"""
function get_ϕ!(ϕ, ρ, κ)
  #V = fill(0.0+im*0.,J) 
  #U = fill(0.0+im*0.,J÷2+1) 
  J = length(ρ)
  # Fourier transform source term
  V = rfft(ρ)

  # Calculate Fourier transform of u

  V[1] =  0.
  #for j in  2:(J÷2+1)
  for j in 2:size(V,1) 
    V[j] = - V[j] / (j-1)^2 / κ^2
  end

  # Inverse Fourier transform to obtain u
  ϕ[:] = irfft(V,J)
end

"""
Solves the Poisson equation in D = 1,2, and 3 dimensions with homogeneous boundary conditions
in a rectangular grid.
The solution is writen in the matrix ϕ, and the data in the matrix ρ, 
Box is the square box where the solution is looked at.
Checked in poisson.ipynb for D = 1,2,3
"""
function get_ϕ_D!(ϕ, ρ, Box)
  #V = fill(0.0+im*0.,J) 
  #U = fill(0.0+im*0.,J÷2+1) 
  J = size(ρ)
  D = length(J)
  
  Box_array = [i for i in Box]

  κ = 2π./(Box_array[2:2:end] - Box_array[1:2:end-1])

  # Fourier transform source term
  V = rfft(ρ)

  # Calculate Fourier transform of u

    k1 = rfftfreq(J[1]).*κ[1]*J[1]

    if D==1
      for i in 2:size(V,1)
          V[i] = - V[i] / k1[i]^2
      end
      V[1] = 0.0
      ϕ[:] = irfft(V,J[1])
    elseif D==2
      k2 = fftfreq(J[2]).*κ[2]*J[2]
      for i in 1:size(V,1)
        for j in 1:size(V,2)
          V[i,j] = - V[i,j] / (k1[i]^2 + k2[j]^2 + eps(1.0))
        end
      end
    V[1,1] =  0.0
    # Inverse Fourier transform to obtain u
    ϕ[:,:] = irfft(V,J[1])
    elseif D==3
      k2 = fftfreq(J[2]).*κ[2]*J[2]
      k3 = fftfreq(J[3]).*κ[3]*J[3]
      for i in 1:size(V,1)
        for j in 1:size(V,2)
          for k in 1:size(V,3)
          V[i,j,k] = - V[i,j,k] / (k1[i]^2 + k2[j]^2 + k3[k]^2 + eps(1.0))
          end
        end
      end
    V[1,1,1] = 0.0
    ϕ[:,:,:] = irfft(V,J[1])
    else
      error("not implemented for D=$D")
    end
  
  
end

"""
get_E_direct!(E,ρ_m)
computes the Electric field from the constraint equation div E = ρ 
using E_k = E_{k-1} + \rho_{k-1/2}*dx
then substracts a constant field so that E_T = 0
"""
function get_E_direct!(E,ρ_m,par_grid)
  N, L, J, dx, order = par_grid
  E[1] = 0.0
  E_T = 0.0
  for j in 2:J
    E[j] = E[j-1] + ρ_m[j-1]*dx
    E_T = E_T + E[j]
  end
  return E .= E .- E_T/J
end
"""
Takes out the mass of the grid function so that the sum is now null
"""
function filter_constant!(E)
  J = length(E)
  V = rfft(E)
  V[1] = 0.0 #extract the first component
  E[:] = irfft(V,J)
end

"""
Computes directly the Electric field in a Box with the potential with homogeneous Dirichlet conditions.
Cheked for 1, 2 and 3 dimensions in poisson.ipynb
"""
function get_E!(E, ρ ,Box)
  J = size(ρ)
  D = length(J)
  
  Box_array = [i for i in Box]

  κ = 2π./(Box_array[2:2:end] - Box_array[1:2:end-1])

  # Fourier transform source term
  V = rfft(ρ)

  # Calculate Fourier transform of u

    k1 = rfftfreq(J[1]).*κ[1]*J[1]

    if D==1
      Ek = Array{ComplexF64,D}(undef,size(V))
      for i in 2:size(V,1)
          Ek[i] = im * V[i] / k1[i]
      end
      Ek[1] = 0.0       
      E[:] = irfft(Ek,J[1])

    elseif D==2
      #@show J, κ
      k2 = fftfreq(J[2]).*κ[2]*J[2]
      Ek1 = Array{ComplexF64,D}(undef,(size(V)))
      Ek2 = Array{ComplexF64,D}(undef,(size(V)))
      for i in 1:size(V,1)
        for j in 1:size(V,2)
          Ek1[i,j] = V[i,j] / (k1[i]^2 + k2[j]^2 + eps(1.0)) * im *k1[i]
          Ek2[i,j] = V[i,j] / (k1[i]^2 + k2[j]^2 + eps(1.0)) * im *k2[j]
        end
      end
      Ek1[1,1] =  0.0 + im*0.0;
      Ek2[1,1] =  0.0 + im*0.0; 
      # Inverse Fourier transform to obtain u
      E[1,:,:] = irfft(Ek1,J[1])
      E[2,:,:] = irfft(Ek2,J[1])
    elseif D==3
      k2 = fftfreq(J[2]).*κ[2]*J[2]
      k3 = fftfreq(J[3]).*κ[3]*J[3]
      Ek1 = Array{ComplexF64,D}(undef,(size(V)))
      Ek2 = Array{ComplexF64,D}(undef,(size(V)))
      Ek3 = Array{ComplexF64,D}(undef,(size(V)))
      for i in 1:size(V,1)
        for j in 1:size(V,2)
          for k in 1:size(V,3)
          Ek1[i,j,k] = V[i,j,k] / (k1[i]^2 + k2[j]^2 + k3[k]^2 + eps(1.0)) * im * k1[i]
          Ek2[i,j,k] = V[i,j,k] / (k1[i]^2 + k2[j]^2 + k3[k]^2 + eps(1.0)) * im * k2[j]
          Ek3[i,j,k] = V[i,j,k] / (k1[i]^2 + k2[j]^2 + k3[k]^2 + eps(1.0)) * im * k3[k]
          end
        end
      end
      Ek1[1,1,1] =  0.0 + im*0.0; Ek2[1,1,1] =  0.0 + im*0.0; Ek3[1,1,1] =  0.0 + im*0.0; 
      E[1,:,:,:] = irfft(Ek1,J[1])
      E[2,:,:,:] = irfft(Ek2,J[1])
      E[3,:,:,:] = irfft(Ek3,J[1])
    else
      error("not implemented for D=$D")
    end
  
  

  
end
