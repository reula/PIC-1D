#=
################################################################################################
Some functions to handle data
################################################################################################
=#
function get_averages(v,par_grid,par_evolv, par_f)
    (N, L, J, dx, order) = par_grid
    (t_i, t_f, M, M_g, dt) = par_evolv
    (θ, nm, κ) = par_f
    Energy_K = zeros(M_g)
    Energy_E = zeros(M_g)
    EField_T = zeros(M_g)
    p_T = zeros(M_g)
    Q_T = zeros(M_g)
    S_T = zeros(M_g)
    #E_E = 0.0
    E_mode = zeros(M_g)
    T = zeros(M_g)
    #P = zeros(J)
    ρ = zeros(J)
    S = zeros(J)
  
    for j in 1:M_g
        (Energy_K[j],Energy_E[j]) = get_energy_rel(v[:,j],(L,N,J))
        EField_T[j] = sum(v[2N+1:end,j])*dx
        p_T[j] = sum(v[N+1:2N])*dx
        get_density!(v[:,j], ρ, par_grid,0.0)
        get_current_rel!(v[:,j], S, par_grid)
        Q_T[j] = get_total_charge(ρ,(J, dx))/L # we divide by L because the density is 1, so the total charge is L, this way we compare with 1.
        S_T[j] = sum(S)/N/Q_T[j]
        T[j] = var(v[N+1:2N,j])
        E_mode = abs(rfft(v[2N+1:end,j])[nm+1])
    end
    return Energy_K, Energy_E, EField_T, p_T, Q_T, S_T, T, E_mode
  end
  
  function get_averages_threads(v,par_grid,par_evolv, par_f)
    
    (N, L, J, dx, order) = par_grid
    (t_i, t_f, M, M_g, dt) = par_evolv
    (θ, nm, κ) = par_f
  
    TS = zeros(J, nthreads())
    Tn = zeros(J, nthreads())
    
  
    par_density = (par_grid, Tn)
    par_current = (par_grid, TS)
  
    Energy_K = zeros(M_g)
    Energy_E = zeros(M_g)
    EField_T = zeros(M_g)
    p_T = zeros(M_g)
    Q_T = zeros(M_g)
    S_T = zeros(M_g)
    #E_E = 0.0
    E_mode = zeros(M_g)
    T = zeros(M_g)
    #P = zeros(J)
    ρ = zeros(J)
    S = zeros(J)
  
    for j in 1:M_g
        (Energy_K[j],Energy_E[j]) = get_energy_rel(v[:,j],(L,N,J))
        EField_T[j] = sum(v[2N+1:end,j])
        p_T[j] = sum(v[N+1:2N])*dx
        get_density_threads!(v[:,j], ρ, par_density,0.0)
        get_current_rel_threads!(v[:,j], S, par_current)
        Q_T[j] = get_total_charge(ρ,(J, dx)) / L # we divide by L because the density is 1, so the total charge is L, this way we compare with 1.
        S_T[j] = sum(S)/N/Q_T[j]
        T[j] = var(v[N+1:2N,j])
        E_mode = abs(rfft(v[2N+1:end,j])[nm+1])
    end
    return Energy_K, Energy_E, EField_T, p_T, Q_T, S_T, T, E_mode
  end
  
  
  function get_local_averages_threads(u,par_grid, par_f)
    
    (N, L, J, dx, order) = par_grid
    (θ, nm, κ) = par_f
  
    TS = zeros(J, nthreads())
    Tn = zeros(J, nthreads())
    
    par_density = (par_grid, Tn)
    par_current = (par_grid, TS)
  
    ρ = zeros(J)
    S = zeros(J)
    Efield = zeros(J)
  
    Energy_K, Energy_E = get_energy_rel(u,(L,N,J))
    Efield = u[2N+1:end]
    EField_T = sum(u[2N+1:end])*dx
    p_T = sum(u[N+1:2N])*dx
    get_density_threads!(u[:], ρ, par_density, 0.0)
    get_current_rel_threads!(u[:], S, par_current)
    Q_T = get_total_charge(ρ,(J, dx)) / L # we divide by L because the density is 1, so the total charge is L, this way we compare with 1.
    S_T = sum(S)/N/Q_T
    #T = var(u[N+1:2N])
    T = get_temperature_rel(u,N)
    E_mode = abs(rfft(u[2N+1:end])[nm+1])
  
    return ρ, S, Efield, Energy_K, Energy_E, EField_T, p_T, Q_T, S_T, T, E_mode
  end
  

  function get_local_averages_threads_D(u,par_grid, par_f)
    
    (N, J, Box, order) = par_grid
    (θ, nm, κ) = par_f
    D = length(J)
  
    TS = zeros(D,J..., nthreads())
    Tn = zeros(J..., nthreads())
    
    par_density = (par_grid, Tn)
    par_current = (par_grid, TS)
  
    ρ = zeros(J)
    S = zeros(D,J...)
    E_field = zeros(D,J...)
    B_field = zeros(J)
  
    Energy_K, Energy_E = get_energy_rel(u,(L,N,J))

    Fuj = view(u,4N+1:4N+3*prod(J))
    Fj = reshape(Fuj,3,J...)

    E_field = Fj[1:2,:,:]
    B_field = Fj[3,:,:]
    E_field_Total = sum(E_field)
    p_T = [sum(u[D+1:2D,end]), sum(u[D+2:2D,end])] /N

    

    get_density_threads_2D!(u[:], ρ, par_density, 0.0)
    get_current_threads_2D!(u[:], S, par_current)
    Q_T = get_total_charge(ρ,prod(J)) / prod(J) # we divide by L because the density is 1, so the total charge is L, this way we compare with 1.
    S_T = [sum(S[1,:,:]),sum(S[2,:,:])]/N/Q_T
    T = get_temperature_rel_D(u,N,D)
    #E_mode = abs(rfft(u[2N+1:end])[nm+1])
  
    return ρ, S, E_field, B_field, Energy_K, Energy_E, E_field_Total, p_T, Q_T, S_T, T#, E_mode
  end

  function load_averages(file_name, j, par_grid, pars_f)
      ρ, S, Efield, Energy_K, Energy_E, EField_T, p_T, Q_T, S_T, T, E_mode = get_local_averages_threads(u,par_grid, pars_f)
      tiempo = @sprintf("%05d", j)
      jldopen(file_name, "a+") do file
          file["n_$(tiempo)"] = ρ
          file["S_$(tiempo)"] = S
          file["Efield_$(tiempo)"] = Efield
          file["Energy_E_$(tiempo)"] = Energy_E
          file["Energy_K_$(tiempo)"] = Energy_K
          file["EField_T_$(tiempo)"] = EField_T
          file["p_T_$(tiempo)"] = p_T
          file["Q_T_$(tiempo)"] = Q_T
          file["S_T_$(tiempo)"] = S_T
          file["T_$(tiempo)"] = T
          file["E_mode_$(tiempo)"] = E_mode
      end
  end
  
  function load_averages_D(file_name, j, par_grid, pars_f)
    ρ, S, E_field, B_field, Energy_K, Energy_E, E_field_T, p_T, Q_T, S_T, T #=, E_mode=# = get_local_averages_threads(u,par_grid, pars_f)
    tiempo = @sprintf("%05d", j)
    jldopen(file_name, "a+") do file
        file["n_$(tiempo)"] = ρ
        file["S_$(tiempo)"] = S
        file["E_field_$(tiempo)"] = E_field
        file["B_field_$(tiempo)"] = B_field
        file["Energy_E_$(tiempo)"] = Energy_E
        file["Energy_K_$(tiempo)"] = Energy_K
        file["EField_T_$(tiempo)"] = E_field_T
        file["p_T_$(tiempo)"] = p_T
        file["Q_T_$(tiempo)"] = Q_T
        file["S_T_$(tiempo)"] = S_T
        file["T_$(tiempo)"] = T
        #file["E_mode_$(tiempo)"] = E_mode
    end
end

  function retrieve_average_data_D(data, par_grid, par_evolv; M_last=nothing)
    (N, J, Box, order) = par_grid
    (t_i, t_f, M, M_g, dt) = par_evolv
    D = length(J)
    #v = zeros(2N+J,M_g)
    n_t = zeros(J...,M_g)
    S_t = zeros(D,J...,M_g)
    E_field_t = zeros(D,J...,M_g)
    E_field_t = zeros(J...,M_g)
    E_field_T = zeros(M_g)
    Energy_K = zeros(M_g)
    Energy_E = zeros(M_g)
    #E_mode = zeros(M_g)
    p_T = zeros(D,M_g)
    Q_T = zeros(M_g)
    S_T = zeros(D,M_g)
    T = zeros(M_g)
    if M_last !== nothing # if we gave some values, then use it.
      M_g = M_last
    end
    for j in 1:M_g
        tiempo = @sprintf("%05d", j)
        n_t[:,j] = data["n_$(tiempo)"]
        S_t[:,j] = data["S_$(tiempo)"]
        E_field_t[:,j] = data["Efield_$(tiempo)"]
        B_field_t[:,j] = data["Efield_$(tiempo)"]
        Energy_K[j] = data["Energy_K_$(tiempo)"]
        Energy_E[j] = data["Energy_E_$(tiempo)"]
        E_field_T[j] = data["E_field_T_$(tiempo)"]
        p_T[j] = data["p_T_$(tiempo)"]
        Q_T[j] = data["Q_T_$(tiempo)"]
        S_T[j] = data["S_T_$(tiempo)"]
        T[j] = data["T_$(tiempo)"]
        #E_mode[j] = data["E_mode_$(tiempo)"]
    end
    return n_t, S_t, E_field_t, B_field_t, (Energy_E,  Energy_K, E_field_T, p_T, Q_T, S_T, T#=, E_mode=#)
  end
  
  function retrieve_average_data(data, par_grid, par_evolv; M_last=nothing)
    (N, L, J, dx, order) = par_grid
    (t_i, t_f, M, M_g, dt) = par_evolv
    #v = zeros(2N+J,M_g)
    n_t = zeros(J,M_g)
    S_t = zeros(J,M_g)
    Efield_t = zeros(J,M_g)
    EField_T = zeros(M_g)
    Energy_K = zeros(M_g)
    Energy_E = zeros(M_g)
    E_mode = zeros(M_g)
    p_T = zeros(M_g)
    Q_T = zeros(M_g)
    S_T = zeros(M_g)
    T = zeros(M_g)
    if M_last !== nothing # if we gave some values, then use it.
      M_g = M_last
    end
    for j in 1:M_g
        tiempo = @sprintf("%05d", j)
        n_t[:,j] = data["n_$(tiempo)"]
        S_t[:,j] = data["S_$(tiempo)"]
        Efield_t[:,j] = data["Efield_$(tiempo)"]
        Energy_K[j] = data["Energy_K_$(tiempo)"]
        Energy_E[j] = data["Energy_E_$(tiempo)"]
        EField_T[j] = data["EField_T_$(tiempo)"]
        p_T[j] = data["p_T_$(tiempo)"]
        Q_T[j] = data["Q_T_$(tiempo)"]
        S_T[j] = data["S_T_$(tiempo)"]
        T[j] = data["T_$(tiempo)"]
        E_mode[j] = data["E_mode_$(tiempo)"]
    end
    return n_t, S_t, Efield_t, (Energy_E,  Energy_K, EField_T, p_T, Q_T, S_T, T, E_mode)
  end

  function retrieve_data(data, par_grid, par_evolv; M_last=nothing)
    (N, L, J, dx, order) = par_grid
    (t_i, t_f, M, M_g, dt) = par_evolv
    v = zeros(2N+J,M_g)

    if M_last !== nothing # if we gave some values, then use it.
        M_g = M_last
    end
    
    for j in 1:M_g
        tiempo = @sprintf("%05d", j)
        v[:,j] = data["u/u_$tiempo"]
    end
    return v
  end

  function retrieve_data(data, par_grid, par_evolv; M_last=nothing, Step=1)
    (N, L, J, dx, order) = par_grid
    (t_i, t_f, M, M_g, dt) = par_evolv
    v = zeros(2N+J,M_g)
    if M_last !== nothing # if we gave some values, then use it.
        M_g = M_last
    end
    M_r = (M_last-1) ÷ Step + 1
    v = zeros(2D*N+3*prod(J),M_r)
    i = 1
    for j in 1:M_r
        tiempo = @sprintf("%05d", i)
        v[:,j] = data["u/u_$tiempo"]
        i = i+Step
    end
    return v
  end
  
  function retrieve_data_D(data, par_grid, par_evolv; M_last=nothing, Step=1)
    (N, J, Box, order) = par_grid
    (t_i, t_f, M, M_g, dt) = par_evolv
    D = length(J)

    if M_last !== nothing # if we gave some values, then use it.
        M_g = M_last
    end
    M_r = (M_last-1) ÷ Step + 1
    v = zeros(2D*N+3*prod(J),M_r)
    i = 1
    
    for j in 1:M_r
      tiempo = @sprintf("%05d", i)
      v[:,j] = data["u/u_$tiempo"]
      i = i+Step
    end
    return v
  end

  function retrieve_meta_data(file_name::String)
    data = load(file_name)
    run_name = data["run_name"]
    par_grid = data["par_grid"]
    par_evolv = data["par_evolv"]
    par_f = data["p_Ini"]
    (N, L, J, dx, order) = par_grid
    (t_i, t_f, M, M_g, dt) = par_evolv
    #@show (θ, nm, κ) = par_f
    n0 = N/L
    x = [(i-1)*dx for i in 1:J]
    dT = dt * (M-1) / (M_g-1)
    t_series = [(i-1)*dT for i in 1:M_g]
    return data, run_name, par_grid, par_evolv, par_f, n0, x, t_series
  end
  
  function retrieve_meta_data_D(file_name::String)
    data = load(file_name)
    run_name = data["run_name"]
    par_grid = data["par_grid"]
    par_evolv = data["par_evolv"]
    par_f = data["p_Ini"]
    (N, J, Box, order) = par_grid
    (t_i, t_f, M, M_g, dt) = par_evolv
    #@show (θ, nm, κ) = par_f
    n0 = N/volume(Box)
    dx = differentials(Box,J)
    
    x = [(i-1)*dx[1] for i in 1:J[1]]
    y = [(i-1)*dx[2] for i in 1:J[2]]
  
    dT = dt * (M-1) / (M_g-1)
    t_series = [t_i+(i-1)*dT for i in 1:M_g]
    return data, run_name, par_grid, par_evolv, par_f, n0, (x,y), t_series
  end
  
  