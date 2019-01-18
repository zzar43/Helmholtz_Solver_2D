
using SparseArrays, LinearAlgebra;

# Make the differential operator A
function make_diff_op(vel, f, fre, h, Nx_pml, Ny_pml, pml_len, pml_coef)

    pml_value = range(0, stop=pml_coef, length=pml_len);
    
    beta = zeros(Nx_pml, Ny_pml);
    for i = 1:pml_len
        beta[pml_len+1-i,:] .= pml_value[i];
        beta[end-pml_len+i,:] .= pml_value[i];
        beta[:,pml_len+1-i] .= pml_value[i];
        beta[:,end-pml_len+i] .= pml_value[i];
    end

    vel_ex = zeros(Nx_pml, Ny_pml);
    vel_ex[pml_len+1:end-pml_len, pml_len+1:end-pml_len] .= vel;
    for i = 1:pml_len
        vel_ex[i,:] = vel_ex[pml_len+1,:];
        vel_ex[end-i+1,:] = vel_ex[end-pml_len,:];
        vel_ex[:,i] = vel_ex[:,pml_len+1];
        vel_ex[:,end-i+1] = vel_ex[:,end-pml_len];
    end

    coef = (1 .+ im*beta[:]) .* (h^2*(2*pi*fre).^2) ./ (vel_ex[:].^2) .- 4;
    vec1 = ones(Nx_pml*Ny_pml-Nx_pml);
    vec2 = ones(Nx_pml*Ny_pml-1);
    A = spdiagm(-Nx_pml => vec1, -1 => vec2, 0 => coef, 1 => vec2, Nx_pml => vec1);
    for i = 1:(Ny_pml-1)
            ind_x = Nx_pml*i+1;
            ind_y = Nx_pml*i;
            A[ind_x,ind_y] = 0;
            A[ind_y,ind_x] = 0;
    end
    A = A ./ h^2;
    
    return A
end

# Make the source vector
function make_source_vec(f, Nx_pml, Ny_pml, pml_len)
    
    f_ex = zeros(Nx_pml, Ny_pml);
    f_ex[pml_len+1:end-pml_len, pml_len+1:end-pml_len] .= f;
    f_ex = f_ex[:];
    
    return f_ex
end

# Compute the wavefield of Au = f.
function compute_wavefield(A, f_ex, Nx, Ny, Nx_pml, Ny_pml, pml_len)
    
    u = A \ f_ex;
    u = reshape(u, Nx_pml, Ny_pml);
    u = u[pml_len+1:end-pml_len, pml_len+1:end-pml_len];
    
    return u
end

# Main solver
function scalar_helmholtz_solver_2d(vel, f, fre, h, Nx, Ny)
    
    # PML
    pml_len = 20;
    pml_coef = maximum(abs.(f));
    Nx_pml = Nx + 2*pml_len;
    Ny_pml = Ny + 2*pml_len;
    
    A = make_diff_op(vel, f, fre, h, Nx_pml, Ny_pml, pml_len, pml_coef)
    
    f_ex = make_source_vec(f, Nx_pml, Ny_pml, pml_len);
    
    u = compute_wavefield(A, f_ex, Nx, Ny, Nx_pml, Ny_pml, pml_len);
    
    return u
end