function [Cp, x_surf, y_surf] = compute_pressure_coefficient(U, mesh, N, gamma, rho_inf, U_inf)
% COMPUTE_PRESSURE_COEFFICIENT  Surface Cp distribution.
%   Cp = (p - p_inf) / (0.5 * rho_inf * U_inf^2)

q_inf = 0.5 * rho_inf * U_inf^2;
p_inf = rho_inf * U_inf^2 / (gamma * 1.4^2);

n_surf = length(mesh.airfoil_ids) * (N+1);
Cp     = zeros(n_surf, 1);
x_surf = zeros(n_surf, 1);
y_surf = zeros(n_surf, 1);

idx = 1;
for e = mesh.airfoil_ids
    for j = 1:N+1
        rho  = U(e,1,j,1);
        rhou = U(e,1,j,2);
        rhov = U(e,1,j,3);
        rhoE = U(e,1,j,4);
        u = rhou/rho; v = rhov/rho; E = rhoE/rho;
        p = (gamma-1)*rho*(E - 0.5*(u^2+v^2));

        Cp(idx)     = (p - p_inf) / q_inf;
        x_surf(idx) = 0.5;
        y_surf(idx) = 0.0;
        idx = idx + 1;
    end
end
end
