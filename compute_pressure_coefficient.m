% =========================================================================
%  compute_pressure_coefficient.m
%  Surface Cp distribution on the NACA 0012 airfoil.
%
%  FIX Applied:
%    - p_inf formula was wrong: used rho_inf*U_inf^2/(gamma*1.4^2)
%      which has a hardcoded 1.4^2.  Corrected to 1/gamma (non-dim).
%    - Surface coordinates were hardcoded to (0.5, 0.0) for ALL nodes,
%      making any Cp-vs-x plot useless.  Corrected to use actual mesh
%      wall coordinates from mesh.wall_x, mesh.wall_y.
%    - Updated to use (n_elem, Np1, Np1, 4) Q-layout and correct
%      mesh field names (airfoil_ei, wall_x, wall_y).
% =========================================================================

function [Cp, x_surf, y_surf] = compute_pressure_coefficient(Q, mesh, gamma, Mach, N)
% Inputs:
%   Q      — (n_elem, Np1, Np1, 4) conservative variables
%   mesh   — mesh struct from generate_mesh_NACA0012
%   gamma  — ratio of specific heats
%   Mach   — freestream Mach number
%   N      — polynomial degree
%
% Outputs:
%   Cp     — pressure coefficient at each surface node
%   x_surf — x-coordinates of surface nodes
%   y_surf — y-coordinates of surface nodes

Np1 = N + 1;

% Non-dimensional reference state
q_inf = 0.5 * 1.0 * Mach^2;      % 0.5 * rho_inf * U_inf^2
p_inf = 1.0 / gamma;              % FIX: was rho_inf*U_inf^2/(gamma*1.4^2)

airfoil_ei = mesh.airfoil_ei;
n_airfoil  = length(airfoil_ei);
n_surf     = n_airfoil * Np1;

Cp     = zeros(n_surf, 1);
x_surf = zeros(n_surf, 1);
y_surf = zeros(n_surf, 1);

idx = 0;
for ii = 1:n_airfoil
    ei = airfoil_ei(ii);
    e  = ei;                       % first radial layer (ej=1)

    for jj = 1:Np1
        idx = idx + 1;

        rho  = max(Q(e, jj, 1, 1), 1e-12);
        u_v  = Q(e, jj, 1, 2) / rho;
        v_v  = Q(e, jj, 1, 3) / rho;
        E_v  = Q(e, jj, 1, 4) / rho;
        p    = (gamma-1) * rho * (E_v - 0.5*(u_v^2 + v_v^2));

        Cp(idx) = (p - p_inf) / q_inf;

        % FIX: use actual mesh surface coordinates (was 0.5, 0.0 for ALL)
        x_surf(idx) = mesh.wall_x(jj, ei);
        y_surf(idx) = mesh.wall_y(jj, ei);
    end
end
end
