% =========================================================================
%  compute_pressure_coefficient.m
%  Surface Cp distribution on the NACA 0012 airfoil.
%
%  Uses mesh wall coordinates and LGL quadrature for accurate
%  surface pressure distribution.
%
%  References:
%    [1] Gregory & O'Reilly (1970) R&M 3726 — NACA 0012 Cp data
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
p_inf = 1.0 / gamma;              % Non-dimensional freestream pressure

airfoil_ei = mesh.airfoil_ei;
wall_elems = mesh.wall_elems;
n_airfoil  = length(airfoil_ei);
n_surf     = n_airfoil * Np1;

Cp     = zeros(n_surf, 1);
x_surf = zeros(n_surf, 1);
y_surf = zeros(n_surf, 1);

idx = 0;
for ii_idx = 1:n_airfoil
    ei = airfoil_ei(ii_idx);
    e  = wall_elems(ei);     % Correct global element ID

    for jj = 1:Np1
        idx = idx + 1;

        rho  = max(Q(e, jj, 1, 1), 1e-12);
        u_v  = Q(e, jj, 1, 2) / rho;
        v_v  = Q(e, jj, 1, 3) / rho;
        E_v  = Q(e, jj, 1, 4) / rho;
        p    = (gamma-1) * rho * (E_v - 0.5*(u_v^2 + v_v^2));

        Cp(idx) = (p - p_inf) / q_inf;

        x_surf(idx) = mesh.wall_x(jj, ei);
        y_surf(idx) = mesh.wall_y(jj, ei);
    end
end
end
