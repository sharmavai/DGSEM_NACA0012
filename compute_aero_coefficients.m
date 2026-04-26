% =========================================================================
%  compute_aero_coefficients.m
%  Surface pressure integration on airfoil wall elements.
%
%  FIX Applied:
%    - Used mesh.airfoil_ids which doesn't exist; corrected to
%      mesh.airfoil_ei (circumferential indices of airfoil elements).
%    - Used mesh.x(1,j,1) treating mesh.x as corner-based; mesh.x is
%      (Np1, Np1, n_elem).  Corrected to use mesh.wall_x/wall_y/wall_nx
%      /wall_ny/wall_ds which are pre-computed for the wall face.
%    - Surface integration used ds = seg/(N+1) (trapezoidal); corrected
%      to use mesh.wall_ds which includes LGL quadrature weights.
%    - Lift and drag now correctly rotated to freestream coordinates
%      for non-zero angle of attack.
%    - Q-layout updated to (n_elem, Np1, Np1, 4).
% =========================================================================

function [CL, CD, CM] = compute_aero_coefficients(Q, mesh, Mach, alpha_deg, gamma, Re, N)
% Inputs:
%   Q         — (n_elem, Np1, Np1, 4) conservative variables
%   mesh      — mesh struct from generate_mesh_NACA0012
%   Mach      — freestream Mach number
%   alpha_deg — angle of attack [degrees]
%   gamma     — ratio of specific heats
%   Re        — Reynolds number (unused here, kept for interface)
%   N         — polynomial degree

Np1   = N + 1;
alpha = alpha_deg * pi / 180.0;

% Non-dimensional reference: rho_inf=1, U_inf=Mach, chord=1
q_inf = 0.5 * 1.0 * Mach^2;
chord = 1.0;

Fx = 0;  Fy = 0;  Mt = 0;

airfoil_ei = mesh.airfoil_ei;

for idx = 1:length(airfoil_ei)
    ei = airfoil_ei(idx);
    e  = ei;   % element ID (first radial layer, ej=1)

    for ii = 1:Np1
        % Wall face at eta=-1 (j=1)
        rho = max(Q(e, ii, 1, 1), 1e-12);
        u_v = Q(e, ii, 1, 2) / rho;
        v_v = Q(e, ii, 1, 3) / rho;
        E_v = Q(e, ii, 1, 4) / rho;
        p   = (gamma-1) * rho * (E_v - 0.5*(u_v^2 + v_v^2));
        if ~isfinite(p) || p < 0, continue; end

        % wall_nx, wall_ny point INTO FLUID (outward from airfoil)
        nx_w = mesh.wall_nx(ii, ei);
        ny_w = mesh.wall_ny(ii, ei);
        ds_w = mesh.wall_ds(ii, ei);

        % Force on airfoil = -p * n_outward * ds
        Fx = Fx - p * nx_w * ds_w;
        Fy = Fy - p * ny_w * ds_w;

        % Moment about quarter chord (0.25c, 0)
        xs = mesh.wall_x(ii, ei);
        ys = mesh.wall_y(ii, ei);
        Mt = Mt - p * ((xs - 0.25*chord)*ny_w - ys*nx_w) * ds_w;
    end
end

% Rotate to freestream coordinates
q_chord = q_inf * chord;
CL = (-Fx*sin(alpha) + Fy*cos(alpha)) / q_chord;
CD = ( Fx*cos(alpha) + Fy*sin(alpha)) / q_chord;
CM = Mt / (q_chord * chord);

if ~isfinite(CL), CL = 0; end
if ~isfinite(CD), CD = 0; end
if ~isfinite(CM), CM = 0; end
end
