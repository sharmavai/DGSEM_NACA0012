% =========================================================================
%  wall_bc.m
%  No-slip adiabatic wall boundary condition (ghost-cell approach)
%
%  FIX #5b Applied:
%    Previous wall BC was approximate — velocity not properly mirrored,
%    pressure gradient not enforced, leading to non-zero wall velocity
%    and incorrect boundary layer behaviour.
%
%    CORRECT ghost-cell construction:
%      - Mirror tangential AND normal velocity components to enforce u_wall=0
%        (average of interior and ghost = 0 at face => ghost = -interior)
%      - Pressure: zero normal gradient => p_ghost = p_interior  (extrapolation)
%      - Density: from pressure and adiabatic condition
%      - Adiabatic: dT/dn = 0  =>  T_ghost = T_interior
%
%  References:
%    [1] Blazek (2015) Sec 8.4 — no-slip wall ghost cell
%    [2] Kopriva (2009) Sec 9.2 — wall BC for compressible NS
% =========================================================================

function Q_ghost = wall_bc(Q_int, gamma)
% Inputs:
%   Q_int   — (4x1) interior conservative state at wall face node
%   gamma   — ratio of specific heats
%
% Output:
%   Q_ghost — (4x1) ghost cell state
%             When averaged with Q_int at the face: u_face=0, v_face=0

    rho_i = Q_int(1);
    u_i   = Q_int(2) / rho_i;
    v_i   = Q_int(3) / rho_i;
    E_i   = Q_int(4);
    p_i   = (gamma-1) * (E_i - 0.5*rho_i*(u_i^2 + v_i^2));
    p_i   = max(p_i, 1e-12);

    % Adiabatic: T_ghost = T_int  =>  rho_ghost = rho_int (same p, same T)
    rho_g = rho_i;
    p_g   = p_i;             % zero normal pressure gradient

    % No-slip: mirror velocity  =>  u_ghost = -u_int, v_ghost = -v_int
    % Then at face: u_face = 0.5*(u_int + u_ghost) = 0  ✓
    u_g   = -u_i;
    v_g   = -v_i;

    E_g   = p_g/(gamma-1) + 0.5*rho_g*(u_g^2 + v_g^2);

    Q_ghost = [rho_g; rho_g*u_g; rho_g*v_g; E_g];
end
