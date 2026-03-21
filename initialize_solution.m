% =========================================================================
%  initialize_solution.m
%  Free-stream initial condition for 2D compressible NS (non-dimensional)
%
%  FIX #4 Applied:
%    - v_inf was 0 regardless of alpha — every AoA run silently solved
%      the alpha=0 case; solver would slowly drift but never reach correct AoA.
%    - CORRECTED: u_inf = M*cos(alpha),  v_inf = M*sin(alpha)
%    - Non-dimensionalisation made explicit:
%        rho_inf = 1,  c_inf = 1  =>  p_inf = 1/gamma,  V_inf = Mach
%    - Q laid out as (4 x Np1*Np1 x n_elem) matching rhs_DGSEM convention
%
%  References:
%    [1] Kopriva (2009) Ch. 9 — non-dimensional NS free-stream conditions
%    [2] Ferrer et al. (2021) — M=0.3, Re=1e4, alpha sweep, same convention
% =========================================================================

function Q = initialize_solution(mesh, Mach, alpha_deg, gamma, N)
% Inputs:
%   mesh      — mesh struct from generate_mesh_NACA0012
%   Mach      — freestream Mach number
%   alpha_deg — angle of attack [degrees]
%   gamma     — ratio of specific heats (1.4)
%   N         — polynomial degree
%
% Output:
%   Q — (4 x Np1^2 x n_elem) conservative variable array
%       Q(1,:,:) = rho
%       Q(2,:,:) = rho*u
%       Q(3,:,:) = rho*v
%       Q(4,:,:) = E  (total energy per unit volume)

    Np1    = N + 1;
    n_elem = mesh.n_elem;
    alpha  = alpha_deg * pi / 180.0;

    %% ── Non-dimensional free-stream state ────────────────────────────────
    % Reference: rho_inf=1, c_inf=1, chord=1
    % => p_inf = rho*c^2/gamma = 1/gamma
    % => V_inf = Mach * c_inf  = Mach
    rho_inf =  1.0;
    p_inf   =  1.0 / gamma;

    % FIX #4: Apply angle of attack correctly to BOTH velocity components
    % WRONG (previous):  u_inf = Mach;  v_inf = 0;
    % CORRECT:
    u_inf   =  Mach * cos(alpha);
    v_inf   =  Mach * sin(alpha);

    % Total energy per unit volume
    E_inf   =  p_inf/(gamma - 1.0) + 0.5*rho_inf*(u_inf^2 + v_inf^2);

    %% ── Fill conservative variable array ─────────────────────────────────
    n_nodes = Np1^2;
    Q = zeros(4, n_nodes, n_elem);

    Q(1, :, :) = rho_inf;
    Q(2, :, :) = rho_inf * u_inf;
    Q(3, :, :) = rho_inf * v_inf;
    Q(4, :, :) = E_inf;

    %% ── Diagnostic output ────────────────────────────────────────────────
    fprintf('[init] M=%.3f  alpha=%.2f deg  u_inf=%.6f  v_inf=%.6f\n', ...
            Mach, alpha_deg, u_inf, v_inf);
    fprintf('[init] rho=%.4f  p=%.6f  E=%.6f\n', rho_inf, p_inf, E_inf);

    % Sanity: recover Mach from velocity
    V_check = sqrt(u_inf^2 + v_inf^2);
    c_check = sqrt(gamma * p_inf / rho_inf);
    M_check = V_check / c_check;
    assert(abs(M_check - Mach) < 1e-12, ...
           'Mach number mismatch: set=%.6f  recovered=%.6f', Mach, M_check);
end