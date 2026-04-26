% =========================================================================
%  initialize_solution.m
%  Free-stream initial condition for 2D compressible NS (non-dimensional)
%
%  FIX Applied:
%    - v_inf was 0 regardless of alpha — every AoA run silently solved
%      the alpha=0 case.  CORRECTED to u_inf=M*cos(alpha), v_inf=M*sin(alpha).
%    - Q layout changed from (4, Np1^2, n_elem) to (n_elem, Np1, Np1, 4)
%      to match rhs_DGSEM, compute_timestep, and all other functions.
%    - Previous layout caused a dimension mismatch crash at first RK step.
%
%  Non-dimensionalisation:
%    rho_inf=1, c_inf=1 => p_inf=1/gamma, V_inf=Mach
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
%   Q — (n_elem, Np1, Np1, 4) conservative variables
%       Q(:,:,:,1) = rho,  Q(:,:,:,2) = rho*u,
%       Q(:,:,:,3) = rho*v,  Q(:,:,:,4) = E

    Np1    = N + 1;
    n_elem = mesh.n_elem;
    alpha  = alpha_deg * pi / 180.0;

    %% ── Non-dimensional free-stream state ────────────────────────────────
    rho_inf = 1.0;
    p_inf   = 1.0 / gamma;

    u_inf   = Mach * cos(alpha);
    v_inf   = Mach * sin(alpha);

    E_inf   = p_inf/(gamma - 1.0) + 0.5*rho_inf*(u_inf^2 + v_inf^2);

    %% ── Fill conservative variable array ─────────────────────────────────
    Q = zeros(n_elem, Np1, Np1, 4);

    Q(:,:,:,1) = rho_inf;
    Q(:,:,:,2) = rho_inf * u_inf;
    Q(:,:,:,3) = rho_inf * v_inf;
    Q(:,:,:,4) = E_inf;

    %% ── Diagnostic output ────────────────────────────────────────────────
    fprintf('[init] M=%.3f  alpha=%.2f deg  u_inf=%.6f  v_inf=%.6f\n', ...
            Mach, alpha_deg, u_inf, v_inf);
    fprintf('[init] rho=%.4f  p=%.6f  E=%.6f\n', rho_inf, p_inf, E_inf);

    V_check = sqrt(u_inf^2 + v_inf^2);
    c_check = sqrt(gamma * p_inf / rho_inf);
    M_check = V_check / c_check;
    assert(abs(M_check - Mach) < 1e-12, ...
           'Mach number mismatch: set=%.6f  recovered=%.6f', Mach, M_check);
end
