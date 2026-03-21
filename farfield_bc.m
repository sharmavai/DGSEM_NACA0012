% =========================================================================
%  farfield_bc.m
%  Far-field boundary condition via Riemann invariants (characteristic BC)
%
%  FIX #5a Applied:
%    Previous code used simple Dirichlet (Q_ghost = Q_freestream), which
%    reflects spurious acoustic waves off the far-field boundary and
%    pollutes the entire domain — especially bad for unsteady Re=1e4 flow.
%
%    CORRECT: Split characteristics at boundary.
%      Subsonic inflow  (Vn < 0): R+ from interior, R- from freestream
%      Subsonic outflow (Vn > 0): R+ from interior, R- from interior (?)
%      => entropy from freestream on inflow, from interior on outflow
%    This is Thompson (1987) / Poinsot & Lele (1992) for Euler;
%    extended straightforwardly to NS ghost-cell approach.
%
%  References:
%    [1] Thompson (1987) J. Comput. Phys. 68 — characteristic BCs
%    [2] Kopriva (2009) Sec 9.2 — far-field BC for compressible flow
%    [3] Ferrer et al. (2021) — same approach for DGSEM NACA0012
% =========================================================================

function Q_ghost = farfield_bc(Q_int, nx, ny, Mach_inf, alpha_deg, gamma)
% Inputs:
%   Q_int     — (4x1) interior conservative state at face node
%   nx, ny    — outward unit normal (pointing OUT of domain, into ghost)
%   Mach_inf  — freestream Mach number
%   alpha_deg — angle of attack [degrees]
%   gamma     — ratio of specific heats
%
% Output:
%   Q_ghost   — (4x1) ghost cell state for Roe flux computation

    alpha   = alpha_deg * pi / 180.0;

    %% ── Freestream state (non-dimensional) ───────────────────────────────
    rho_inf = 1.0;
    p_inf   = 1.0 / gamma;
    u_inf   = Mach_inf * cos(alpha);
    v_inf   = Mach_inf * sin(alpha);
    c_inf   = sqrt(gamma * p_inf / rho_inf);   % = 1.0 by construction
    Vn_inf  = u_inf*nx + v_inf*ny;

    %% ── Interior state ───────────────────────────────────────────────────
    rho_i = Q_int(1);
    u_i   = Q_int(2) / rho_i;
    v_i   = Q_int(3) / rho_i;
    E_i   = Q_int(4);
    p_i   = (gamma-1) * (E_i - 0.5*rho_i*(u_i^2 + v_i^2));
    p_i   = max(p_i, 1e-12);   % pressure floor
    c_i   = sqrt(gamma * p_i / rho_i);
    Vn_i  = u_i*nx + v_i*ny;

    %% ── Riemann invariants ───────────────────────────────────────────────
    % R+ = Vn + 2c/(gamma-1)  travels at lambda = Vn + c  (from interior)
    % R- = Vn - 2c/(gamma-1)  travels at lambda = Vn - c
    Rplus  =  Vn_i   + 2.0*c_i   / (gamma-1);   % from interior
    Rminus =  Vn_inf - 2.0*c_inf / (gamma-1);   % from freestream

    Vn_b = 0.5*(Rplus + Rminus);
    c_b  = 0.25*(gamma-1)*(Rplus - Rminus);

    if c_b < 0
        % Unphysical — fall back to freestream (very rare, only if BL too small)
        Q_ghost = [rho_inf; rho_inf*u_inf; rho_inf*v_inf; ...
                   p_inf/(gamma-1)+0.5*rho_inf*(u_inf^2+v_inf^2)];
        return;
    end

    %% ── Entropy: from freestream (inflow) or interior (outflow) ─────────
    if Vn_i <= 0.0   % inflow: characteristics enter from freestream
        s_b = p_inf / rho_inf^gamma;
    else              % outflow: entropy advected from interior
        s_b = p_i   / rho_i^gamma;
    end

    %% ── Recover primitive variables ──────────────────────────────────────
    rho_b = (c_b^2 / (gamma * s_b))^(1.0/(gamma-1));
    p_b   = s_b * rho_b^gamma;

    % Tangential velocity: from freestream (maintains correct flow angle)
    % Decompose freestream into normal + tangential
    % t = (-ny, nx)  (90-deg rotation of n)
    tx = -ny;  ty = nx;
    Vt_inf = u_inf*tx + v_inf*ty;

    u_b = Vn_b*nx + Vt_inf*tx;
    v_b = Vn_b*ny + Vt_inf*ty;
    E_b = p_b/(gamma-1) + 0.5*rho_b*(u_b^2 + v_b^2);

    Q_ghost = [rho_b; rho_b*u_b; rho_b*v_b; E_b];
end
