% =========================================================================
%  compute_timestep.m
%  CFL-based time step for DGSEM (inviscid + viscous limits)
%
%  FIX Applied:
%    - Function signature was (U, mesh, N, CFL, gamma) — 5 args
%      but all callers passed 7 args (U, mesh, CFL, CFL_v, gamma, mu, N).
%      Corrected to match caller convention.
%    - Mesh indexing was mesh.x(i,j,corner) assuming a corner-based layout.
%      mesh.x is actually (Np1, Np1, n_elem).  Corrected to use
%      mesh.h_min(e) which is pre-computed per element.
%    - Added viscous CFL (Fourier) limit:  dt_vis = CFL_v * h^2 / (nu*(2N+1)^2)
%    - Loops over all elements correctly using mesh.n_elem.
% =========================================================================

function dt = compute_timestep(Q, mesh, CFL, CFL_v, gamma, mu, N)
% Inputs:
%   Q     — (n_elem, Np1, Np1, 4) conservative variables
%   mesh  — mesh struct
%   CFL   — inviscid CFL number
%   CFL_v — viscous CFL number (Fourier number)
%   gamma — ratio of specific heats
%   mu    — dynamic viscosity (non-dim: 1/Re)
%   N     — polynomial degree

Np1    = N + 1;
n_elem = mesh.n_elem;
dt     = inf;

for e = 1:n_elem
    % Max wave speed in element
    lam_max = 0;
    for jj = 1:Np1
        for ii = 1:Np1
            rho = max(Q(e, ii, jj, 1), 1e-12);
            u   = Q(e, ii, jj, 2) / rho;
            v   = Q(e, ii, jj, 3) / rho;
            E   = Q(e, ii, jj, 4) / rho;
            p   = max((gamma-1)*rho*(E - 0.5*(u^2+v^2)), 1e-12);
            c   = sqrt(gamma * p / rho);
            lam_max = max(lam_max, sqrt(u^2 + v^2) + c);
        end
    end
    lam_max = max(lam_max, 1e-12);

    h = mesh.h_min(e);

    % Inviscid CFL limit
    dt_inv = CFL * h / (lam_max * (2*N+1));

    % Viscous (Fourier) limit: dt <= CFL_v * h^2 / (nu * (2N+1)^2)
    if mu > 0
        dt_vis = CFL_v * h^2 / (mu * (2*N+1)^2);
    else
        dt_vis = inf;
    end

    dt = min(dt, min(dt_inv, dt_vis));
end

if ~isfinite(dt) || dt <= 0
    error('Invalid time step: dt = %e', dt);
end
end
