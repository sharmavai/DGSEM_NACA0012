% =========================================================================
%  compute_timestep.m
%  CFL-based time step for DGSEM (inviscid + viscous limits)
%
%  Takes the minimum of:
%    dt_inv = CFL * h / (lambda_max * (2N+1))         [inviscid]
%    dt_vis = CFL_v * h^2 / (nu * (2N+1)^2)           [viscous Fourier]
%
%  References:
%    [1] Kopriva (2009) Ch. 8 — CFL limits for DGSEM
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
%
% Output:
%   dt    — stable time step

Np1    = N + 1;
n_elem = mesh.n_elem;
dt     = inf;

for e = 1:n_elem
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

    % Viscous (Fourier) limit
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
