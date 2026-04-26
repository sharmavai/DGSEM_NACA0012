% =========================================================================
%  rk4_step.m
%  Classical four-stage Runge-Kutta time integrator.
%
%  FIX: This function was referenced in main_NACA0012.m but did not exist.
% =========================================================================

function Qnew = rk4_step(Q, dt, mesh, gamma, mu, k_cond)
% Inputs:
%   Q      — (n_elem, Np1, Np1, 4) conservative variables
%   dt     — time step
%   mesh   — mesh struct
%   gamma  — ratio of specific heats
%   mu     — dynamic viscosity (1/Re)
%   k_cond — thermal conductivity (mu*gamma/(Pr*(gamma-1)))
%
% Output:
%   Qnew   — updated solution

k1 = rhs_DGSEM(Q,              mesh, gamma, mu, k_cond);
k2 = rhs_DGSEM(Q + 0.5*dt*k1,  mesh, gamma, mu, k_cond);
k3 = rhs_DGSEM(Q + 0.5*dt*k2,  mesh, gamma, mu, k_cond);
k4 = rhs_DGSEM(Q + dt*k3,      mesh, gamma, mu, k_cond);

Qnew = Q + (dt/6.0) * (k1 + 2*k2 + 2*k3 + k4);
end
