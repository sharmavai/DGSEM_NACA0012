%% Diagnostic — Initialization and Timestep Check
%  FIX: Updated calling conventions to match corrected function signatures.
clear; close all; clc;

Re = 1e4; Mach = 0.3; alpha = 0.0;
gamma = 1.4; Pr = 0.72;
mu    = 1.0/Re;
k_cond = mu*gamma/(Pr*(gamma-1));

N = 3;
nelem_airfoil = 40; nelem_wake = 20; nelem_radial = 15;
R_ff = 30.0; y1 = 5.0/Re;

mesh = generate_mesh_NACA0012(nelem_airfoil, nelem_wake, nelem_radial, R_ff, y1, N);
mesh.Mach = Mach;  mesh.alpha_deg = alpha;

Q = initialize_solution(mesh, Mach, alpha, gamma, N);

fprintf('NaN: %d  Inf: %d\n', any(isnan(Q(:))), any(isinf(Q(:))));
fprintf('rho : [%.3e, %.3e]\n', min(Q(:,:,:,1),[],'all'), max(Q(:,:,:,1),[],'all'));
fprintf('rhoE: [%.3e, %.3e]\n', min(Q(:,:,:,4),[],'all'), max(Q(:,:,:,4),[],'all'));

CFL = 0.1; CFL_v = 0.1;
dt = compute_timestep(Q, mesh, CFL, CFL_v, gamma, mu, N);
fprintf('dt = %.6e\n', dt);
fprintf('Estimated iters (~10 convective times): %.0f\n', 10/Mach/dt);
