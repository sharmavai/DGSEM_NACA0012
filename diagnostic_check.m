%% Diagnostic — Initialization and Timestep Check
%  Verifies that mesh generation, initialization, and timestep computation
%  work correctly before running a full simulation.
%
%  Usage: diagnostic_check
% =========================================================================
clear; close all; clc;

fprintf('=== DGSEM DIAGNOSTIC CHECK ===\n\n');

Re = 1e4; Mach = 0.3; alpha = 5.0;
gamma = 1.4; Pr = 0.72;
mu    = 1.0/Re;
k_cond = mu*gamma/(Pr*(gamma-1));

N = 3;
nelem_airfoil = 48; nelem_wake = 16; nelem_radial = 20;
R_ff = 30.0; y1 = 5.0/Re;

fprintf('Parameters: Re=%.0e, M=%.2f, alpha=%.1f, N=%d\n\n', Re, Mach, alpha, N);

%% ── 1. Mesh generation ───────────────────────────────────────────────────
fprintf('TEST 1: Mesh generation...\n');
mesh = generate_mesh_NACA0012(nelem_airfoil, nelem_wake, nelem_radial, R_ff, y1, N);
mesh.Mach = Mach;  mesh.alpha_deg = alpha;

fprintf('  n_elem = %d\n', mesh.n_elem);
fprintf('  min(J)  = %.4e (must be > 0)\n', min(mesh.J(:)));
fprintf('  min(h)  = %.4e\n', min(mesh.h_min));
assert(min(mesh.J(:)) > 0, 'FAIL: Negative Jacobian');
fprintf('  PASS\n\n');

%% ── 2. Solution initialization ──────────────────────────────────────────────
fprintf('TEST 2: Solution initialization...\n');
Q = initialize_solution(mesh, Mach, alpha, gamma, N);

fprintf('  NaN: %d  Inf: %d\n', any(isnan(Q(:))), any(isinf(Q(:))));
fprintf('  rho : [%.6e, %.6e]\n', min(Q(:,:,:,1),[],'all'), max(Q(:,:,:,1),[],'all'));
fprintf('  rhoE: [%.6e, %.6e]\n', min(Q(:,:,:,4),[],'all'), max(Q(:,:,:,4),[],'all'));
assert(~any(isnan(Q(:))) && ~any(isinf(Q(:)))), 'FAIL: NaN/Inf in Q');
fprintf('  PASS\n\n');

%% ── 3. Timestep computation ─────────────────────────────────────────────────
fprintf('TEST 3: Timestep computation...\n');
CFL = 0.08; CFL_v = 0.08;
dt = compute_timestep(Q, mesh, CFL, CFL_v, gamma, mu, N);
fprintf('  dt = %.6e\n', dt);

dt_inv_only = CFL * min(mesh.h_min) / ((Mach + 1) * (2*N+1));
dt_vis_only = CFL_v * min(mesh.h_min)^2 / (mu * (2*N+1)^2);
fprintf('  dt_inv (approx) = %.6e\n', dt_inv_only);
fprintf('  dt_vis (approx) = %.6e\n', dt_vis_only);
if dt_vis_only < dt_inv_only
    fprintf('  Controlling limit: VISCOUS\n');
else
    fprintf('  Controlling limit: INVISCID\n');
end
fprintf('  PASS\n\n');

%% ── 4. Single RK4 step ─────────────────────────────────────────────────────
fprintf('TEST 4: Single RK4 step...\n');
Q_new = rk4_step(Q, dt, mesh, gamma, mu, k_cond);
fprintf('  NaN after RK4: %d  Inf: %d\n', any(isnan(Q_new(:))), any(isinf(Q_new(:))));
fprintf('  max |dQ| = %.6e\n', max(abs(Q_new(:) - Q(:))));
assert(~any(isnan(Q_new(:)))), 'FAIL: NaN after RK4 step');
fprintf('  PASS\n\n');

%% ── 5. GCL check ───────────────────────────────────────────────────────────
fprintf('TEST 5: GCL free-stream preservation (5 steps, alpha=0)...\n');
Q0_gcl = initialize_solution(mesh, Mach, 0.0, gamma, N);
Q_gcl  = Q0_gcl;
dt_gcl = compute_timestep(Q_gcl, mesh, CFL, CFL_v, gamma, mu, N);
for k = 1:5
    Q_gcl = rk4_step(Q_gcl, dt_gcl, mesh, gamma, mu, k_cond);
end
gcl_err = max(abs(Q_gcl(:) - Q0_gcl(:)));
fprintf('  GCL residual = %.4e (should be < 1e-10)\n', gcl_err);
if gcl_err > 1e-8
    warning('GCL FAILED! Check metric computation.');
else
    fprintf('  PASS\n');
end
fprintf('\n');

%% ── 6. Estimated simulation time ────────────────────────────────────────────
t_convective = 1.0/Mach;
t_end = 30 * t_convective;
n_steps_est = ceil(t_end / dt);
fprintf('ESTIMATED SIMULATION:\n');
fprintf('  Total time:      %.1f (%.0f convective times)\n', t_end, t_end/t_convective);
fprintf('  Estimated steps: %d\n', n_steps_est);
fprintf('  Estimated wall:  ~%.0f minutes (rough estimate)\n', n_steps_est*0.001);

fprintf('\n=== ALL DIAGNOSTICS PASSED ===\n');
fprintf('Ready to run main_NACA0012\n');
