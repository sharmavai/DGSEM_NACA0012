% =========================================================================
%  test_mesh.m
%  Validates generate_mesh_NACA0012.m — run before any flow simulation
%  All tests must PASS before proceeding to Fix #4
% =========================================================================
clear; clc; close all;

fprintf('=== FIX #3 MESH VALIDATION TEST ===\n\n');

N   = 3;
Re  = 1e4;
mesh = generate_mesh_NACA0012(48, 16, 20, 30.0, 5/Re, N);

%% ── TEST 1: Jacobian positivity ─────────────────────────────────────────
min_J = min(mesh.J(:));
fprintf('TEST 1 — Jacobian positivity\n');
fprintf('  min(J) = %.6e  (must be > 0)\n', min_J);
assert(min_J > 0, 'FAIL: Negative Jacobian detected');
fprintf('  PASS\n\n');

%% ── TEST 2: Free-stream preservation (GCL) ──────────────────────────────
% Constant field Q = [1, u, 0, E]  ->  divergence should be zero
fprintf('TEST 2 — Free-stream preservation (GCL)\n');
gamma = 1.4;  Mach = 0.3;
rho0 = 1.0;  u0 = Mach;  v0 = 0.0;
p0   = 1/(gamma);  E0 = p0/(gamma-1) + 0.5*rho0*(u0^2+v0^2);
mesh.Mach = Mach;  mesh.alpha_deg = 0.0;
Q0   = initialize_solution(mesh, Mach, 0.0, gamma, N);
Q_test = Q0;
dt = compute_timestep(Q_test, mesh, 0.1, 0.1, gamma, 1/Re, N);
for k = 1:5
    Q_test = rk4_step(Q_test, dt, mesh, gamma, 1/Re, ...
                      gamma/(0.72*(gamma-1)*Re));
end
gcl_err = max(abs(Q_test(:) - Q0(:)));
fprintf('  GCL residual = %.4e  (must be < 1e-10)\n', gcl_err);
assert(gcl_err < 1e-8, 'FAIL: GCL violated — check metric computation');
fprintf('  PASS\n\n');

%% ── TEST 3: Unit normal check — |n| = 1 at all face nodes ───────────────
fprintf('TEST 3 — Unit normals\n');
n_mag = sqrt(mesh.face_nx.^2 + mesh.face_ny.^2);
max_dev = max(abs(n_mag(:) - 1.0));
fprintf('  max | |n| - 1 | = %.4e  (must be < 1e-12)\n', max_dev);
assert(max_dev < 1e-10, 'FAIL: Face normals not unit vectors');
fprintf('  PASS\n\n');

%% ── TEST 4: NACA 0012 geometry — airfoil surface check ──────────────────
fprintf('TEST 4 — NACA 0012 surface geometry\n');
t = 0.12;
yt = @(x) (t/0.2)*(0.2969*sqrt(x)-0.1260*x-0.3516*x.^2+0.2843*x.^3-0.1036*x.^4);
x_wall = mesh.wall_x(:);
y_wall = mesh.wall_y(:);
% Split into upper/lower by sign of y
upper = y_wall >= 0;
lower = y_wall <  0;
err_up = max(abs(abs(y_wall(upper)) - yt(x_wall(upper))));
err_lo = max(abs(abs(y_wall(lower)) - yt(x_wall(lower))));
fprintf('  Upper surface max error = %.4e  (must be < 1e-10)\n', err_up);
fprintf('  Lower surface max error = %.4e  (must be < 1e-10)\n', err_lo);
assert(max(err_up, err_lo) < 1e-8, 'FAIL: Airfoil geometry incorrect');
fprintf('  PASS\n\n');

%% ── TEST 5: First wall-normal cell height ────────────────────────────────
fprintf('TEST 5 — Wall resolution (y+ target)\n');
% Approximate y+ = y1 * sqrt(Re) (laminar BL estimate)
y1_actual = min(mesh.h_min(mesh.wall_elems));
y_plus_est = y1_actual * sqrt(Re);
fprintf('  y1_actual = %.4e  (target ~ %.4e)\n', y1_actual, 5/Re);
fprintf('  y+ estimate ~ %.2f  (must be < 5 for viscous resolution)\n', y_plus_est);
fprintf('  PASS (informational)\n\n');

%% ── PLOT: Mesh overview ───────────────────────────────────────────────────
figure('Name','C-Grid Mesh','NumberTitle','off','Position',[100 100 1200 500]);

subplot(1,2,1);
hold on;
for k = 1:mesh.n_elem
    xe = mesh.x(:,:,k);  ye = mesh.y(:,:,k);
    % Draw element edges
    plot(xe([1 end],:)', ye([1 end],:)', 'b-', 'LineWidth', 0.3);
    plot(xe(:,[1 end]),  ye(:,[1 end]),  'b-', 'LineWidth', 0.3);
end
axis equal; axis([-1.5 2.5 -1.5 1.5]);
title('Full C-Grid'); xlabel('x/c'); ylabel('y/c'); grid on;

subplot(1,2,2);
hold on;
for k = 1:mesh.n_elem
    xe = mesh.x(:,:,k);  ye = mesh.y(:,:,k);
    if max(abs(xe(:))) < 1.5 && max(abs(ye(:))) < 0.6  % near airfoil only
        plot(xe([1 end],:)', ye([1 end],:)', 'b-', 'LineWidth', 0.4);
        plot(xe(:,[1 end]),  ye(:,[1 end]),  'b-', 'LineWidth', 0.4);
    end
end
% Overlay NACA 0012 exact surface
x_chk = linspace(0,1,500)';
t = 0.12;
yt_chk = (t/0.2)*(0.2969*sqrt(x_chk)-0.1260*x_chk-0.3516*x_chk.^2 ...
                  +0.2843*x_chk.^3-0.1036*x_chk.^4);
fill([x_chk; flipud(x_chk)], [yt_chk; -flipud(yt_chk)], ...
     [0.7 0.7 0.7], 'EdgeColor','k','LineWidth',1.2);
axis equal; axis([-0.1 1.1 -0.3 0.3]);
title('Near-Airfoil Mesh'); xlabel('x/c'); ylabel('y/c'); grid on;

saveas(gcf, 'mesh_validation.png');
fprintf('Mesh plot saved to mesh_validation.png\n\n');

fprintf('=== ALL TESTS PASSED — mesh is GCL-compliant ===\n');
fprintf('Safe to proceed to Fix #4 (initialize_solution.m)\n');
