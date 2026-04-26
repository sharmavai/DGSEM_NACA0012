% =========================================================================
%  main_NACA0012.m
%  DGSEM Solver — 2D Compressible Navier-Stokes, NACA 0012 Airfoil
%
%  FIX #1 Applied:
%    - CFL corrected from 0.2 to 0.10 (RK4+DGSEM N=3 stability limit ~0.143)
%    - Added separate viscous CFL parameter CFL_v
%    - Simulation time set to 30 convective units (Re=1e4 is unsteady;
%      needs ~20c/V transient + 10c/V averaging window per Ferrer et al. 2021)
%    - Convergence monitoring with time-averaged CL, CD, CM
%    - Output every 50 steps for diagnostics
%
%  References:
%    [1] Kopriva (2009) Implementing Spectral Methods — Ch. 8 (CFL limits)
%    [2] Ferrer et al. (2021) DGSEM iLES NACA0012 Re=1e4, Springer NNFM 143
%    [3] Hesthaven & Warburton (2008) Nodal DG Methods — Ch. 6 (RK stability)
%
%  Stability note (from [1], Ch. 8):
%    CFL_max(RK4, DGSEM) = 1/(2N+1)  =>  N=3: CFL_max ~ 0.143
%    Using CFL = 0.10 provides a safe margin and handles mesh skewness.
% =========================================================================

clear; clc; close all;

%% ── 1. PHYSICAL PARAMETERS ──────────────────────────────────────────────
Re    = 1e4;          % Reynolds number  (Re=5000 for Swanson steady benchmark)
Mach  = 0.3;          % Freestream Mach number
alpha = 5.0;          % Angle of attack [degrees]  (0 for symmetric cases)
gamma = 1.4;          % Ratio of specific heats
Pr    = 0.72;         % Prandtl number (air)

% Non-dimensional viscosity (mu_inf = 1/Re in non-dim formulation)
mu    = 1.0 / Re;
% Thermal conductivity from Pr: k = mu * cp / Pr = mu*gamma/(Pr*(gamma-1))
k_cond = mu * gamma / (Pr * (gamma - 1));

%% ── 2. NUMERICAL PARAMETERS ──────────────────────────────────────────────
N     = 3;            % Polynomial degree (P4 elements; 4x4 nodes per element)

% ─── FIX #1: CFL CORRECTION ──────────────────────────────────────────────
% PREVIOUS (WRONG):  CFL = 0.2   <-- exceeds stability limit, causes blow-up
% CORRECTED:
CFL   = 0.10;         % Inviscid CFL   [stability limit = 1/(2*N+1) = 0.143]
CFL_v = 0.10;         % Viscous  CFL   (Fourier number; often more restrictive)
% Both are passed to compute_timestep.m which takes min(dt_inv, dt_visc)
% ─────────────────────────────────────────────────────────────────────────

%% ── 3. MESH PARAMETERS ───────────────────────────────────────────────────
nelem_airfoil = 48;   % Elements wrapping airfoil  (min 48 for Re=1e4 BL)
nelem_wake    = 16;   % Elements in wake extension
nelem_radial  = 20;   % Elements in wall-normal direction
R_farfield    = 30.0; % Far-field radius [chord lengths]  (>=30c recommended)

% First wall-normal cell height for Re=1e4:
%   y1 ~ 5*c/Re * 0.2 ~ 1e-4  (ensures y+ < 1 at wall)
y1_wall = 5.0 / Re;   % passed to mesh generator for radial stretching

%% ── 4. TIME INTEGRATION PARAMETERS ──────────────────────────────────────
% Re=1e4 flow is UNSTEADY (laminar vortex shedding, St ~ 0.15)
% Shedding period T_s ~ c/V / St ~ 1.0/0.15 ~ 6.7 convective times
% Strategy (Ferrer et al. 2021):
%   - Run 20 convective times to kill transient
%   - Time-average over final 10 convective times (> 1 shedding cycle)
t_convective  = 1.0 / Mach;   % One convective time = c/V_inf  [non-dim]
t_transient   = 20.0 * t_convective;   % Transient washout period
t_average     = 10.0 * t_convective;   % Averaging window
t_end         = t_transient + t_average;  % Total simulation time

% NOTE: For the Swanson steady benchmark (Re=5000, alpha=0), the flow IS
% steady. Use t_end = 50 convective times and stop when residual < 1e-9.
% Uncomment for steady case:
% t_end        = 50.0 * t_convective;
% res_tol      = 1e-9;   % steady convergence tolerance

res_tol       = 1e-10;  % residual tolerance for early termination (steady only)

%% ── 5. OUTPUT PARAMETERS ─────────────────────────────────────────────────
output_interval  = 50;    % Print monitor every N steps
save_interval    = 500;   % Save solution snapshot every N steps
fig_interval     = 200;   % Update live convergence plot every N steps

%% ── 6. MESH GENERATION ───────────────────────────────────────────────────
fprintf('Generating C-grid mesh...\n');
fprintf('  Airfoil elements : %d\n', nelem_airfoil);
fprintf('  Radial  elements : %d\n', nelem_radial);
fprintf('  Far-field radius : %.1f c\n', R_farfield);

mesh = generate_mesh_NACA0012(nelem_airfoil, nelem_wake, nelem_radial, ...
                               R_farfield, y1_wall, N);

% Store freestream parameters in mesh (needed by rhs_DGSEM for BCs)
mesh.Mach      = Mach;
mesh.alpha_deg = alpha;

% ── Mandatory free-stream preservation (GCL) check ──────────────────────
fprintf('Running GCL / free-stream preservation test...\n');
Q_fs  = initialize_solution(mesh, Mach, 0.0, gamma, N);  % alpha=0 for GCL test
Q_gcltest = Q_fs;
dt_gcltest = compute_timestep(Q_gcltest, mesh, CFL, CFL_v, gamma, mu, N);
for gclstep = 1:5
    Q_gcltest = rk4_step(Q_gcltest, dt_gcltest, mesh, gamma, mu, k_cond);
end
gcl_residual = max(abs(Q_gcltest(:) - Q_fs(:)));
fprintf('  GCL residual (should be < 1e-10): %e\n', gcl_residual);
if gcl_residual > 1e-8
    warning('GCL FAILED! Metric computation may be incorrect. Check generate_mesh_NACA0012.m');
end

%% ── 7. INITIALIZATION ────────────────────────────────────────────────────
fprintf('\nInitializing solution: M=%.2f, Re=%.0e, alpha=%.1f deg\n', ...
        Mach, Re, alpha);
Q = initialize_solution(mesh, Mach, alpha, gamma, N);

%% ── 8. CONVERGENCE HISTORY STORAGE ──────────────────────────────────────
% Pre-allocate (generous upper bound)
max_steps     = ceil(t_end / (CFL * min(mesh.h_min(:)) / (1 + Mach))) + 1000;
conv_hist     = zeros(max_steps, 8);
% Columns: [step, time, CL, CD, CM, Res_rho, Res_rhoE, dt]

% Time-averaging accumulators
CL_avg = 0;  CD_avg = 0;  CM_avg = 0;
n_avg  = 0;

%% ── 9. MAIN TIME-MARCH LOOP ──────────────────────────────────────────────
fprintf('\nStarting time march: t_end = %.2f (%.0f convective times)\n', ...
        t_end, t_end / t_convective);
fprintf('Averaging starts at t = %.2f (after transient)\n\n', t_transient);
fprintf('%8s %10s %10s %10s %10s %12s\n', ...
        'Step', 'Time', 'CL', 'CD', 'CM', 'Residual');
fprintf('%s\n', repmat('-', 1, 65));

time    = 0.0;
step    = 0;
Q_prev  = Q;

while time < t_end
    step = step + 1;

    % ── Time step (inviscid + viscous limits) ──────────────────────────
    dt = compute_timestep(Q, mesh, CFL, CFL_v, gamma, mu, N);
    dt = min(dt, t_end - time);   % don't overshoot end time

    % ── RK4 update ────────────────────────────────────────────────────
    Q_prev = Q;
    Q = rk4_step(Q, dt, mesh, gamma, mu, k_cond);

    time = time + dt;

    % ── Residual (L-inf norm of dQ/dt) ────────────────────────────────
    % Q is (n_elem, Np1, Np1, 4) — use (:) to get scalar residual
    res_rho  = max(abs(Q(:,:,:,1) - Q_prev(:,:,:,1)), [], 'all') / dt;
    res_rhoE = max(abs(Q(:,:,:,4) - Q_prev(:,:,:,4)), [], 'all') / dt;
    residual = max(res_rho, res_rhoE);

    % ── Check for NaN / blow-up ────────────────────────────────────────
    if ~isfinite(residual) || ~isfinite(dt)
        fprintf('\n*** SOLVER DIVERGED at step=%d, t=%.4f ***\n', step, time);
        fprintf('    Check: CFL=%.3f, dt=%.3e, residual=%.3e\n', CFL, dt, residual);
        fprintf('    Likely causes: metric errors, incorrect BC, or bad mesh.\n');
        break;
    end

    % ── Aerodynamic coefficients ───────────────────────────────────────
    if mod(step, output_interval) == 0
        [CL, CD, CM] = compute_aero_coefficients(Q, mesh, Mach, alpha, gamma, Re, N);

        conv_hist(step, :) = [step, time, CL, CD, CM, res_rho, res_rhoE, dt];

        fprintf('%8d %10.4f %10.5f %10.6f %10.5f %12.4e\n', ...
                step, time, CL, CD, CM, residual);

        % ── Time averaging (post-transient only) ──────────────────────
        if time > t_transient
            CL_avg = CL_avg + CL;
            CD_avg = CD_avg + CD;
            CM_avg = CM_avg + CM;
            n_avg  = n_avg  + 1;
        end

        % ── Early exit for steady cases ────────────────────────────────
        if residual < res_tol && time > 5.0 * t_convective
            fprintf('\nSteady convergence reached (res = %.3e < %.3e)\n', ...
                    residual, res_tol);
            break;
        end
    end

    % ── Save snapshot ──────────────────────────────────────────────────
    if mod(step, save_interval) == 0
        fname = sprintf('snapshot_step%06d.mat', step);
        save(fname, 'Q', 'time', 'step', 'mesh', 'Mach', 'Re', 'alpha');
    end
end

%% ── 10. POST-PROCESSING ──────────────────────────────────────────────────
fprintf('\n%s\n', repmat('=', 1, 65));
fprintf('SIMULATION COMPLETE:  steps=%d,  t=%.4f\n', step, time);

if n_avg > 0
    CL_mean = CL_avg / n_avg;
    CD_mean = CD_avg / n_avg;
    CM_mean = CM_avg / n_avg;
    fprintf('\nTIME-AVERAGED RESULTS (over %.1f convective times, N=%d samples):\n', ...
            t_average / t_convective, n_avg);
    fprintf('  <CL> = %+.5f\n', CL_mean);
    fprintf('  <CD> =  %.5f\n', CD_mean);
    fprintf('  <CM> = %+.5f\n', CM_mean);

    % Validation targets (Ferrer et al. 2021, M=0.3, Re=1e4):
    fprintf('\nValidation targets (Ferrer 2021, M=0.3, Re=1e4):\n');
    fprintf('  alpha=0:  <CL>~0.000,  <CD>~0.040-0.050\n');
    fprintf('  alpha=5:  <CL>~0.55-0.65, <CD>~0.045-0.060\n');
    fprintf('  alpha=10: <CL>~0.95-1.10, <CD>~0.090-0.130\n');
else
    fprintf('  (No time averaging — simulation ended before t_transient)\n');
end

%% ── 11. SAVE RESULTS & CONVERGENCE HISTORY ───────────────────────────────
conv_hist = conv_hist(1:step, :);  % trim unused rows
result_file = sprintf('result_Re%.0e_M%.2f_a%.1f.mat', Re, Mach, alpha);
if n_avg > 0
    save(result_file, 'Q', 'mesh', 'conv_hist', ...
         'CL_mean', 'CD_mean', 'CM_mean', ...
         'Re', 'Mach', 'alpha', 'N', 'CFL', 'time');
else
    save(result_file, 'Q', 'mesh', 'conv_hist', ...
         'Re', 'Mach', 'alpha', 'N', 'CFL', 'time');
end
fprintf('\nResults saved to: %s\n', result_file);

%% ── 12. CONVERGENCE PLOT ─────────────────────────────────────────────────
idx_valid = find(conv_hist(:,1) > 0);
if numel(idx_valid) > 2
    t_plot  = conv_hist(idx_valid, 2);
    CL_plot = conv_hist(idx_valid, 3);
    CD_plot = conv_hist(idx_valid, 4);
    res_plt = max(conv_hist(idx_valid, 6), conv_hist(idx_valid, 7));

    figure('Name','Convergence History','NumberTitle','off');

    subplot(3,1,1);
    plot(t_plot, CL_plot, 'b-', 'LineWidth', 1.2); hold on;
    if n_avg > 0
        yline(CL_mean, 'r--', sprintf('<CL>=%.4f', CL_mean), 'LineWidth', 1.2);
        xline(t_transient, 'k:', 'Averaging start');
    end
    xlabel('Convective time  t \cdot V_\infty/c'); ylabel('C_L');
    title(sprintf('DGSEM NACA0012: Re=%.0e, M=%.2f, \\alpha=%.1f°, N=%d', ...
                  Re, Mach, alpha, N));
    grid on;

    subplot(3,1,2);
    plot(t_plot, CD_plot, 'r-', 'LineWidth', 1.2); hold on;
    if n_avg > 0
        yline(CD_mean, 'k--', sprintf('<CD>=%.5f', CD_mean), 'LineWidth', 1.2);
    end
    xlabel('Convective time'); ylabel('C_D'); grid on;

    subplot(3,1,3);
    semilogy(t_plot, res_plt, 'k-', 'LineWidth', 1.2);
    xlabel('Convective time'); ylabel('Residual ||\deltaQ/\deltat||_\infty');
    grid on;

    saveas(gcf, sprintf('convergence_Re%.0e_M%.2f_a%.1f.png', Re, Mach, alpha));
end

fprintf('\nDone.\n');