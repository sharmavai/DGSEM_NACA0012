% =========================================================================
%  main_NACA0012.m
%  DGSEM Solver — 2D Compressible Navier-Stokes, NACA 0012 Airfoil
%
%  Complete solver driver with:
%    - BR2 viscous scheme for proper Navier-Stokes computation
%    - CFL corrected for RK4+DGSEM N=3 stability (limit ~0.143)
%    - Separate viscous CFL parameter (Fourier number)
%    - Time-averaged CL, CD, CM for unsteady cases
%    - Optional grid convergence study (h-refinement)
%    - GCL free-stream preservation check
%
%  Validation case: M=0.3, Re=1e4, alpha=5° (Ferrer et al. 2021)
%
%  References:
%    [1] Kopriva (2009) Implementing Spectral Methods — Ch. 8
%    [2] Ferrer et al. (2021) DGSEM iLES NACA0012 Re=1e4, Springer NNFM 143
%    [3] Hesthaven & Warburton (2008) Nodal DG Methods — Ch. 6
%    [4] Bassi & Rebay (1997) J. Comput. Phys. — BR2 scheme
%    [5] Gregory & O'Reilly (1970) R&M 3726 — NACA 0012 experimental data
% =========================================================================

clear; clc; close all;

%% ── 1. PHYSICAL PARAMETERS ──────────────────────────────────────────────
Re    = 1e4;          % Reynolds number
Mach  = 0.3;          % Freestream Mach number
alpha = 5.0;          % Angle of attack [degrees]
gamma = 1.4;          % Ratio of specific heats
Pr    = 0.72;         % Prandtl number (air)

% Non-dimensional viscosity (mu_inf = 1/Re)
mu    = 1.0 / Re;
% Thermal conductivity from Pr: k = mu * cp / Pr = mu*gamma/(Pr*(gamma-1))
k_cond = mu * gamma / (Pr * (gamma - 1));

%% ── 2. NUMERICAL PARAMETERS ──────────────────────────────────────────────
N     = 3;            % Polynomial degree (P4 elements; 4x4 nodes per element)
CFL   = 0.08;         % Inviscid CFL   [stability limit = 1/(2*N+1) = 0.143]
CFL_v = 0.08;         % Viscous  CFL   (Fourier number)

%% ── 3. GRID CONVERGENCE STUDY ────────────────────────────────────────────
% Run multiple mesh levels for h-refinement convergence analysis.
% Set do_convergence_study = true for paper-publishable grid study.
do_convergence_study = true;

% Mesh levels: [n_airfoil, n_wake, n_radial]
if do_convergence_study
    mesh_levels = {
        struct('n_airfoil', 24, 'n_wake', 8,  'n_radial', 10, 'label', 'L1 (coarse)');
        struct('n_airfoil', 48, 'n_wake', 16, 'n_radial', 20, 'label', 'L2 (medium)');
        struct('n_airfoil', 72, 'n_wake', 24, 'n_radial', 30, 'label', 'L3 (fine)');
    };
else
    mesh_levels = {
        struct('n_airfoil', 48, 'n_wake', 16, 'n_radial', 20, 'label', 'L2');
    };
end

R_farfield = 30.0;     % Far-field radius [chord lengths]
y1_wall    = 5.0/Re;   % First wall-normal cell height

%% ── 4. TIME INTEGRATION PARAMETERS ────────────────────────────────────────
t_convective  = 1.0 / Mach;   % One convective time = c/V_inf
t_transient   = 20.0 * t_convective;   % Transient washout
t_average     = 10.0 * t_convective;   % Averaging window
t_end         = t_transient + t_average;
res_tol       = 1e-10;

%% ── 5. OUTPUT PARAMETERS ──────────────────────────────────────────────────
output_interval  = 50;
save_interval    = 500;
fig_interval     = 200;

%% ── 6. RESULTS STORAGE ───────────────────────────────────────────────────
results_CL = zeros(numel(mesh_levels), 1);
results_CD = zeros(numel(mesh_levels), 1);
results_CM = zeros(numel(mesh_levels), 1);
results_Nelem = zeros(numel(mesh_levels), 1);
results_h    = zeros(numel(mesh_levels), 1);
results_time = zeros(numel(mesh_levels), 1);

%% ── 7. MAIN LOOP OVER MESH LEVELS ───────────────────────────────────────
for lvl = 1:numel(mesh_levels)
    ml = mesh_levels{lvl};

    fprintf('\n%s\n', repmat('=', 1, 70));
    fprintf('GRID LEVEL %d: %s\n', lvl, ml.label);
    fprintf('  n_airfoil=%d  n_wake=%d  n_radial=%d\n', ...
            ml.n_airfoil, ml.n_wake, ml.n_radial);
    fprintf('%s\n', repmat('=', 1, 70));

    %% ── MESH GENERATION ──────────────────────────────────────────────────
    fprintf('Generating C-grid mesh...\n');
    mesh = generate_mesh_NACA0012(ml.n_airfoil, ml.n_wake, ml.n_radial, ...
                                   R_farfield, y1_wall, N);
    mesh.Mach      = Mach;
    mesh.alpha_deg = alpha;

    results_Nelem(lvl) = mesh.n_elem;
    results_h(lvl)    = mean(mesh.h_min);

    %% ── GCL CHECK ────────────────────────────────────────────────────────
    fprintf('Running GCL / free-stream preservation test...\n');
    Q_fs  = initialize_solution(mesh, Mach, 0.0, gamma, N);
    Q_gcltest = Q_fs;
    dt_gcltest = compute_timestep(Q_gcltest, mesh, CFL, CFL_v, gamma, mu, N);
    for gclstep = 1:5
        Q_gcltest = rk4_step(Q_gcltest, dt_gcltest, mesh, gamma, mu, k_cond);
    end
    gcl_residual = max(abs(Q_gcltest(:) - Q_fs(:)));
    fprintf('  GCL residual: %e (should be < 1e-10)\n', gcl_residual);
    if gcl_residual > 1e-8
        warning('GCL FAILED! Metric computation may be incorrect.');
    end

    %% ── INITIALIZATION ──────────────────────────────────────────────────
    fprintf('\nInitializing: M=%.2f, Re=%.0e, alpha=%.1f deg\n', Mach, Re, alpha);
    Q = initialize_solution(mesh, Mach, alpha, gamma, N);

    %% ── CONVERGENCE HISTORY ──────────────────────────────────────────────
    max_steps = ceil(t_end / (CFL * min(mesh.h_min(:)) / (1 + Mach))) + 1000;
    conv_hist = zeros(max_steps, 8);

    CL_avg = 0;  CD_avg = 0;  CM_avg = 0;
    n_avg  = 0;

    %% ── TIME MARCH ──────────────────────────────────────────────────────
    fprintf('\nStarting time march: t_end = %.2f (%.0f convective times)\n', ...
            t_end, t_end / t_convective);
    fprintf('%8s %10s %10s %10s %10s %12s\n', ...
            'Step', 'Time', 'CL', 'CD', 'CM', 'Residual');
    fprintf('%s\n', repmat('-', 1, 65));

    time    = 0.0;
    step    = 0;
    Q_prev  = Q;
    t_start_loop = tic;

    while time < t_end
        step = step + 1;

        dt = compute_timestep(Q, mesh, CFL, CFL_v, gamma, mu, N);
        dt = min(dt, t_end - time);

        Q_prev = Q;
        Q = rk4_step(Q, dt, mesh, gamma, mu, k_cond);
        time = time + dt;

        res_rho  = max(abs(Q(:,:,:,1) - Q_prev(:,:,:,1)), [], 'all') / dt;
        res_rhoE = max(abs(Q(:,:,:,4) - Q_prev(:,:,:,4)), [], 'all') / dt;
        residual = max(res_rho, res_rhoE);

        if ~isfinite(residual) || ~isfinite(dt)
            fprintf('\n*** SOLVER DIVERGED at step=%d, t=%.4f ***\n', step, time);
            fprintf('    CFL=%.3f, dt=%.3e, residual=%.3e\n', CFL, dt, residual);
            break;
        end

        if mod(step, output_interval) == 0
            [CL, CD, CM] = compute_aero_coefficients(Q, mesh, Mach, alpha, gamma, Re, N);
            conv_hist(step, :) = [step, time, CL, CD, CM, res_rho, res_rhoE, dt];

            fprintf('%8d %10.4f %10.5f %10.6f %10.5f %12.4e\n', ...
                    step, time, CL, CD, CM, residual);

            if time > t_transient
                CL_avg = CL_avg + CL;
                CD_avg = CD_avg + CD;
                CM_avg = CM_avg + CM;
                n_avg  = n_avg  + 1;
            end

            if residual < res_tol && time > 5.0 * t_convective
                fprintf('\nSteady convergence reached (res = %.3e)\n', residual);
                break;
            end
        end

        if mod(step, save_interval) == 0
            fname = sprintf('snapshot_step%06d.mat', step);
            save(fname, 'Q', 'time', 'step', 'mesh', 'Mach', 'Re', 'alpha');
        end
    end

    elapsed = toc(t_start_loop);
    fprintf('\nWall time: %.1f seconds\n', elapsed);

    %% ── POST-PROCESSING ─────────────────────────────────────────────────
    fprintf('\n%s\n', repmat('=', 1, 65));
    fprintf('SIMULATION COMPLETE:  steps=%d,  t=%.4f\n', step, time);

    results_time(lvl) = elapsed;

    if n_avg > 0
        CL_mean = CL_avg / n_avg;
        CD_mean = CD_avg / n_avg;
        CM_mean = CM_avg / n_avg;
        results_CL(lvl) = CL_mean;
        results_CD(lvl) = CD_mean;
        results_CM(lvl) = CM_mean;

        fprintf('\nTIME-AVERAGED RESULTS (%.1f conv. times, N=%d samples):\n', ...
                t_average / t_convective, n_avg);
        fprintf('  <CL> = %+.5f\n', CL_mean);
        fprintf('  <CD> =  %.5f\n', CD_mean);
        fprintf('  <CM> = %+.5f\n', CM_mean);
    else
        fprintf('  (No averaging — simulation ended before transient)\n');
    end

    %% ── SAVE RESULTS ────────────────────────────────────────────────────
    conv_hist = conv_hist(1:step, :);
    result_file = sprintf('result_Re%.0e_M%.2f_a%.1f_L%d.mat', Re, Mach, alpha, lvl);
    if n_avg > 0
        save(result_file, 'Q', 'mesh', 'conv_hist', ...
             'CL_mean', 'CD_mean', 'CM_mean', ...
             'Re', 'Mach', 'alpha', 'N', 'CFL', 'time');
    else
        save(result_file, 'Q', 'mesh', 'conv_hist', ...
             'Re', 'Mach', 'alpha', 'N', 'CFL', 'time');
    end
    fprintf('Results saved to: %s\n', result_file);

    %% ── CONVERGENCE PLOT ──────────────────────────────────────────────────
    idx_valid = find(conv_hist(:,1) > 0);
    if numel(idx_valid) > 2
        t_plot  = conv_hist(idx_valid, 2);
        CL_plot = conv_hist(idx_valid, 3);
        CD_plot = conv_hist(idx_valid, 4);
        res_plt = max(conv_hist(idx_valid, 6), conv_hist(idx_valid, 7));

        figure('Name', sprintf('Level %d Convergence', lvl), ...
               'NumberTitle', 'off', 'Position', [100 100 800 600]);

        subplot(3,1,1);
        plot(t_plot, CL_plot, 'b-', 'LineWidth', 1.2); hold on;
        if n_avg > 0
            yline(CL_mean, 'r--', sprintf('<C_L>=%.4f', CL_mean), 'LineWidth', 1.2);
            xline(t_transient, 'k:', 'Averaging start');
        end
        xlabel('Convective time  t \cdot V_\infty/c'); ylabel('C_L');
        title(sprintf('DGSEM NACA0012: Re=%.0e, M=%.2f, \\alpha=%.1f°, N=%d, %s', ...
                      Re, Mach, alpha, N, ml.label));
        grid on;

        subplot(3,1,2);
        plot(t_plot, CD_plot, 'r-', 'LineWidth', 1.2); hold on;
        if n_avg > 0
            yline(CD_mean, 'k--', sprintf('<C_D>=%.5f', CD_mean), 'LineWidth', 1.2);
        end
        xlabel('Convective time'); ylabel('C_D'); grid on;

        subplot(3,1,3);
        semilogy(t_plot, res_plt, 'k-', 'LineWidth', 1.2);
        xlabel('Convective time'); ylabel('Residual ||\deltaQ/\deltat||_\infty');
        grid on;

        saveas(gcf, sprintf('convergence_Re%.0e_M%.2f_a%.1f_L%d.png', ...
                            Re, Mach, alpha, lvl));
    end
end

%% ── 8. GRID CONVERGENCE ANALYSIS ─────────────────────────────────────────
if do_convergence_study && all(results_CL ~= 0)
    fprintf('\n\n%s\n', repmat('=', 1, 70));
    fprintf('GRID CONVERGENCE STUDY RESULTS\n');
    fprintf('%s\n', repmat('=', 1, 70));

    fprintf('\n%6s %8s %12s %12s %12s %10s\n', ...
            'Level', 'N_elem', '<C_L>', '<C_D>', '<C_M>', 'h_avg');
    fprintf('%s\n', repmat('-', 1, 65));
    for lvl = 1:numel(mesh_levels)
        fprintf('%6d %8d %+12.5f %12.6f %+12.5f %10.3e\n', ...
                lvl, results_Nelem(lvl), results_CL(lvl), results_CD(lvl), ...
                results_CM(lvl), results_h(lvl));
    end

    % Estimate convergence rates (Richardson extrapolation)
    fprintf('\nConvergence rates (estimated):\n');
    for lvl = 2:min(numel(mesh_levels), 3)
        h_ratio = results_h(lvl-1) / results_h(lvl);
        fprintf('  L%d -> L%d: h ratio = %.2f\n', lvl-1, lvl, h_ratio);

        if abs(results_CL(lvl-1)) > 1e-10 && abs(results_CL(lvl)) > 1e-10
            p_CL = log(abs(results_CL(lvl-1) - results_CL(min(lvl+1,numel(mesh_levels)))) / ...
                       abs(results_CL(lvl) - results_CL(min(lvl+1,numel(mesh_levels))))) / ...
                       log(h_ratio);
            fprintf('    C_L order ~ %.2f\n', abs(p_CL));
        end
        if abs(results_CD(lvl-1)) > 1e-10 && abs(results_CD(lvl)) > 1e-10
            p_CD = log(abs(results_CD(lvl-1) - results_CD(min(lvl+1,numel(mesh_levels)))) / ...
                       abs(results_CD(lvl) - results_CD(min(lvl+1,numel(mesh_levels))))) / ...
                       log(h_ratio);
            fprintf('    C_D order ~ %.2f\n', abs(p_CD));
        end
    end

    % Validation targets
    fprintf('\nValidation targets (Ferrer et al. 2021, M=0.3, Re=1e4):\n');
    fprintf('  alpha=0:  <CL>~0.000,  <CD>~0.040-0.050\n');
    fprintf('  alpha=5:  <CL>~0.55-0.65, <CD>~0.045-0.060\n');
    fprintf('  alpha=10: <CL>~0.95-1.10, <CD>~0.090-0.130\n');

    %% ── Grid convergence plot ────────────────────────────────────────────
    figure('Name', 'Grid Convergence', 'NumberTitle', 'off', ...
           'Position', [200 200 1000 400]);

    subplot(1,2,1);
    for lvl = 1:numel(mesh_levels)
        ml = mesh_levels{lvl};
        semilogy(results_h(lvl), abs(results_CL(lvl)), 'bo-', ...
                 'LineWidth', 1.5, 'MarkerSize', 8, 'MarkerFaceColor', 'b');
        hold on;
        text(results_h(lvl)*1.3, abs(results_CL(lvl))*1.2, ml.label, 'FontSize', 9);
    end
    xlabel('Average cell size h'); ylabel('|C_L|'); grid on;
    title('Lift coefficient convergence');

    subplot(1,2,2);
    for lvl = 1:numel(mesh_levels)
        ml = mesh_levels{lvl};
        semilogy(results_h(lvl), abs(results_CD(lvl)), 'rs-', ...
                 'LineWidth', 1.5, 'MarkerSize', 8, 'MarkerFaceColor', 'r');
        hold on;
        text(results_h(lvl)*1.3, abs(results_CD(lvl))*1.2, ml.label, 'FontSize', 9);
    end
    xlabel('Average cell size h'); ylabel('C_D'); grid on;
    title('Drag coefficient convergence');

    saveas(gcf, sprintf('grid_convergence_Re%.0e_M%.2f_a%.1f.png', Re, Mach, alpha));
    fprintf('\nGrid convergence plot saved.\n');
end

%% ── 9. Cp DISTRIBUTION (finest mesh) ─────────────────────────────────────
if do_convergence_study && exist('Q', 'var') && exist('mesh', 'var')
    [Cp, x_surf, y_surf] = compute_pressure_coefficient(Q, mesh, gamma, Mach, N);

    figure('Name', 'Cp Distribution', 'NumberTitle', 'off', ...
           'Position', [300 300 800 400]);
    plot(x_surf, Cp, 'b.-', 'LineWidth', 1.0, 'MarkerSize', 4);
    xlabel('x/c'); ylabel('C_p'); grid on;
    title(sprintf('Surface Cp: Re=%.0e, M=%.2f, \\alpha=%.1f°', Re, Mach, alpha));
    ylim([max(min(Cp)*1.1, -3), max(Cp)*1.1]);
    saveas(gcf, sprintf('Cp_Re%.0e_M%.2f_a%.1f.png', Re, Mach, alpha));
    fprintf('Cp distribution saved.\n');
end

fprintf('\nDone.\n');
