%% DGSEM Solver — NACA 0012 Viscous Flow (non-dimensional)
clear; close all; clc;

% Non-dimensional reference state: rho_inf=1, chord=1
% c_inf = sqrt(gamma*p_inf/rho_inf) = 1  =>  U_inf = Mach
Re    = 1e4;
Mach  = 0.3;
alpha = 0.0;       % degrees
gamma = 1.4;
Pr    = 0.72;
R     = 1.0;

rho_inf = 1.0;
p_inf   = 1.0 / gamma;   % gives c_inf = 1
U_inf   = Mach;           % Ma = U/c = Mach/1
T_inf   = 1.0;
chord   = 1.0;
mu      = 1.0 / Re;       % Re = rho*U*c/mu = 1

N             = 3;
nelem_airfoil = 24;
nelem_wake    = 12;
nelem_radial  = 10;
CFL           = 0.1;
max_iter      = 100000;
tol           = 1e-10;
out_freq      = 100;

fprintf('Re=%.2e  M=%.2f  alpha=%.1f  N=%d\n', Re, Mach, alpha, N);

[xi, w] = LGL_quadrature(N);
D       = derivative_matrix(xi, N);
mesh    = generate_mesh_NACA0012(nelem_airfoil, nelem_wake, nelem_radial, chord);
U       = initialize_solution(mesh, N, rho_inf, U_inf, alpha, p_inf, gamma);

fprintf('Elements: %d\n', mesh.nelem);

alpha_rad = alpha * pi / 180;
u_init = U_inf * cos(alpha_rad);
v_init = U_inf * sin(alpha_rad);
E_init = p_inf / ((gamma-1)*rho_inf) + 0.5*(u_init^2 + v_init^2);

res_hist = zeros(max_iter, 1);
time     = 0;

fprintf('\n%-8s %-10s %-12s %-10s %-10s %-10s\n', ...
        'Iter','Time','Residual','CL','CD','CM');
fprintf('%s\n', repmat('-',1,62));

for iter = 1:max_iter
    U_old = U;
    dt    = compute_timestep(U, mesh, N, CFL, gamma);

    rhs = @(Ui) rhs_DGSEM(Ui, mesh, N, D, w, gamma, mu, Pr, R, ...
                           U_inf, alpha, p_inf, T_inf);
    k1 = rhs(U);
    k2 = rhs(U + 0.5*dt*k1);
    k3 = rhs(U + 0.5*dt*k2);
    k4 = rhs(U +     dt*k3);
    U  = U + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);

    % NaN guard: reset any bad element to free-stream
    bad = squeeze(any(any(any(~isfinite(U), 4), 3), 2));
    if any(bad)
        fprintf('WARNING: NaN in %d elements at iter %d\n', sum(bad), iter);
        U(bad,:,:,1) = rho_inf;
        U(bad,:,:,2) = rho_inf * u_init;
        U(bad,:,:,3) = rho_inf * v_init;
        U(bad,:,:,4) = rho_inf * E_init;
    end

    dU  = U(:,:,:,1) - U_old(:,:,:,1);
    res = sqrt(sum(dU(:).^2) / numel(dU)) / dt;
    res_hist(iter) = res;
    time = time + dt;

    if mod(iter, out_freq) == 0 || iter == 1
        [CL,CD,CM] = compute_aero_coefficients(U, mesh, N, gamma, ...
                                                rho_inf, U_inf, chord);
        fprintf('%-8d %-10.4f %-12.3e %-10.6f %-10.6f %-10.6f\n', ...
                iter, time, res, CL, CD, CM);
    end

    if iter > 200 && isfinite(res) && res < tol
        fprintf('\nConverged at iteration %d\n', iter); break
    end
    if iter > 50 && (~isfinite(res) || res > 1e10)
        fprintf('\nDiverged at iteration %d — reduce CFL\n', iter); break
    end
end

[CL,CD,CM]    = compute_aero_coefficients(U, mesh, N, gamma, rho_inf, U_inf, chord);
[Cp, x_s, ~] = compute_pressure_coefficient(U, mesh, N, gamma, rho_inf, U_inf);

fprintf('\nCL=%.6f  CD=%.6f  CM=%.6f\n', CL, CD, CM);

if ~exist('figures','dir'), mkdir('figures'); end

figure;
semilogy(1:iter, res_hist(1:iter), 'b-', 'LineWidth', 1.5);
grid on; xlabel('Iteration'); ylabel('Residual');
title(sprintf('Convergence  Re=%.0e  M=%.2f', Re, Mach));
saveas(gcf, 'figures/residual.png');

figure;
plot(x_s, Cp, 'b-o', 'MarkerSize', 3, 'LineWidth', 1.5);
set(gca,'YDir','reverse'); grid on; xlim([0 1]);
xlabel('x/c'); ylabel('C_p');
title(sprintf('Pressure Distribution  Re=%.0e', Re));
saveas(gcf, 'figures/Cp.png');