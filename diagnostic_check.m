%% Diagnostic — Initialization and Timestep Check
clear; close all; clc;

Re = 1e4; Mach = 0.3; alpha = 0.0;
gamma = 1.4; Pr = 0.72; R = 287.05;

T_inf   = 288.15; p_inf = 101325;
rho_inf = p_inf/(R*T_inf);
c_inf   = sqrt(gamma*R*T_inf);
U_inf   = Mach*c_inf;
chord   = 1.0;
mu      = rho_inf*U_inf*chord/Re;

fprintf('rho=%.4f kg/m^3  U=%.2f m/s  Re=%.2e  mu=%.4e Pa·s\n', ...
        rho_inf, U_inf, Re, mu);

N = 3; nelem_airfoil=40; nelem_wake=20; nelem_radial=15;
mesh = generate_mesh_NACA0012(nelem_airfoil, nelem_wake, nelem_radial, chord);
U    = initialize_solution(mesh, N, rho_inf, U_inf, alpha, p_inf, gamma);

fprintf('NaN: %d  Inf: %d\n', any(isnan(U(:))), any(isinf(U(:))));
fprintf('rho : [%.3e, %.3e]\n', min(U(:,:,:,1),[],'all'), max(U(:,:,:,1),[],'all'));
fprintf('rhoE: [%.3e, %.3e]\n', min(U(:,:,:,4),[],'all'), max(U(:,:,:,4),[],'all'));

dt = compute_timestep(U, mesh, N, 0.1, gamma);
fprintf('dt = %.6e s\n', dt);
fprintf('Estimated iters (~10 flow-through times): %.0f\n', 10*(chord/U_inf)/dt);
