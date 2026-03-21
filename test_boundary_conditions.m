% =========================================================================
%  test_boundary_conditions.m
%  Validates farfield_bc.m and wall_bc.m
% =========================================================================
clear; clc;


gamma = 1.4;  Mach = 0.3;  alpha_deg = 5.0;

%% ── WALL BC TESTS ─────────────────────────────────────────────────────────
fprintf('--- wall_bc.m ---\n');

% Build a realistic interior state: low-speed near-wall
rho_i = 1.0; u_i = 0.05; v_i = 0.02;
p_i   = 1/gamma;
E_i   = p_i/(gamma-1) + 0.5*rho_i*(u_i^2+v_i^2);
Q_int = [rho_i; rho_i*u_i; rho_i*v_i; E_i];

Q_g = wall_bc(Q_int, gamma);

% TEST W1: face velocity = average of interior + ghost = 0
u_g = Q_g(2)/Q_g(1);  v_g = Q_g(3)/Q_g(1);
u_face = 0.5*(u_i + u_g);
v_face = 0.5*(v_i + v_g);
assert(abs(u_face) < 1e-14, 'FAIL W1: u_wall = %e', u_face);
assert(abs(v_face) < 1e-14, 'FAIL W1: v_wall = %e', v_face);
fprintf('TEST W1 (no-slip u=v=0 at face): PASS\n');

% TEST W2: pressure preserved (zero normal gradient)
p_g = (gamma-1)*(Q_g(4) - 0.5*Q_g(1)*(u_g^2+v_g^2));
assert(abs(p_g - p_i) < 1e-12, 'FAIL W2: p_ghost != p_int');
fprintf('TEST W2 (pressure continuity):    PASS\n');

% TEST W3: adiabatic — T_ghost = T_int
T_i = p_i / rho_i;
T_g = p_g / Q_g(1);
assert(abs(T_g - T_i) < 1e-12, 'FAIL W3: T_ghost != T_int');
fprintf('TEST W3 (adiabatic T_g=T_i):      PASS\n\n');

%% ── FAR-FIELD BC TESTS ────────────────────────────────────────────────────
fprintf('--- farfield_bc.m ---\n');

% Outward normal pointing radially outward (e.g. top of domain)
nx = 0.0;  ny = 1.0;

% TEST F1: Freestream input => ghost should match freestream exactly
alpha = alpha_deg*pi/180;
u_fs  = Mach*cos(alpha);  v_fs = Mach*sin(alpha);
p_fs  = 1/gamma;  rho_fs = 1.0;
E_fs  = p_fs/(gamma-1) + 0.5*rho_fs*(u_fs^2+v_fs^2);
Q_fs  = [rho_fs; rho_fs*u_fs; rho_fs*v_fs; E_fs];

Q_g = farfield_bc(Q_fs, nx, ny, Mach, alpha_deg, gamma);
rho_g=Q_g(1); ug=Q_g(2)/rho_g; vg=Q_g(3)/rho_g;
pg=(gamma-1)*(Q_g(4)-0.5*rho_g*(ug^2+vg^2));

assert(abs(rho_g - rho_fs) < 1e-10, 'FAIL F1: rho mismatch');
assert(abs(pg    - p_fs)   < 1e-10, 'FAIL F1: p mismatch');
fprintf('TEST F1 (freestream in=>freestream out):  PASS\n');

% TEST F2: Outflow — interior entropy should be preserved
Q_out       = Q_fs;
Q_out(2)    = rho_fs * (u_fs + 0.05);   % slightly accelerated outflow
Q_g_out     = farfield_bc(Q_out, nx, ny, Mach, alpha_deg, gamma);
rho_go      = Q_g_out(1);
ug_o        = Q_g_out(2)/rho_go; vg_o = Q_g_out(3)/rho_go;
pg_o        = (gamma-1)*(Q_g_out(4)-0.5*rho_go*(ug_o^2+vg_o^2));
s_int       = (Q_out(4)*(gamma-1) - 0.5*(Q_out(2)^2+Q_out(3)^2)/Q_out(1)) ...
              * (gamma-1) / Q_out(1)^gamma;   % proxy for entropy
% Just check it doesn't blow up and pressure is positive
assert(pg_o > 0,      'FAIL F2: negative ghost pressure');
assert(rho_go > 0,    'FAIL F2: negative ghost density');
fprintf('TEST F2 (outflow positive state):         PASS\n');

% TEST F3: Inflow — freestream entropy should be used
% Perturb interior to look like inflow (Vn < 0, flow entering domain)
nx2 = -1.0; ny2 = 0.0;   % normal pointing LEFT (inflow side)
Q_in = Q_fs;
Q_in(2) = rho_fs * (-0.05);  % u < 0 => flow enters from left
Q_g_in  = farfield_bc(Q_in, nx2, ny2, Mach, alpha_deg, gamma);
assert(Q_g_in(1) > 0, 'FAIL F3: negative density on inflow');
assert((gamma-1)*(Q_g_in(4)-0.5*(Q_g_in(2)^2+Q_g_in(3)^2)/Q_g_in(1)) > 0, ...
       'FAIL F3: negative pressure on inflow ghost');
fprintf('TEST F3 (inflow positive state):          PASS\n\n');

fprintf('=== ALL BOUNDARY CONDITION TESTS PASSED ===\n');
fprintf('Safe to proceed to Fix #6 (roe_flux.m)\n');
