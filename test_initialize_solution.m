% =========================================================================
%  test_initialize_solution.m
%  Validates initialize_solution.m  — run before any flow simulation
%  FIX: Updated indexing for new Q layout (n_elem, Np1, Np1, 4)
% =========================================================================
clear; clc;

gamma = 1.4;  Mach = 0.3;  N = 3;

% Minimal stub mesh
mesh.n_elem = 2;

%% TEST 1 — zero AoA: v_inf must be 0
Q = initialize_solution(mesh, Mach, 0.0, gamma, N);
rho = Q(1,1,1,1);  ru = Q(1,1,1,2);  rv = Q(1,1,1,3);  E = Q(1,1,1,4);
u = ru/rho;  v = rv/rho;
assert(abs(v) < 1e-14,       'FAIL: v_inf != 0 at alpha=0');
assert(abs(u - Mach) < 1e-14,'FAIL: u_inf != Mach at alpha=0');
fprintf('TEST 1 (alpha=0):   u=%.6f  v=%.2e  PASS\n', u, v);

%% TEST 2 — alpha=5 deg: both components non-zero and correct
alpha_deg = 5.0;
Q = initialize_solution(mesh, Mach, alpha_deg, gamma, N);
u2 = Q(1,1,1,2)/Q(1,1,1,1);
v2 = Q(1,1,1,3)/Q(1,1,1,1);
u_ref = Mach*cos(alpha_deg*pi/180);
v_ref = Mach*sin(alpha_deg*pi/180);
assert(abs(u2-u_ref)<1e-12, 'FAIL: u_inf wrong at alpha=5');
assert(abs(v2-v_ref)<1e-12, 'FAIL: v_inf wrong at alpha=5');
fprintf('TEST 2 (alpha=5):   u=%.6f  v=%.6f  PASS\n', u2, v2);

%% TEST 3 — total velocity magnitude = Mach (for all AoA)
for alpha_deg = [0 2 5 10 -5]
    Q  = initialize_solution(mesh, Mach, alpha_deg, gamma, N);
    rho = Q(1,1,1,1); u=Q(1,1,1,2)/rho; v=Q(1,1,1,3)/rho;
    V   = sqrt(u^2+v^2);
    assert(abs(V-Mach)<1e-12,'FAIL: |V| != Mach at alpha=%.1f',alpha_deg);
end
fprintf('TEST 3 (|V|=Mach):  PASS for alpha = 0,2,5,10,-5 deg\n');

%% TEST 4 — energy consistency: recover p from E
Q  = initialize_solution(mesh, Mach, 5.0, gamma, N);
rho=Q(1,1,1,1); u=Q(1,1,1,2)/rho; v=Q(1,1,1,3)/rho; E=Q(1,1,1,4);
p_rec = (gamma-1)*(E - 0.5*rho*(u^2+v^2));
p_ref = 1/gamma;
assert(abs(p_rec-p_ref)<1e-12,'FAIL: pressure recovery error = %e',abs(p_rec-p_ref));
fprintf('TEST 4 (pressure):  p=%.6f  PASS\n', p_rec);

%% TEST 5 — uniform field: every node identical
Q = initialize_solution(mesh, Mach, 5.0, gamma, N);
assert(max(abs(diff(Q(:,:,:,1)),1,'all'))<1e-15,'FAIL: rho not uniform');
assert(max(abs(diff(Q(:,:,:,2)),1,'all'))<1e-15,'FAIL: rho*u not uniform');
fprintf('TEST 5 (uniform):   PASS\n\n');

fprintf('=== ALL TESTS PASSED ===\n');
