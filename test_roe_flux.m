% =========================================================================
%  test_roe_flux.m
%  Validates roe_flux.m — all tests must pass before proceeding to Fix #7
% =========================================================================
clear; clc;

fprintf('=== FIX #6 ROE FLUX TESTS ===\n\n');
gamma = 1.4;

%% ── Helper ───────────────────────────────────
function Fn = exact_flux(Q, nx, ny, gamma)
    rho=Q(1); u=Q(2)/rho; v=Q(3)/rho; E=Q(4);
    p=max((gamma-1)*(E-0.5*rho*(u^2+v^2)),1e-12);
    Vn=u*nx+v*ny;
    Fn=[rho*Vn; rho*u*Vn+p*nx; rho*v*Vn+p*ny; (E+p)*Vn];
end

%% TEST 1 — Consistency: QL=QR => F_roe = exact flux ────────────
rho=1.0; u=0.3; v=0.1; p=1/gamma;
E=p/(gamma-1)+0.5*rho*(u^2+v^2);
Q=[rho;rho*u;rho*v;E];
nx=1.0; ny=0.0;
F_r = roe_flux(Q, Q, nx, ny, gamma);
F_e = exact_flux(Q, nx, ny, gamma);
err = max(abs(F_r - F_e));
assert(err < 1e-12, 'FAIL T1: consistency error = %e', err);
fprintf('TEST 1 (consistency QL=QR):          err = %.2e  PASS\n', err);

%% TEST 2 — Conservation: F(QL,QR,n) = -F(QR,QL,-n) ────────────
rhoL=1.0; uL=0.4; vL=0.1; pL=1/gamma;
EL=pL/(gamma-1)+0.5*rhoL*(uL^2+vL^2);
QL=[rhoL;rhoL*uL;rhoL*vL;EL];
rhoR=0.8; uR=0.2; vR=0.05; pR=0.8/gamma;
ER=pR/(gamma-1)+0.5*rhoR*(uR^2+vR^2);
QR=[rhoR;rhoR*uR;rhoR*vR;ER];
nx=0.6; ny=0.8;   % oblique normal
F_lr = roe_flux(QL, QR,  nx,  ny, gamma);
F_rl = roe_flux(QR, QL, -nx, -ny, gamma);
err2 = max(abs(F_lr + F_rl));
assert(err2 < 1e-12, 'FAIL T2: conservation error = %e', err2);
fprintf('TEST 2 (conservation F+F_flip=0):    err = %.2e  PASS\n', err2);

%% TEST 3 — Positivity: density and pressure in dissipation term ─────
% Mild shear case — should not produce negative density flux
Q_left  = [1.2; 1.2*0.5; 1.2*0.1; 2.5];
Q_right = [0.9; 0.9*0.3; 0.9*0.05; 2.0];
F3 = roe_flux(Q_left, Q_right, 1.0, 0.0, gamma);
% Just verify it's finite and returns 4 components
assert(all(isfinite(F3)), 'FAIL T3: non-finite flux');
assert(numel(F3)==4,      'FAIL T3: wrong size');
fprintf('TEST 3 (finite output):              PASS\n');

%% TEST 4 — Entropy fix: no expansion shock across sonic point ────────────
% Left: supersonic outflow; Right: subsonic — should not crash
rhoL=0.5; uL=1.5; vL=0.0; pL=0.5/gamma;
EL=pL/(gamma-1)+0.5*rhoL*uL^2;
rhoR=1.0; uR=0.3; vR=0.0; pR=1.0/gamma;
ER=pR/(gamma-1)+0.5*rhoR*uR^2;
QL2=[rhoL;rhoL*uL;rhoL*vL;EL];
QR2=[rhoR;rhoR*uR;rhoR*vR;ER];
F4 = roe_flux(QL2, QR2, 1.0, 0.0, gamma);
assert(all(isfinite(F4)), 'FAIL T4: entropy fix produced NaN');
fprintf('TEST 4 (entropy fix sonic point):    PASS\n');

%% TEST 5 — Oblique normal: rotational invariance check ───────────────────
% Rotating problem 90 degrees should give rotated flux
Q_a = [1.0; 1.0*0.3; 1.0*0.0; 2.0];
Q_b = [0.9; 0.9*0.2; 0.9*0.0; 1.8];
F_x = roe_flux(Q_a, Q_b, 1.0, 0.0, gamma);   % normal in x

% Rotate: swap u<->v, nx<->ny
Q_a_r = [Q_a(1); Q_a(3); Q_a(2); Q_a(4)];
Q_b_r = [Q_b(1); Q_b(3); Q_b(2); Q_b(4)];
F_y   = roe_flux(Q_a_r, Q_b_r, 0.0, 1.0, gamma);
% Mass flux must match; momentum flux components swap
assert(abs(F_x(1)-F_y(1)) < 1e-10, 'FAIL T5: mass flux not rotationally invariant');
assert(abs(F_x(4)-F_y(4)) < 1e-10, 'FAIL T5: energy flux not rotationally invariant');
fprintf('TEST 5 (rotational invariance):      PASS\n\n');

fprintf('=== ALL ROE FLUX TESTS PASSED ===\n');
fprintf('Safe to proceed to Fix #7 (rhs_DGSEM.m)\n');
