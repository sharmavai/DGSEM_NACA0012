% =========================================================================
%  roe_flux.m
%  Roe approximate Riemann solver for 2D compressible Euler/NS equations
%
%  FIX #6 Applied:
%    - Roe averages: H_roe was missing sqrt-density weighting
%    - c_roe: guarded against negative argument (can occur near shocks)
%    - Right eigenvector matrix R_mat: previous version had wrong columns
%      for 2D — tangential velocity rows were incorrect
%    - Wave strengths alpha: used ad-hoc formula; replaced with correct
%      Roe (1981) decomposition using pressure and normal-velocity jumps
%    - Harten-Hyman entropy fix: applied to lambda_1 and lambda_4 only
%      (acoustic waves); lambda_2,3 (shear/entropy) left unchanged
%    - Normal flux F_L, F_R computed via local helper (not assumed available)
%
%  References:
%    [1] Roe (1981) J. Comput. Phys. 43 — original Roe scheme
%    [2] Toro (2009) Ch. 11 — 2D Roe eigensystem, wave strengths
%    [3] Harten & Hyman (1983) — entropy fix
%    [4] Blazek (2015) Sec 5.3 — implementation details
% =========================================================================

function F_roe = roe_flux(QL, QR, nx, ny, gamma)
% Inputs:
%   QL, QR  — (4x1) left/right conservative states  [rho,rhou,rhov,E]
%   nx, ny  — outward unit normal FROM left element
%   gamma   — ratio of specific heats
%
% Output:
%   F_roe   — (4x1) numerical flux in normal direction

    %% ── 1. Primitive variables ───────────────────────────────────────────
    rhoL = QL(1);  uL = QL(2)/rhoL;  vL = QL(3)/rhoL;
    EL   = QL(4);
    pL   = (gamma-1)*(EL - 0.5*rhoL*(uL^2+vL^2));
    pL   = max(pL, 1e-12);
    HL   = (EL + pL) / rhoL;
    VnL  = uL*nx + vL*ny;

    rhoR = QR(1);  uR = QR(2)/rhoR;  vR = QR(3)/rhoR;
    ER   = QR(4);
    pR   = (gamma-1)*(ER - 0.5*rhoR*(uR^2+vR^2));
    pR   = max(pR, 1e-12);
    HR   = (ER + pR) / rhoR;
    VnR  = uR*nx + vR*ny;

    %% ── 2. Roe averages (Roe 1981, eq. 5) ───────────────────────────────
    sqL  = sqrt(rhoL);
    sqR  = sqrt(rhoR);
    denom = sqL + sqR;

    u_r  = (sqL*uL + sqR*uR) / denom;
    v_r  = (sqL*vL + sqR*vR) / denom;
    H_r  = (sqL*HL + sqR*HR) / denom;   % FIX: was missing sqrt weighting
    Vn_r = u_r*nx + v_r*ny;

    c2_r = (gamma-1)*(H_r - 0.5*(u_r^2 + v_r^2));
    c2_r = max(c2_r, 1e-10);            % guard against sqrt(negative)
    c_r  = sqrt(c2_r);

    %% ── 3. Eigenvalues ───────────────────────────────────────────────────
    lam1 = Vn_r - c_r;
    lam2 = Vn_r;
    lam3 = Vn_r;
    lam4 = Vn_r + c_r;

    %% ── 4. Harten-Hyman entropy fix (acoustic waves only) ───────────────
    % Prevents expansion shocks at sonic transitions
    % Fix applied to lambda_1 (Vn-c) and lambda_4 (Vn+c)
    cL = sqrt(max((gamma-1)*(HL - 0.5*(uL^2+vL^2)), 1e-10));
    cR = sqrt(max((gamma-1)*(HR - 0.5*(uR^2+vR^2)), 1e-10));

    eps1 = max(0.0, (Vn_r - c_r) - (VnL - cL));
    eps4 = max(0.0, (VnR + cR)   - (Vn_r + c_r));

    if abs(lam1) < eps1 && eps1 > 0
        lam1 = 0.5*(lam1^2/eps1 + eps1);
    end
    if abs(lam4) < eps4 && eps4 > 0
        lam4 = 0.5*(lam4^2/eps4 + eps4);
    end

    %% ── 5. Wave strengths (Toro 2009 eq. 11.75) ─────────────────────────
    % Jump in conserved variables
    drho = rhoR - rhoL;
    du   = uR   - uL;
    dv   = vR   - vL;
    dp   = pR   - pL;
    dVn  = VnR  - VnL;

    % Tangential unit vector: t = (-ny, nx)
    tx = -ny;  ty = nx;
    dVt = du*tx + dv*ty;   % jump in tangential velocity

    % alpha_k = wave strength for k-th characteristic field
    % FIX: wave strengths must use Roe-averaged density (sqL*sqR),
    %       not rhoR.  Previous code had rhoR in a1/a4 — violates Roe
    %       property P (exact resolution of isolated discontinuities).
    rho_r = sqL * sqR;                          % Roe-averaged density
    a1 = (dp - rho_r*c_r*dVn) / (2.0*c_r^2);   % left-running acoustic
    a4 = (dp + rho_r*c_r*dVn) / (2.0*c_r^2);   % right-running acoustic
    a2 = drho - dp/c_r^2;                       % entropy wave
    a3 = rho_r * dVt;                           % shear wave

    %% ── 6. Right eigenvectors (columns, Toro 2009 Table 11.1) ───────────
    % R = [r1 | r2 | r3 | r4]
    % r1: left-running acoustic  (lam = Vn - c)
    r1 = [1;
          u_r - c_r*nx;
          v_r - c_r*ny;
          H_r - Vn_r*c_r];

    % r2: entropy wave  (lam = Vn)
    r2 = [1;
          u_r;
          v_r;
          0.5*(u_r^2 + v_r^2)];

    % r3: shear wave  (lam = Vn)
    r3 = [0;
          tx;
          ty;
          u_r*tx + v_r*ty];

    % r4: right-running acoustic  (lam = Vn + c)
    r4 = [1;
          u_r + c_r*nx;
          v_r + c_r*ny;
          H_r + Vn_r*c_r];

    %% ── 7. Roe flux assembly ─────────────────────────────────────────────
    % F_roe = 0.5*(F_L + F_R) - 0.5 * sum_k |lam_k| * a_k * r_k
    FL = normal_flux(QL, nx, ny, gamma);
    FR = normal_flux(QR, nx, ny, gamma);

    dissipation = abs(lam1)*a1*r1 + abs(lam2)*a2*r2 + ...
                  abs(lam3)*a3*r3 + abs(lam4)*a4*r4;

    F_roe = 0.5*(FL + FR) - 0.5*dissipation;
end

% ── Normal inviscid flux  F*nx + G*ny ────────────────────────────────────
function Fn = normal_flux(Q, nx, ny, gamma)
    rho = Q(1);  u = Q(2)/rho;  v = Q(3)/rho;
    E   = Q(4);
    p   = max((gamma-1)*(E - 0.5*rho*(u^2+v^2)), 1e-12);
    Vn  = u*nx + v*ny;
    Fn  = [rho*Vn;
           rho*u*Vn + p*nx;
           rho*v*Vn + p*ny;
           (E+p)*Vn];
end