% =========================================================================
%  LGL_quadrature.m
%  Legendre-Gauss-Lobatto nodes and weights
%
%  BUGFIX (N=4 weight sum wrong):
%    Initial guess was -cos((j+0.5)*pi/N)  — WRONG for N>=4.
%    For N=4, j=1: guess=-0.383, but nearest LGL node is -0.655.
%    Newton on P''_4 at -0.383 takes a step of ~-10 and diverges/
%    converges to wrong node => wrong P_N value => wrong weight.
%
%    CORRECT (Kopriva 2009, Algorithm 24):
%    Initial guess = -cos(j*pi/N)  (standard Chebyshev nodes)
%    For N=4: j=1:-0.707, j=2:0, j=3:+0.707 — all within 0.05 of true nodes.
%
%  References:
%    [1] Kopriva (2009) Algorithm 24 p.67
% =========================================================================

function [xi, w] = LGL_quadrature(N)

    Np1 = N + 1;

    %% ── 1. Endpoints ─────────────────────────────────────────────────────
    xi      = zeros(Np1, 1);
    xi(1)   = -1.0;
    xi(end) =  1.0;

    %% ── 2. Interior nodes — Newton on P'_N(x) = 0 ───────────────────────
    % FIXED initial guess: -cos(j*pi/N)  (Kopriva Algorithm 24)
    for j = 1:N-1
        x = -cos(j * pi / N);          % <-- KEY FIX (was (j+0.5)*pi/N)
        for iter = 1:100
            [~, ~, dPn, d2Pn] = leg_derivs(x, N);
            delta = -dPn / d2Pn;
            x     =  x + delta;
            if abs(delta) < 1e-15, break; end
        end
        xi(j+1) = x;
    end

    %% ── 3. Weights (Kopriva eq. 3.18) ────────────────────────────────────
    w = zeros(Np1, 1);
    for i = 1:Np1
        [~, Pn, ~, ~] = leg_derivs(xi(i), N);
        w(i) = 2.0 / (N * Np1 * Pn^2);
    end

    assert(abs(sum(w) - 2.0) < 1e-12, ...
           'LGL weights sum = %.15f, expected 2.0  (N=%d)', sum(w), N);
    assert(xi(1) == -1.0 && xi(end) == 1.0, 'Endpoint nodes not exactly +-1');
end

% ── P_{N-1}, P_N, P'_N, P''_N at scalar x ────────────────────────────────
function [Pnm1, Pn, dPn, d2Pn] = leg_derivs(x, N)
    if N==0, Pnm1=0; Pn=1;  dPn=0;  d2Pn=0; return; end
    if N==1, Pnm1=1; Pn=x;  dPn=1;  d2Pn=0; return; end
    Pkm2=1; Pkm1=x; dPkm2=0; dPkm1=1; d2Pkm2=0; d2Pkm1=0;
    for k = 2:N
        Pk   = ((2*k-1)*x*Pkm1   - (k-1)*Pkm2)                   / k;
        dPk  = ((2*k-1)*(Pkm1  + x*dPkm1)  - (k-1)*dPkm2)        / k;
        d2Pk = ((2*k-1)*(2*dPkm1 + x*d2Pkm1) - (k-1)*d2Pkm2)     / k;
        Pkm2=Pkm1; Pkm1=Pk; dPkm2=dPkm1; dPkm1=dPk; d2Pkm2=d2Pkm1; d2Pkm1=d2Pk;
    end
    Pnm1=Pkm2; Pn=Pkm1; dPn=dPkm1; d2Pn=d2Pkm1;
    if abs(d2Pn) < 1e-30, d2Pn = 1e-30; end
end