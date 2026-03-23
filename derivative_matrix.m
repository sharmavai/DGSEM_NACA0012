% =========================================================================
%  derivative_matrix.m
%  DGSEM spectral derivative matrix via direct Lagrange formula
%
%  Uses Kopriva (2009) eq. 3.44 — no Vandermonde inversion needed.
%  Calls LGL_quadrature for [xi, w] only (2-output form).
%
%  SBP:  M*D + (M*D)' = B
%  References: [1] Kopriva (2009) eq. 3.44, 3.55
% =========================================================================

function [D, M, B] = derivative_matrix(N)

    Np1 = N + 1;
    [xi, w] = LGL_quadrature(N);   % 2-output — no Vandermonde

    %% ── 1. Evaluate P_N at all LGL nodes ────────────────────────────────
    PN = zeros(Np1, 1);
    for i = 1:Np1
        [~, Pn, ~, ~] = leg_derivs(xi(i), N);
        PN(i) = Pn;
    end

    %% ── 2. Off-diagonal entries (Kopriva eq. 3.44) ───────────────────────
    D = zeros(Np1);
    for i = 1:Np1
        for j = 1:Np1
            if i ~= j
                D(i,j) = PN(i) / (PN(j) * (xi(i) - xi(j)));
            end
        end
    end

    %% ── 3. Diagonal — negative-sum trick ─────────────────────────────────
    for i = 1:Np1
        D(i,i) = -sum(D(i,:));
    end

    %% ── 4. Exact endpoints (Kopriva eq. 3.55) ────────────────────────────
    D(1,   1)   = -N*(N+1)/4.0;
    D(end, end) =  N*(N+1)/4.0;

    %% ── 5. SBP verification ──────────────────────────────────────────────
    M = diag(w);
    B = diag([-1; zeros(N-1,1); 1]);

    sbp_err = norm(M*D + (M*D)' - B, 'fro');
    con_err = max(abs(D * ones(Np1,1)));
    lin_err = max(abs(D * xi - ones(Np1,1)));

    assert(sbp_err < 1e-10, 'SBP violated: %e  (N=%d)', sbp_err, N);
    assert(con_err < 1e-12, 'Constant preservation: %e', con_err);
    assert(lin_err < 1e-12, 'Linear exactness: %e', lin_err);

    fprintf('  [D N=%d] SBP=%.2e  const=%.2e  linear=%.2e\n', ...
            N, sbp_err, con_err, lin_err);
end

% ── Local copy of leg_derivs (mirrors LGL_quadrature helper) ─────────────
function [Pnm1, Pn, dPn, d2Pn] = leg_derivs(x, N)
    if N==0, Pnm1=0; Pn=1;  dPn=0;  d2Pn=0; return; end
    if N==1, Pnm1=1; Pn=x;  dPn=1;  d2Pn=0; return; end
    Pkm2=1; Pkm1=x; dPkm2=0; dPkm1=1; d2Pkm2=0; d2Pkm1=0;
    for k=2:N
        Pk   = ((2*k-1)*x*Pkm1   - (k-1)*Pkm2)                   /k;
        dPk  = ((2*k-1)*(Pkm1  + x*dPkm1)  - (k-1)*dPkm2)        /k;
        d2Pk = ((2*k-1)*(2*dPkm1 + x*d2Pkm1) - (k-1)*d2Pkm2)     /k;
        Pkm2=Pkm1; Pkm1=Pk; dPkm2=dPkm1; dPkm1=dPk; d2Pkm2=d2Pkm1; d2Pkm1=d2Pk;
    end
    Pnm1=Pkm2; Pn=Pkm1; dPn=dPkm1; d2Pn=d2Pkm1;
    if abs(d2Pn)<1e-30, d2Pn=1e-30; end
end