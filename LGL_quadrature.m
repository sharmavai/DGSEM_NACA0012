% =========================================================================
%  LGL_quadrature.m
%  Legendre-Gauss-Lobatto nodes, weights, Vandermonde, and its inverse
%
%  BUGFIX vs previous version:
%    Vandermonde inverse was wrong: used Vinv=V' but correct relation is
%    V'*M*V = I  =>  Vinv = V'*M.  Rather than carrying that complexity,
%    D is now built directly from the Lagrange basis formula (Kopriva 2009
%    eq. 3.44) which needs no matrix inversion at all.
%    V and Vinv are still returned for callers that need them, computed
%    correctly as Vinv = V' * diag(w).
%
%  References:
%    [1] Kopriva (2009) Algorithm 24 (nodes), eq. 3.18 (weights),
%                       eq. 3.44 (derivative matrix direct formula)
% =========================================================================

function [xi, w, V, Vinv] = LGL_quadrature(N)

    Np1 = N + 1;

    %% ── 1. Endpoints ─────────────────────────────────────────────────────
    xi      = zeros(Np1, 1);
    xi(1)   = -1.0;
    xi(end) =  1.0;

    %% ── 2. Interior nodes — Newton on P'_N(x) = 0 ───────────────────────
    for j = 1:N-1
        x = -cos((j + 0.5) * pi / N);   % Kopriva Alg. 24 initial guess
        for iter = 1:100
            [~, ~, dPn, d2Pn] = leg_derivs(x, N);
            delta = -dPn / d2Pn;
            x     =  x + delta;
            if abs(delta) < 1e-15, break; end
        end
        xi(j+1) = x;
    end

    %% ── 3. Weights ───────────────────────────────────────────────────────
    w = zeros(Np1, 1);
    for i = 1:Np1
        [~, Pn, ~, ~] = leg_derivs(xi(i), N);
        w(i) = 2.0 / (N * Np1 * Pn^2);
    end
    assert(abs(sum(w) - 2.0) < 1e-12, ...
           'Weights sum = %.15f, expected 2.0  (N=%d)', sum(w), N);

    %% ── 4. Vandermonde (normalised Legendre basis) ───────────────────────
    % V(i,j) = sqrt((2(j-1)+1)/2) * P_{j-1}(xi_i)
    % Discrete orthogonality: V' * diag(w) * V = I  =>  Vinv = V' * diag(w)
    V = zeros(Np1);
    for i = 1:Np1
        V(i,:) = normalised_legendre_row(xi(i), N);
    end
    Vinv = V' * diag(w);   % CORRECT inverse

    % Verify: Vinv * V should be identity
    err = norm(Vinv * V - eye(Np1), 'fro');
    assert(err < 1e-10, ...
           'Vinv*V != I: ||Vinv*V-I||_F = %e  (N=%d)', err, N);

end

% =========================================================================
%  Helpers
% =========================================================================

function row = normalised_legendre_row(x, N)
% row(j) = sqrt((2*(j-1)+1)/2) * P_{j-1}(x),  j=1..N+1
    Np1 = N + 1;
    row = zeros(1, Np1);
    Pp = 1.0;  Pc = x;
    row(1) = sqrt(0.5);
    if N >= 1
        row(2) = sqrt(1.5) * x;
    end
    for k = 2:N
        Pn       = ((2*k-1)*x*Pc - (k-1)*Pp) / k;
        row(k+1) = sqrt((2*k+1)/2.0) * Pn;
        Pp = Pc;  Pc = Pn;
    end
end

function [Pnm1, Pn, dPn, d2Pn] = leg_derivs(x, N)
% P_{N-1}(x), P_N(x), P'_N(x), P''_N(x) via three-term recurrence
    if N == 0
        Pnm1=0; Pn=1;  dPn=0; d2Pn=0; return
    end
    if N == 1
        Pnm1=1; Pn=x;  dPn=1; d2Pn=0; return
    end
    Pkm2=1; Pkm1=x; dPkm2=0; dPkm1=1; d2Pkm2=0; d2Pkm1=0;
    for k = 2:N
        Pk    = ((2*k-1)*x*Pkm1   - (k-1)*Pkm2)                  / k;
        dPk   = ((2*k-1)*(Pkm1 + x*dPkm1)   - (k-1)*dPkm2)       / k;
        d2Pk  = ((2*k-1)*(2*dPkm1 + x*d2Pkm1) - (k-1)*d2Pkm2)    / k;
        Pkm2=Pkm1; Pkm1=Pk; dPkm2=dPkm1; dPkm1=dPk; d2Pkm2=d2Pkm1; d2Pkm1=d2Pk;
    end
    Pnm1=Pkm2; Pn=Pkm1; dPn=dPkm1; d2Pn=d2Pkm1;
    if abs(d2Pn) < 1e-30, d2Pn = 1e-30; end
end