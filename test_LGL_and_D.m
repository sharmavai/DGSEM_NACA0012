% =========================================================================
%  test_LGL_and_D.m
%  Validates LGL_quadrature.m and derivative_matrix.m
%  Run this BEFORE any simulation. All assertions must pass.
% =========================================================================
clear; clc;

fprintf('=== FIX #2 VALIDATION TEST ===\n\n');

for N = [2 3 4 5]
    fprintf('--- N = %d ---\n', N);
    [xi, w, V, Vinv] = LGL_quadrature(N);
    [D, M, B]        = derivative_matrix(N);

    % 1. Endpoints exact
    assert(xi(1) == -1.0 && xi(end) == 1.0, 'Endpoint error');

    % 2. Weight sum
    assert(abs(sum(w) - 2.0) < 1e-13, 'Weight sum error');

    % 3. SBP
    sbp = norm(M*D + (M*D)' - B, 'fro');
    fprintf('  SBP error     : %e  (must be < 1e-10)\n', sbp);
    assert(sbp < 1e-10);

    % 4. Constant preservation
    assert(max(abs(D*ones(N+1,1))) < 1e-12, 'Constant preservation');

    % 5. Polynomial exactness: D applied to x^k should give k*x^(k-1)
    for k = 1:N
        poly_err = max(abs(D*(xi.^k) - k*xi.^(k-1)));
        assert(poly_err < 1e-10, sprintf('Poly degree %d exactness failed: %e', k, poly_err));
    end
    fprintf('  Poly exactness: PASS (degrees 1 to %d)\n', N);

    % 6. Integration exactness: integrate x^k over [-1,1]
    for k = 0:2:2*N-1   % LGL exact for polynomials up to degree 2N-1
        exact = 2/(k+1);
        numerical = sum(w .* xi.^k);
        quad_err = abs(numerical - exact);
        assert(quad_err < 1e-12, sprintf('Quadrature degree %d failed: %e', k, quad_err));
    end
    fprintf('  Quadrature    : PASS (exact up to degree %d)\n', 2*N-1);
    fprintf('\n');
end

fprintf('=== ALL TESTS PASSED ===\n');
fprintf('LGL_quadrature.m and derivative_matrix.m are verified correct.\n');
fprintf('Safe to proceed to Fix #3 (mesh generation + metrics).\n');
