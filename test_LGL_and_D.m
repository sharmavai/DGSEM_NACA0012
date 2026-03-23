% =========================================================================
%  test_LGL_and_D.m  — validates LGL_quadrature.m and derivative_matrix.m
% =========================================================================
clear; clc;
fprintf('=== FIX #2 VALIDATION TEST ===\n\n');

for N = [2 3 4 5]
    fprintf('--- N = %d ---\n', N);
    [xi, w]      = LGL_quadrature(N);
    [D, M, B]    = derivative_matrix(N);

    % 1. Endpoints exact
    assert(xi(1)==-1.0 && xi(end)==1.0, 'Endpoint error');

    % 2. Weight sum
    assert(abs(sum(w)-2.0) < 1e-12, 'Weight sum error: %.15f', sum(w));

    % 3. SBP
    sbp = norm(M*D + (M*D)' - B, 'fro');
    assert(sbp < 1e-10, 'SBP error: %e', sbp);

    % 4. Constant preservation
    assert(max(abs(D*ones(N+1,1))) < 1e-12, 'Constant preservation failed');

    % 5. Polynomial exactness degrees 1..N
    for k = 1:N
        err = max(abs(D*(xi.^k) - k*xi.^(k-1)));
        assert(err < 1e-10, 'Poly degree %d: %e', k, err);
    end
    fprintf('  Poly exactness: PASS (1..%d)\n', N);

    % 6. Quadrature exactness (exact up to degree 2N-1)
    for k = 0:2:2*N-2
        err = abs(sum(w.*xi.^k) - 2/(k+1));
        assert(err < 1e-12, 'Quadrature degree %d: %e', k, err);
    end
    fprintf('  Quadrature:     PASS (up to degree %d)\n\n', 2*N-2);
end

fprintf('=== ALL TESTS PASSED ===\n');