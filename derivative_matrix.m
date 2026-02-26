function D = derivative_matrix(xi, N)
% DERIVATIVE_MATRIX  Spectral derivative matrix at LGL points.
%   D(i,j) = d(l_j)/dxi evaluated at xi_i.

D = zeros(N+1, N+1);

LNp = zeros(N+1, 1);
for i = 1:N+1
    [~, LNp(i)] = legendre_poly(N, xi(i));
end

for i = 1:N+1
    for j = 1:N+1
        if i ~= j
            D(i,j) = LNp(i) / (LNp(j) * (xi(i) - xi(j)));
        end
    end
end

for i = 1:N+1
    D(i,i) = -sum(D(i,:));
end

D(1,   1  ) = -N*(N+1)/4;
D(N+1, N+1) =  N*(N+1)/4;
end

function [LN, LNp] = legendre_poly(N, x)
if N == 0
    LN = 1; LNp = 0;
elseif N == 1
    LN = x; LNp = 1;
else
    L0 = 1; L1 = x; L0p = 0; L1p = 1;
    for k = 2:N
        LN  = ((2*k-1)*x*L1 - (k-1)*L0) / k;
        LNp = L0p + (2*k-1)*L1;
        L0  = L1;  L1  = LN;
        L0p = L1p; L1p = LNp;
    end
end
end
