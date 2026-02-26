function [xi, w] = LGL_quadrature(N)
% LGL_QUADRATURE  Legendre-Gauss-Lobatto points and weights on [-1,1]
%   [xi, w] = LGL_QUADRATURE(N)  returns N+1 points and weights.

xi = zeros(N+1, 1);
w  = zeros(N+1, 1);

xi(1)   = -1.0;
xi(N+1) =  1.0;

if N == 3
    xi(2) = -sqrt(1/5);
    xi(3) =  sqrt(1/5);
    w = [1/6; 5/6; 5/6; 1/6];
    return;
end

if N > 1
    for i = 2:N
        xi(i) = -cos(pi * (i-1) / N);
    end
    for i = 2:N
        for iter = 1:100
            [LN, LNp] = legendre_poly(N, xi(i));
            LNpp = legendre_poly_second_deriv(N, xi(i));
            xi(i) = xi(i) - ((1 - xi(i)^2) * LNp) / ...
                    (-2*xi(i)*LNp + (1 - xi(i)^2)*LNpp);
        end
    end
end

for i = 1:N+1
    [LN, ~] = legendre_poly(N, xi(i));
    w(i) = 2 / (N*(N+1)*LN^2);
end
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

function LNpp = legendre_poly_second_deriv(N, x)
[LN, LNp] = legendre_poly(N, x);
LNpp = (2*x*LNp - N*(N+1)*LN) / (1 - x^2);
end
