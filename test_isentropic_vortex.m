%% Isentropic Vortex — DGSEM Verification
clear; close all; clc;

N       = 3;
nelem_x = 10;
nelem_y = 10;
domain  = [-5, 5, -5, 5];
gamma   = 1.4;
beta    = 5.0;
x_c     = 0.0; y_c = 0.0;
U_conv  = 1.0;
CFL     = 0.1;
t_final = 2.0;

[xi, w] = LGL_quadrature(N);
D       = derivative_matrix(N);

dx = (domain(2)-domain(1))/nelem_x;
dy = (domain(4)-domain(3))/nelem_y;

mesh.nelem = nelem_x * nelem_y;
mesh.x     = zeros(nelem_y, nelem_x, 4);
mesh.y     = zeros(nelem_y, nelem_x, 4);
mesh.airfoil_ids  = [];
mesh.farfield_ids = [];

for j = 1:nelem_y
    for i = 1:nelem_x
        x0 = domain(1)+(i-1)*dx; y0 = domain(3)+(j-1)*dy;
        mesh.x(j,i,:) = [x0, x0+dx, x0+dx, x0];
        mesh.y(j,i,:) = [y0, y0,    y0+dy,  y0+dy];
    end
end

U = zeros(mesh.nelem, N+1, N+1, 4);
elem = 0;
for j = 1:nelem_y
    for i = 1:nelem_x
        elem = elem+1;
        for jj = 1:N+1
            for ii = 1:N+1
                x = bilinear_map(xi(ii), xi(jj), mesh.x(j,i,:));
                y = bilinear_map(xi(ii), xi(jj), mesh.y(j,i,:));
                [rho,u,v,p] = vortex_exact(x, y, 0, x_c, y_c, U_conv, beta, gamma);
                E = p/((gamma-1)*rho) + 0.5*(u^2+v^2);
                U(elem,ii,jj,:) = [rho; rho*u; rho*v; rho*E];
            end
        end
    end
end

dt   = CFL*dx / ((U_conv+1.0)*(2*N+1));
time = 0; iter = 0;
fprintf('Integrating to t=%.2f  dt=%.6f\n', t_final, dt);

while time < t_final
    iter = iter+1;
    if time+dt > t_final, dt = t_final-time; end
    % RK4 — uncomment rhs_DGSEM calls for full run
    % k1 = rhs_DGSEM(U,...); k2=...; k3=...; k4=...;
    % U = U + (dt/6)*(k1+2*k2+2*k3+k4);
    time = time+dt;
end

U_ex = zeros(size(U));
elem = 0;
for j = 1:nelem_y
    for i = 1:nelem_x
        elem = elem+1;
        for jj = 1:N+1
            for ii = 1:N+1
                x = bilinear_map(xi(ii), xi(jj), mesh.x(j,i,:));
                y = bilinear_map(xi(ii), xi(jj), mesh.y(j,i,:));
                [rho,u,v,p] = vortex_exact(x, y, t_final, x_c, y_c, U_conv, beta, gamma);
                E = p/((gamma-1)*rho) + 0.5*(u^2+v^2);
                U_ex(elem,ii,jj,:) = [rho; rho*u; rho*v; rho*E];
            end
        end
    end
end

names = {'density','x-momentum','y-momentum','energy'};
for k = 1:4
    err = sqrt(sum((U(:,:,:,k)-U_ex(:,:,:,k)).^2,'all')) / ...
          sqrt(sum(U_ex(:,:,:,k).^2,'all'));
    fprintf('L2 error %-12s: %.6e\n', names{k}, err);
end

% ── helpers ────────────────────────────────────────────────────────────────
function v = bilinear_map(xi, eta, coords)
v = 0.25*((1-xi)*(1-eta)*coords(1) + (1+xi)*(1-eta)*coords(2) + ...
          (1+xi)*(1+eta)*coords(3) + (1-xi)*(1+eta)*coords(4));
end

function [rho, u, v, p] = vortex_exact(x, y, t, xc, yc, Uc, beta, gamma)
xct = xc + Uc*t;
r2  = (x-xct)^2 + (y-yc)^2;
f   = exp((1-r2)/2);
dT  = -(gamma-1)*beta^2/(8*gamma*pi^2)*exp(1-r2);
T   = 1 + dT;
rho = T^(1/(gamma-1));
p   = rho^gamma;
u   = Uc - (beta/(2*pi))*(y-yc)*f;
v   =      (beta/(2*pi))*(x-xct)*f;
end
