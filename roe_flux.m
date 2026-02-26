function F_star = roe_flux(UL, UR, nx, ny, gamma)
% ROE_FLUX  Roe approximate Riemann solver with entropy fix.

rhoL = UL(1); uL = UL(2)/rhoL; vL = UL(3)/rhoL; EL = UL(4)/rhoL;
pL = (gamma-1)*rhoL*(EL - 0.5*(uL^2+vL^2)); HL = EL + pL/rhoL;

rhoR = UR(1); uR = UR(2)/rhoR; vR = UR(3)/rhoR; ER = UR(4)/rhoR;
pR = (gamma-1)*rhoR*(ER - 0.5*(uR^2+vR^2)); HR = ER + pR/rhoR;

sqL = sqrt(rhoL); sqR = sqrt(rhoR);
u_r = (sqL*uL + sqR*uR)/(sqL+sqR);
v_r = (sqL*vL + sqR*vR)/(sqL+sqR);
H_r = (sqL*HL + sqR*HR)/(sqL+sqR);
c_r = sqrt((gamma-1)*(H_r - 0.5*(u_r^2+v_r^2)));

Vn = u_r*nx + v_r*ny;
rho_r = sqL*sqR;
eps = 0.1*c_r;

l1 = entropy_fix(Vn - c_r, eps);
l2 = abs(Vn);
l4 = entropy_fix(Vn + c_r, eps);

FL = normal_flux(UL, nx, ny, gamma);
FR = normal_flux(UR, nx, ny, gamma);

dU = UR - UL;
drho = dU(1); drhou = dU(2); drhov = dU(3); drhoE = dU(4);

dp   = (gamma-1)*((H_r - u_r^2 - v_r^2)*drho + u_r*drhou + v_r*drhov - drhoE);
dVn  = (drhou*nx + drhov*ny)/rho_r - Vn*drho/rho_r;

a1 = (dp - rho_r*c_r*dVn)/(2*c_r^2);
a2 = drho - dp/c_r^2;
a5 = (dp + rho_r*c_r*dVn)/(2*c_r^2);

K1 = l1*a1*[1; u_r-c_r*nx; v_r-c_r*ny; H_r-c_r*Vn];
K2 = l2*a2*[1; u_r; v_r; 0.5*(u_r^2+v_r^2)];
K5 = l4*a5*[1; u_r+c_r*nx; v_r+c_r*ny; H_r+c_r*Vn];

F_star = 0.5*(FL+FR) - 0.5*(K1+K2+K5);
end

function la = entropy_fix(lambda, eps)
if abs(lambda) >= eps
    la = abs(lambda);
else
    la = (lambda^2 + eps^2)/(2*eps);
end
end

function Fn = normal_flux(U, nx, ny, gamma)
rho = U(1); u = U(2)/rho; v = U(3)/rho; E = U(4)/rho;
p = (gamma-1)*rho*(E - 0.5*(u^2+v^2)); H = E + p/rho;
Vn = u*nx + v*ny;
Fn = [rho*Vn; rho*u*Vn+p*nx; rho*v*Vn+p*ny; rho*H*Vn];
end
