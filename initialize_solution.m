function U = initialize_solution(mesh, N, rho_inf, U_inf, alpha, p_inf, gamma)
% INITIALIZE_SOLUTION  Uniform free-stream initial condition.
%   U: [nelem, N+1, N+1, 4]  conservative variables (rho, rhou, rhov, rhoE)

alpha_rad = alpha * pi / 180;
u_inf = U_inf * cos(alpha_rad);
v_inf = U_inf * sin(alpha_rad);
E_inf = p_inf/((gamma-1)*rho_inf) + 0.5*(u_inf^2 + v_inf^2);

U = zeros(mesh.nelem, N+1, N+1, 4);
U(:,:,:,1) = rho_inf;
U(:,:,:,2) = rho_inf * u_inf;
U(:,:,:,3) = rho_inf * v_inf;
U(:,:,:,4) = rho_inf * E_inf;
end
