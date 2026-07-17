% =========================================================================
%  test_isentropic_vortex.m
%  Isentropic vortex verification test for DGSEM on Cartesian grid.
%
%  Uses exact solution of Euler equations to verify spatial accuracy
%  of the inviscid DGSEM operator (rhs_DGSEM.m).
%  Expected L2 errors for N=3: O(10^-4) to O(10^-6) with 4th-order
%  spatial convergence under h-refinement.
%
%  BR2 viscous terms are bypassed by setting mu=0.
%
%  References:
%    [1] Hesthaven & Warburton (2008) — isentropic vortex test case
% =========================================================================
clear; close all; clc;

fprintf('=== ISENTROPIC VORTEX VERIFICATION TEST ===\n\n');

%% ── Parameters ────────────────────────────────────────────────────────────
N       = 3;
gamma   = 1.4;
beta    = 5.0;          % vortex strength
x_c     = 0.0; y_c = 0.0;  % initial vortex center
U_conv  = 1.0;          % convection velocity
t_final = 2.0;
CFL     = 0.08;

%% ── Grid refinement levels ────────────────────────────────────────────────
nelem_list = [6 8 10 14];
domain    = [-5, 5, -5, 5];
n_levels  = numel(nelem_list);

var_names = {'rho', 'rho*u', 'rho*v', 'E'};

L2_errors = zeros(n_levels, 4);
h_vals    = zeros(n_levels, 1);

for lvl = 1:n_levels
    nelem_x = nelem_list(lvl);
    nelem_y = nelem_list(lvl);
    dx = (domain(2)-domain(1))/nelem_x;
    dy = (domain(4)-domain(3))/nelem_y;
    h_vals(lvl) = max(dx, dy);

    fprintf('--- Grid %dx%d (h=%.3f) ---\n', nelem_x, nelem_y, h_vals(lvl));

    %% ── Build mesh struct (Cartesian, curvilinear-compatible) ───────────
    [xi_lgl, w_lgl] = LGL_quadrature(N);
    [Dmat, ~, ~]    = derivative_matrix(N);
    Np1 = N + 1;
    n_elem = nelem_x * nelem_y;

    x_elem = zeros(Np1, Np1, n_elem);
    y_elem = zeros(Np1, Np1, n_elem);
    J      = zeros(Np1, Np1, n_elem);
    xi_x   = zeros(Np1, Np1, n_elem);
    xi_y   = zeros(Np1, Np1, n_elem);
    eta_x  = zeros(Np1, Np1, n_elem);
    eta_y  = zeros(Np1, Np1, n_elem);

    elem = 0;
    for jj = 1:nelem_y
        for ii = 1:nelem_x
            elem = elem + 1;
            x0 = domain(1) + (ii-1)*dx;
            y0 = domain(3) + (jj-1)*dy;
            xe_ref = zeros(Np1, Np1);
            ye_ref = zeros(Np1, Np1);
            for a = 1:Np1
                for b = 1:Np1
                    xe_ref(a,b) = x0 + 0.5*(xi_lgl(a)+1)*dx;
                    ye_ref(a,b) = y0 + 0.5*(xi_lgl(b)+1)*dy;
                end
            end
            x_elem(:,:,elem) = xe_ref;
            y_elem(:,:,elem) = ye_ref;

            x_xi_e  = Dmat * xe_ref;
            x_eta_e = xe_ref * Dmat';
            y_xi_e  = Dmat * ye_ref;
            y_eta_e = ye_ref * Dmat';
            Je = x_xi_e .* y_eta_e - x_eta_e .* y_xi_e;
            J(:,:,elem)  = Je;
            xi_x(:,:,elem)  =  y_eta_e ./ Je;
            xi_y(:,:,elem)  = -x_eta_e ./ Je;
            eta_x(:,:,elem) = -y_xi_e  ./ Je;
            eta_y(:,:,elem) =  x_xi_e  ./ Je;
        end
    end

    h_min = zeros(n_elem, 1);
    for k = 1:n_elem
        h_min(k) = min(J(:,:,k),[],'all') * 2;
    end

    % Build minimal mesh struct
    mesh.n_elem  = n_elem;
    mesh.nxi_c   = nelem_x;
    mesh.neta    = nelem_y;
    mesh.nxi_airfoil = nelem_x;
    mesh.nwake   = 0;
    mesh.N       = N;
    mesh.x       = x_elem;
    mesh.y       = y_elem;
    mesh.J       = J;
    mesh.xi_x    = xi_x;
    mesh.xi_y    = xi_y;
    mesh.eta_x   = eta_x;
    mesh.eta_y   = eta_y;
    mesh.D       = Dmat;
    mesh.w_lgl   = w_lgl;
    mesh.h_min   = h_min;
    mesh.Mach    = 0.0;
    mesh.alpha_deg = 0.0;
    mesh.wall_elems = (1:nelem_x)';
    mesh.ff_elems   = (nelem_x*(nelem_y-1)+1:nelem_x*nelem_y)';
    mesh.airfoil_ei = 1:nelem_x;
    mesh.wall_nx = zeros(Np1, nelem_x);
    mesh.wall_ny = -ones(Np1, nelem_x);
    mesh.wall_ds = ones(Np1, nelem_x) * dy * w_lgl;
    mesh.wall_x  = x_elem(:,1,1:nelem_x);
    mesh.wall_y  = y_elem(:,1,1:nelem_x);

    %% ── Initialize with exact vortex solution ────────────────────────────
    U = zeros(n_elem, Np1, Np1, 4);
    elem = 0;
    for jj = 1:nelem_y
        for ii = 1:nelem_x
            elem = elem + 1;
            for a = 1:Np1
                for b = 1:Np1
                    x = x_elem(a,b,elem);
                    y = y_elem(a,b,elem);
                    [rho, u, v, p] = vortex_exact(x, y, 0, x_c, y_c, U_conv, beta, gamma);
                    E = p/((gamma-1)*rho) + 0.5*(u^2+v^2);
                    U(elem,a,b,:) = [rho; rho*u; rho*v; rho*E];
                end
            end
        end
    end

    %% ── Time march with RK4 + rhs_DGSEM ─────────────────────────────────
    mu_v = 0;  k_v = 0;  % Inviscid test
    dt   = CFL * h_vals(lvl) / ((U_conv + 1.0) * (2*N+1));
    time = 0;
    iter = 0;

    while time < t_final
        if time + dt > t_final
            dt = t_final - time;
        end

        k1 = rhs_DGSEM(U,              mesh, gamma, mu_v, k_v);
        k2 = rhs_DGSEM(U + 0.5*dt*k1,  mesh, gamma, mu_v, k_v);
        k3 = rhs_DGSEM(U + 0.5*dt*k2,  mesh, gamma, mu_v, k_v);
        k4 = rhs_DGSEM(U + dt*k3,      mesh, gamma, mu_v, k_v);
        U  = U + (dt/6.0) * (k1 + 2*k2 + 2*k3 + k4);

        time = time + dt;
        iter = iter + 1;
    end

    %% ── Compute L2 errors against exact solution ───────────────────────
    U_ex = zeros(size(U));
    elem = 0;
    for jj = 1:nelem_y
        for ii = 1:nelem_x
            elem = elem + 1;
            for a = 1:Np1
                for b = 1:Np1
                    x = x_elem(a,b,elem);
                    y = y_elem(a,b,elem);
                    [rho, u, v, p] = vortex_exact(x, y, t_final, x_c, y_c, U_conv, beta, gamma);
                    E = p/((gamma-1)*rho) + 0.5*(u^2+v^2);
                    U_ex(elem,a,b,:) = [rho; rho*u; rho*v; rho*E];
                end
            end
        end
    end

    fprintf('  Steps: %d, L2 errors:\n', iter);
    for k = 1:4
        % L2 norm with Jacobian and quadrature weighting
        num = 0; den = 0;
        for e = 1:n_elem
            diff_k = squeeze(U(e,:,:,k) - U_ex(e,:,:,k));
            ex_k   = squeeze(U_ex(e,:,:,k));
            Je = J(:,:,e);
            num = num + sum((w_lgl(:) * w_lgl(:)' .* Je) .* diff_k.^2, 'all');
            den = den + sum((w_lgl(:) * w_lgl(:)' .* Je) .* ex_k.^2, 'all');
        end
        L2_errors(lvl, k) = sqrt(num) / sqrt(den);
        fprintf('    %-8s: %.6e\n', var_names{k}, L2_errors(lvl, k));
    end
end

%% ── Convergence rate estimation ──────────────────────────────────────────
fprintf('\n=== CONVERGENCE RATES ===\n');
fprintf('%s\n', repmat('-', 1, 60));
header = sprintf('%6s %12s %12s %12s %12s', 'h', 'rho', 'rho*u', 'rho*v', 'E');
fprintf('%s\n', header);

for lvl = 1:n_levels
    row = sprintf('%6.3f %12.4e %12.4e %12.4e %12.4e', h_vals(lvl), ...
                 L2_errors(lvl,1), L2_errors(lvl,2), L2_errors(lvl,3), L2_errors(lvl,4));
    fprintf('%s\n', row);
end

if n_levels >= 3
    fprintf('\nEstimated orders of convergence:\n');
    for k = 1:4
        p_avg = 0;
        count = 0;
        for lvl = 2:n_levels-1
            r_h = h_vals(lvl) / h_vals(lvl+1);
            if L2_errors(lvl,k) > 1e-14 && L2_errors(lvl+1,k) > 1e-14
                p = log(L2_errors(lvl,k) / L2_errors(lvl+1,k)) / log(r_h);
                p_avg = p_avg + p;
                count = count + 1;
            end
        end
        if count > 0
            fprintf('  %-8s: %.2f (expected ~4 for N=3)\n', var_names{k}, p_avg/count);
        else
            fprintf('  %-8s: (errors too small for rate estimation)\n', var_names{k});
        end
    end
end

fprintf('\n=== VORTEX VERIFICATION COMPLETE ===\n');

%% ── Helper functions ──────────────────────────────────────────────────────
function [rho, u, v, p] = vortex_exact(x, y, t, xc, yc, Uc, beta, gamma)
    xct = xc + Uc * t;
    r2  = (x - xct)^2 + (y - yc)^2;
    f   = exp((1 - r2) / 2);
    dT  = -(gamma-1) * beta^2 / (8 * gamma * pi^2) * exp(1 - r2);
    T   = 1 + dT;
    rho = T^(1/(gamma-1));
    p   = rho^gamma;
    u   = Uc - (beta / (2*pi)) * (y - yc) * f;
    v   =        (beta / (2*pi)) * (x - xct) * f;
end
