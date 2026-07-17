% =========================================================================
%  rhs_DGSEM.m
%  Spatial residual for 2D compressible NAVIER-STOKES equations via DGSEM
%  (strong-form SBP-SAT on curvilinear C-grid with BR2 viscous scheme)
%
%  INVISCID:
%    Volume integral: D_xi(Fh) + D_eta(Gh)  where Fh, Gh are contravariant
%    fluxes J*(xi_x*F + xi_y*G), J*(eta_x*F + eta_y*G)
%    Surface: Lax-Friedrichs (Rusanov) numerical flux at all faces
%
%  VISCOUS (BR2 scheme, Bassi & Rebay 1997, extended by Cockburn & Shu):
%    Volume: D_xi(Fv_xi) + D_eta(Fv_eta)  with curvilinear metric terms
%    Surface: BR2 lifting operator r_eta^* with penalty parameter eta_B = 1
%    Viscous contravariant fluxes use velocity gradients computed at LGL nodes
%
%  SAT signs (SBP framework):
%    Right face  (xi=+1, i=Np1):  +sJ*(Fnum - F_int) / w(Np1)
%    Left face   (xi=-1, i=1):   -sJ*(Fnum - F_int) / w(1)
%    Right face  (eta=+1, j=Np1): +sJ*(Fnum - F_int) / w(Np1)
%    Left face   (eta=-1, j=1):  -sJ*(Fnum - F_int) / w(1)
%    (F_int is the flux evaluated with the INTERIOR state at the face)
%
%  Q layout: (n_elem, Np1, Np1, 4)
%    Q(e,i,j,1:4) = [rho, rho*u, rho*v, E]
%    i = xi-node (circumferential), j = eta-node (radial)
%
%  References:
%    [1] Kopriva (2009) Ch. 8 — DGSEM with metric terms
%    [2] Hesthaven & Warburton (2008) — SBP-SAT framework
%    [3] Bassi & Rebay (1997) J. Comput. Phys. — BR2 scheme
%    [4] Cockburn & Shu (1998) Math. Comp. — BR2 lifting operators
%    [5] Ferrer et al. (2021) DGSEM iLES NACA0012, NNFM 143
% =========================================================================

function dUdt = rhs_DGSEM(U, mesh, gamma, mu, k_cond)

n_elem = mesh.n_elem;
nxi_c  = mesh.nxi_c;
neta   = mesh.neta;
N      = mesh.N;
Np1    = N + 1;
D      = mesh.D;
w      = mesh.w_lgl;

Mach      = mesh.Mach;
alpha_deg = mesh.alpha_deg;

% BR2 penalty parameter
eta_BR2 = 1.0;

dUdt = zeros(size(U));

% ── Precompute viscous gradient lifting terms ─────────────────────────────
% For BR2 we need to compute velocity gradients, then lift the jump of
% gradient traces at element interfaces. We store: r_xi, r_eta (lifting
% of gradient jumps) per element.
%
% STEP 1: Compute velocity gradients in each element
du_dxi  = zeros(Np1, Np1, n_elem);   dv_dxi  = zeros(Np1, Np1, n_elem);
du_deta = zeros(Np1, Np1, n_elem);   dv_deta = zeros(Np1, Np1, n_elem);

for e = 1:n_elem
    Qe = squeeze(U(e,:,:,:));
    rho = max(Qe(:,:,1), 1e-12);
    uu = Qe(:,:,2) ./ rho;
    vv = Qe(:,:,3) ./ rho;

    % Physical gradients via metric terms:
    % du/dx = xi_x*du/dxi + eta_x*du/deta
    % du_dxi and du_deta are spectral derivatives of u in reference space
    du_dxi_ref  = D * uu;          % (Np1, Np1) — du/dxi
    du_deta_ref = (uu * D')';      % (Np1, Np1) — du/deta
    dv_dxi_ref  = D * vv;
    dv_deta_ref = (vv * D')';

    % Store reference-space derivatives (we'll apply metric terms later)
    du_dxi(:,:,e)  = du_dxi_ref;
    du_deta(:,:,e) = du_deta_ref;
    dv_dxi(:,:,e)  = dv_dxi_ref;
    dv_deta(:,:,e) = dv_deta_ref;
end

% STEP 2: BR2 lifting of gradient jumps at element interfaces
% r_xi: jump of gradient traces at xi-faces, weighted by penalty
% r_eta: jump of gradient traces at eta-faces, weighted by penalty
% These are stored per element: rx(e,:,:,:), ry(e,:,:,:)
% but since BR2 uses the SAME penalty at both sides, we accumulate
% contributions from each face into the neighbor.

% Lifting operators (per element, reference-space gradient corrections)
rxi_lift = zeros(Np1, Np1, n_elem);  % correction to d/dxi component
reta_lift = zeros(Np1, Np1, n_elem); % correction to d/deta component
ryi_lift = zeros(Np1, Np1, n_elem);
ryeta_lift = zeros(Np1, Np1, n_elem);

for e = 1:n_elem
    ei = mod(e-1, nxi_c) + 1;
    ej = floor((e-1) / nxi_c) + 1;

    Qe = squeeze(U(e,:,:,:));
    xe = mesh.x(:,:,e);
    ye = mesh.y(:,:,e);

    % ── xi=+1 face (i=Np1) ──────────────────────────────────────────────
    ei_p = mod(ei, nxi_c) + 1;
    nb   = ei_p + (ej-1)*nxi_c;
    Qnb  = squeeze(U(nb,:,:,:));

    % Reference-space face tangent vector: (dx/deta, dy/deta) at i=Np1
    x_eta_f = (xe(Np1,:) * D')';
    y_eta_f = (ye(Np1,:) * D')';

    for j = 1:Np1
        sJ = sqrt(y_eta_f(j)^2 + x_eta_f(j)^2);
        if sJ < 1e-30, continue; end

        % Average reference gradient at face
        du_avg = 0.5*(du_dxi(Np1,j,e) + du_dxi(1,j,nb));
        dv_avg = 0.5*(dv_dxi(Np1,j,e) + dv_dxi(1,j,nb));

        % Jump (interior - neighbor) of gradient trace
        jump_u = du_dxi(Np1,j,e) - du_dxi(1,j,nb);
        jump_v = dv_dxi(Np1,j,e) - dv_dxi(1,j,nb);

        % BR2 lifting: add to both elements
        pen = eta_BR2 / w(Np1);
        rxi_lift(Np1,j,e)   = rxi_lift(Np1,j,e)   + pen * jump_u;
        rxi_lift(1,j,nb)     = rxi_lift(1,j,nb)     + pen * jump_u;
        ryi_lift(Np1,j,e)    = ryi_lift(Np1,j,e)    + pen * jump_v;
        ryi_lift(1,j,nb)     = ryi_lift(1,j,nb)     + pen * jump_v;
    end

    % ── eta=+1 face (j=Np1) ─────────────────────────────────────────────
    x_xi_f = D * xe(:,Np1);
    y_xi_f = D * ye(:,Np1);

    if ej < neta
        nb = ei + ej*nxi_c;
        for i = 1:Np1
            sJ = sqrt(y_xi_f(i)^2 + x_xi_f(i)^2);
            if sJ < 1e-30, continue; end

            jump_u = du_deta(i,Np1,e) - du_deta(i,1,nb);
            jump_v = dv_deta(i,Np1,e) - dv_deta(i,1,nb);

            pen = eta_BR2 / w(Np1);
            reta_lift(i,Np1,e) = reta_lift(i,Np1,e) + pen * jump_u;
            reta_lift(i,1,nb)   = reta_lift(i,1,nb)   + pen * jump_u;
            ryeta_lift(i,Np1,e) = ryeta_lift(i,Np1,e) + pen * jump_v;
            ryeta_lift(i,1,nb)   = ryeta_lift(i,1,nb)   + pen * jump_v;
        end
    else
        % Far-field boundary: zero-gradient lift (no jump correction)
        % r_lift already initialized to zero — do nothing
    end

    % ── eta=-1 face (j=1) ────────────────────────────────────────────────
    if ej > 1
        nb = ei + (ej-2)*nxi_c;
        for i = 1:Np1
            sJ = sqrt(y_xi_f(i)^2 + x_xi_f(i)^2);
            if sJ < 1e-30, continue; end

            jump_u = du_deta(i,1,e) - du_deta(i,Np1,nb);
            jump_v = dv_deta(i,1,e) - dv_deta(i,Np1,nb);

            pen = eta_BR2 / w(1);
            reta_lift(i,1,e)     = reta_lift(i,1,e)     + pen * jump_u;
            reta_lift(i,Np1,nb)  = reta_lift(i,Np1,nb)  + pen * jump_u;
            ryeta_lift(i,1,e)    = ryeta_lift(i,1,e)    + pen * jump_v;
            ryeta_lift(i,Np1,nb) = ryeta_lift(i,Np1,nb) + pen * jump_v;
        end
    else
        % Wall boundary (ej=1): lifting for no-slip condition
        for i = 1:Np1
            % At wall, the gradient jump includes the ghost state.
            % Ghost velocity = -interior velocity (mirror), so
            % d(u)/dxi|_ghost = d(-u)/dxi|_int = -du/dxi|_int
            % But the lifting only needs the jump in the trace of
            % the COMMON gradient (not the velocity itself).
            % For BR2 at a solid wall, we set the trace of the gradient
            % from the ghost side to zero (no-slip enforced weakly).
            % The jump is then: (du_deta|_int - 0) = du_deta|_int
            jump_u = du_deta(i,1,e);
            jump_v = dv_deta(i,1,e);

            pen = eta_BR2 / w(1);
            reta_lift(i,1,e)   = reta_lift(i,1,e)   + pen * jump_u;
            ryeta_lift(i,1,e)   = ryeta_lift(i,1,e)   + pen * jump_v;
        end
    end
end

% ── MAIN ELEMENT LOOP: assemble full residual ──────────────────────────
for e = 1:n_elem
    ei = mod(e-1, nxi_c) + 1;
    ej = floor((e-1) / nxi_c) + 1;

    Qe = squeeze(U(e,:,:,:));
    xe = mesh.x(:,:,e);
    ye = mesh.y(:,:,e);
    Je = mesh.J(:,:,e);

    % ── Primitive variables ──────────────────────────────────────────────
    rho = max(Qe(:,:,1), 1e-12);
    uu  = Qe(:,:,2) ./ rho;
    vv  = Qe(:,:,3) ./ rho;
    EE  = Qe(:,:,4) ./ rho;
    pp  = max((gamma-1)*rho.*(EE - 0.5*(uu.^2 + vv.^2)), 1e-12);
    HH  = EE + pp./rho;

    % ── INVISCID physical fluxes ─────────────────────────────────────────
    F = zeros(Np1, Np1, 4);  G = zeros(Np1, Np1, 4);
    F(:,:,1) = rho.*uu;            G(:,:,1) = rho.*vv;
    F(:,:,2) = rho.*uu.*uu + pp;   G(:,:,2) = rho.*uu.*vv;
    F(:,:,3) = rho.*uu.*vv;        G(:,:,3) = rho.*vv.*vv + pp;
    F(:,:,4) = rho.*uu.*HH;        G(:,:,4) = rho.*vv.*HH;

    % ── INVISCID contravariant fluxes ────────────────────────────────────
    xi_xe  = mesh.xi_x(:,:,e);    xi_ye  = mesh.xi_y(:,:,e);
    eta_xe = mesh.eta_x(:,:,e);   eta_ye = mesh.eta_y(:,:,e);

    Fh = zeros(Np1, Np1, 4);
    Gh = zeros(Np1, Np1, 4);
    for k = 1:4
        Fh(:,:,k) = Je .* (xi_xe.*F(:,:,k)  + xi_ye.*G(:,:,k));
        Gh(:,:,k) = Je .* (eta_xe.*F(:,:,k) + eta_ye.*G(:,:,k));
    end

    % ── VISCOUS gradients in physical space ──────────────────────────────
    % Corrected gradient = reference gradient + lifting correction
    % Then: du/dx = xi_x*du_dxi_corrected + eta_x*du_deta_corrected
    % etc.

    du_dxi_c  = du_dxi(:,:,e)  + rxi_lift(:,:,e);
    du_deta_c = du_deta(:,:,e) + reta_lift(:,:,e);
    dv_dxi_c  = dv_dxi(:,:,e)  + ryi_lift(:,:,e);
    dv_deta_c = dv_deta(:,:,e) + ryeta_lift(:,:,e);

    % Physical velocity gradients
    du_dx = xi_xe .* du_dxi_c + eta_xe .* du_deta_c;
    du_dy = xi_ye .* du_dxi_c + eta_ye .* du_deta_c;
    dv_dx = xi_xe .* dv_dxi_c + eta_xe .* dv_deta_c;
    dv_dy = xi_ye .* dv_dxi_c + eta_ye .* dv_deta_c;

    % Velocity divergence
    div_u = du_dx + dv_dy;

    % Temperature (for energy viscous flux)
    T_field = pp ./ rho;

    % ── VISCOUS contravariant fluxes ──────────────────────────────────────
    % Physical viscous fluxes (Newtonian fluid):
    %   Fv = [0; tau_xx; tau_xy; tau_xx*u + tau_xy*v + k*dT/dx]
    %   Gv = [0; tau_xy; tau_yy; tau_xy*u + tau_yy*v + k*dT/dy]
    % where tau_ij = mu*(du_i/dx_j + du_j/dx_i)  [for div(u) = 0 form,
    % or tau_ij = mu*(du_i/dx_j + du_j/dx_i - 2/3*div(u)*delta_ij)]
    %
    % For compressible NS, we use the Stokes hypothesis:
    %   tau_xx = mu*(2*du/dx - 2/3*div(u))
    %   tau_yy = mu*(2*dv/dy - 2/3*div(u))
    %   tau_xy = mu*(du/dy + dv/dx)

    tau_xx = mu .* (2*du_dx - (2.0/3.0)*div_u);
    tau_yy = mu .* (2*dv_dy - (2.0/3.0)*div_u);
    tau_xy = mu .* (du_dy + dv_dx);

    % Temperature gradients
    dT_dx = xi_xe .* ((D * T_field)) + eta_xe .* ((T_field * D')');
    dT_dy = xi_ye .* ((D * T_field)) + eta_ye .* ((T_field * D')');

    % Viscous energy flux
    theta_x = tau_xx .* uu + tau_xy .* vv + k_cond .* dT_dx;
    theta_y = tau_xy .* uu + tau_yy .* vv + k_cond .* dT_dy;

    % Viscous physical fluxes Fv, Gv
    Fv = zeros(Np1, Np1, 4);  Gv = zeros(Np1, Np1, 4);
    Fv(:,:,2) = tau_xx;    Gv(:,:,2) = tau_xy;
    Fv(:,:,3) = tau_xy;    Gv(:,:,3) = tau_yy;
    Fv(:,:,4) = theta_x;   Gv(:,:,4) = theta_y;

    % Viscous contravariant fluxes
    Fvh = zeros(Np1, Np1, 4);
    Gvh = zeros(Np1, Np1, 4);
    for k = 2:4
        Fvh(:,:,k) = Je .* (xi_xe .* Fv(:,:,k)  + xi_ye .* Gv(:,:,k));
        Gvh(:,:,k) = Je .* (eta_xe .* Fv(:,:,k) + eta_ye .* Gv(:,:,k));
    end

    % ── Volume integral: inviscid + viscous ────────────────────────────
    Res = zeros(Np1, Np1, 4);
    for k = 1:4
        Res_invisc = D * Fh(:,:,k) + (D * Gh(:,:,k)')';
        Res_visc   = zeros(Np1, Np1);
        if k >= 2
            Res_visc = D * Fvh(:,:,k) + (D * Gvh(:,:,k)')';
        end
        Res(:,:,k) = Res_invisc - Res_visc;
    end

    % Max wave speed in element
    cc    = sqrt(gamma * pp ./ rho);
    lam_e = max(sqrt(uu(:).^2 + vv(:).^2) + cc(:));

    % ── FACE 1: xi = +1 (i=Np1), circumferential +1 ────────────────────
    ei_p = mod(ei, nxi_c) + 1;
    nb   = ei_p + (ej-1)*nxi_c;
    Qnb  = squeeze(U(nb,:,:,:));
    y_eta_f = (ye(Np1,:) * D')';   x_eta_f = (xe(Np1,:) * D')';
    for j = 1:Np1
        sJ = sqrt(y_eta_f(j)^2 + x_eta_f(j)^2);
        if sJ < 1e-30, continue; end
        nx =  y_eta_f(j)/sJ;   ny = -x_eta_f(j)/sJ;
        QL4 = squeeze(Qe(Np1,j,:));
        QR4 = squeeze(Qnb(1,j,:));
        lam = max(lam_e, ws(QR4,gamma));
        Fn_L = nflux(QL4,nx,ny,gamma);
        Fn_R = nflux(QR4,nx,ny,gamma);
        Fn_num = 0.5*(Fn_L+Fn_R) - 0.5*lam*(QR4(:)-QL4(:));
        % SAT sign: right face (i=Np1) => +sJ/w(Np1) * (Fnum - F_int)
        for k=1:4
            Res(Np1,j,k) = Res(Np1,j,k) + sJ*(Fn_num(k)-Fn_L(k))/w(Np1);
        end
        % Viscous surface flux (BR2): use corrected gradient trace
        if mu > 0
            % Average corrected gradient at face from both sides
            % Physical gradient of u in normal direction (outward from e)
            % n = (y_eta, -x_eta)/sJ
            % d(u)/dn = (y_eta*du/dx - x_eta*du/dy)/sJ (approx)
            % For BR2, the surface flux uses the common gradient
            % evaluated with the lifting corrections
            ux_L = du_dx(Np1,j); uy_L = du_dy(Np1,j);
            vx_L = dv_dx(Np1,j); vy_L = dv_dy(Np1,j);

            % Neighbor gradient (compute on-the-fly for viscous face term)
            du_dxi_nb = du_dxi(:,:,nb) + rxi_lift(:,:,nb);
            du_deta_nb = du_deta(:,:,nb) + reta_lift(:,:,nb);
            dv_dxi_nb = dv_dxi(:,:,nb) + ryi_lift(:,:,nb);
            dv_deta_nb = dv_deta(:,:,nb) + ryeta_lift(:,:,nb);
            xi_xn = mesh.xi_x(:,:,nb); xi_yn = mesh.xi_y(:,:,nb);
            eta_xn = mesh.eta_x(:,:,nb); eta_yn = mesh.eta_y(:,:,nb);
            ux_R = xi_xn(1,j)*du_dxi_nb(1,j) + eta_xn(1,j)*du_deta_nb(1,j);
            uy_R = xi_yn(1,j)*du_dxi_nb(1,j) + eta_yn(1,j)*du_deta_nb(1,j);
            vx_R = xi_xn(1,j)*dv_dxi_nb(1,j) + eta_xn(1,j)*dv_deta_nb(1,j);
            vy_R = xi_yn(1,j)*dv_dxi_nb(1,j) + eta_yn(1,j)*dv_deta_nb(1,j);

            % Average (common) gradient at face
            ux_avg = 0.5*(ux_L + ux_R); uy_avg = 0.5*(uy_L + uy_R);
            vx_avg = 0.5*(vx_L + vx_R); vy_avg = 0.5*(vy_L + vy_R);

            % No-slip if neighbor is wall (not applicable at xi faces)

            % Normal gradient components
            dudn = ux_avg*nx + uy_avg*ny;
            dvdn = vx_avg*nx + vy_avg*ny;

            % Physical normal viscous flux at face
            tau_nn = mu*(2*dudn*nx*nx + 2*dvdn*ny*ny ...
                       - (2.0/3.0)*(ux_avg+vy_avg) ...
                       + (ux_avg*ny + uy_avg*nx)*ny ...
                       + (vx_avg*ny + vy_avg*nx)*nx);
            % Simplified: tau . n = [tau_xx*nx + tau_xy*ny, tau_xy*nx + tau_yy*ny]
            Fvnx = (mu*(2*ux_avg - (2.0/3.0)*(ux_avg+vy_avg)))*nx ...
                  + mu*(uy_avg + vx_avg)*ny;
            Fvny = mu*(ux_avg + vy_avg)*nx ...
                  + (mu*(2*vy_avg - (2.0/3.0)*(ux_avg+vy_avg)))*ny;

            % Viscous momentum flux SAT: subtract from residual
            Res(Np1,j,2) = Res(Np1,j,2) - sJ*Fvnx/w(Np1);
            Res(Np1,j,3) = Res(Np1,j,3) - sJ*Fvny/w(Np1);
        end
    end

    % ── FACE 2: xi = -1 (i=1), circumferential -1 ───────────────────────
    ei_m = mod(ei-2, nxi_c) + 1;
    nb   = ei_m + (ej-1)*nxi_c;
    Qnb  = squeeze(U(nb,:,:,:));
    y_eta_f = (ye(1,:) * D')';   x_eta_f = (xe(1,:) * D')';
    for j = 1:Np1
        sJ = sqrt(y_eta_f(j)^2 + x_eta_f(j)^2);
        if sJ < 1e-30, continue; end
        nx = y_eta_f(j)/sJ;  ny = -x_eta_f(j)/sJ;
        % SAT sign: left face (i=1) => -sJ/w(1) * (Fnum - F_right)
        % where F_right is the flux with the RIGHT state = interior at i=1
        Q_left  = squeeze(Qnb(Np1,j,:));
        Q_right = squeeze(Qe(1,j,:));
        lam = max(lam_e, ws(Q_left,gamma));
        Fn_L = nflux(Q_left,nx,ny,gamma);
        Fn_R = nflux(Q_right,nx,ny,gamma);
        Fn_num = 0.5*(Fn_L+Fn_R) - 0.5*lam*(Q_right(:)-Q_left(:));
        for k=1:4
            Res(1,j,k) = Res(1,j,k) - sJ*(Fn_num(k)-Fn_R(k))/w(1);
        end
        % Viscous surface flux
        if mu > 0
            ux_R = du_dx(1,j); uy_R = du_dy(1,j);
            vx_R = dv_dx(1,j); vy_R = dv_dy(1,j);
            du_dxi_nb = du_dxi(:,:,nb) + rxi_lift(:,:,nb);
            du_deta_nb = du_deta(:,:,nb) + reta_lift(:,:,nb);
            dv_dxi_nb = dv_dxi(:,:,nb) + ryi_lift(:,:,nb);
            dv_deta_nb = dv_deta(:,:,nb) + ryeta_lift(:,:,nb);
            xi_xn = mesh.xi_x(:,:,nb); xi_yn = mesh.xi_y(:,:,nb);
            eta_xn = mesh.eta_x(:,:,nb); eta_yn = mesh.eta_y(:,:,nb);
            ux_L = xi_xn(Np1,j)*du_dxi_nb(Np1,j) + eta_xn(Np1,j)*du_deta_nb(Np1,j);
            uy_L = xi_yn(Np1,j)*du_dxi_nb(Np1,j) + eta_yn(Np1,j)*du_deta_nb(Np1,j);
            vx_L = xi_xn(Np1,j)*dv_dxi_nb(Np1,j) + eta_xn(Np1,j)*dv_deta_nb(Np1,j);
            vy_L = xi_yn(Np1,j)*dv_dxi_nb(Np1,j) + eta_yn(Np1,j)*dv_deta_nb(Np1,j);
            ux_avg = 0.5*(ux_L+ux_R); uy_avg = 0.5*(uy_L+uy_R);
            vx_avg = 0.5*(vx_L+vx_R); vy_avg = 0.5*(vy_L+vy_R);
            Fvnx = (mu*(2*ux_avg-(2.0/3.0)*(ux_avg+vy_avg)))*nx ...
                  + mu*(uy_avg+vx_avg)*ny;
            Fvny = mu*(ux_avg+vy_avg)*nx ...
                  + (mu*(2*vy_avg-(2.0/3.0)*(ux_avg+vy_avg)))*ny;
            Res(1,j,2) = Res(1,j,2) + sJ*Fvnx/w(1);
            Res(1,j,3) = Res(1,j,3) + sJ*Fvny/w(1);
        end
    end

    % ── FACE 3: eta = +1 (j=Np1), radial outward ────────────────────────
    y_xi_f = D * ye(:,Np1);   x_xi_f = D * xe(:,Np1);
    if ej < neta
        nb  = ei + ej*nxi_c;
        Qnb = squeeze(U(nb,:,:,:));
        for i = 1:Np1
            sJ = sqrt(y_xi_f(i)^2 + x_xi_f(i)^2);
            if sJ < 1e-30, continue; end
            nx = -y_xi_f(i)/sJ;  ny = x_xi_f(i)/sJ;
            QL4 = squeeze(Qe(i,Np1,:));
            QR4 = squeeze(Qnb(i,1,:));
            lam = max(lam_e, ws(QR4,gamma));
            Fn_L = nflux(QL4,nx,ny,gamma);
            Fn_R = nflux(QR4,nx,ny,gamma);
            Fn_num = 0.5*(Fn_L+Fn_R) - 0.5*lam*(QR4(:)-QL4(:));
            % Right face (j=Np1) => +sJ/w(Np1) * (Fnum - F_int)
            for k=1:4
                Res(i,Np1,k) = Res(i,Np1,k) + sJ*(Fn_num(k)-Fn_L(k))/w(Np1);
            end
            if mu > 0
                ux_L = du_dx(i,Np1); uy_L = du_dy(i,Np1);
                vx_L = dv_dx(i,Np1); vy_L = dv_dy(i,Np1);
                du_dxi_nb = du_dxi(:,:,nb)+rxi_lift(:,:,nb);
                du_deta_nb = du_deta(:,:,nb)+reta_lift(:,:,nb);
                dv_dxi_nb = dv_dxi(:,:,nb)+ryi_lift(:,:,nb);
                dv_deta_nb = dv_deta(:,:,nb)+ryeta_lift(:,:,nb);
                xi_xn = mesh.xi_x(:,:,nb); xi_yn = mesh.xi_y(:,:,nb);
                eta_xn = mesh.eta_x(:,:,nb); eta_yn = mesh.eta_y(:,:,nb);
                ux_R = xi_xn(i,1)*du_dxi_nb(i,1)+eta_xn(i,1)*du_deta_nb(i,1);
                uy_R = xi_yn(i,1)*du_dxi_nb(i,1)+eta_yn(i,1)*du_deta_nb(i,1);
                vx_R = xi_xn(i,1)*dv_dxi_nb(i,1)+eta_xn(i,1)*dv_deta_nb(i,1);
                vy_R = xi_yn(i,1)*dv_dxi_nb(i,1)+eta_yn(i,1)*dv_deta_nb(i,1);
                ux_avg = 0.5*(ux_L+ux_R); uy_avg = 0.5*(uy_L+uy_R);
                vx_avg = 0.5*(vx_L+vx_R); vy_avg = 0.5*(vy_L+vy_R);
                Fvnx = (mu*(2*ux_avg-(2.0/3.0)*(ux_avg+vy_avg)))*nx ...
                      + mu*(uy_avg+vx_avg)*ny;
                Fvny = mu*(ux_avg+vy_avg)*nx ...
                      + (mu*(2*vy_avg-(2.0/3.0)*(ux_avg+vy_avg)))*ny;
                Res(i,Np1,2) = Res(i,Np1,2) - sJ*Fvnx/w(Np1);
                Res(i,Np1,3) = Res(i,Np1,3) - sJ*Fvny/w(Np1);
            end
        end
    else
        % Farfield boundary
        for i = 1:Np1
            sJ = sqrt(y_xi_f(i)^2 + x_xi_f(i)^2);
            if sJ < 1e-30, continue; end
            nx = -y_xi_f(i)/sJ;  ny = x_xi_f(i)/sJ;
            QL4 = squeeze(Qe(i,Np1,:));
            QR4 = farfield_bc(QL4, nx, ny, Mach, alpha_deg, gamma);
            lam = max(lam_e, ws(QR4,gamma));
            Fn_L = nflux(QL4,nx,ny,gamma);
            Fn_R = nflux(QR4,nx,ny,gamma);
            Fn_num = 0.5*(Fn_L+Fn_R) - 0.5*lam*(QR4(:)-QL4(:));
            for k=1:4
                Res(i,Np1,k) = Res(i,Np1,k) + sJ*(Fn_num(k)-Fn_L(k))/w(Np1);
            end
            % Viscous: far-field — zero viscous flux (freestream gradient ~0)
            % Residual already correct without viscous surface term
        end
    end

    % ── FACE 4: eta = -1 (j=1), radial inward ──────────────────────────
    y_xi_f = D * ye(:,1);   x_xi_f = D * xe(:,1);
    if ej > 1
        nb  = ei + (ej-2)*nxi_c;
        Qnb = squeeze(U(nb,:,:,:));
        for i = 1:Np1
            sJ = sqrt(y_xi_f(i)^2 + x_xi_f(i)^2);
            if sJ < 1e-30, continue; end
            nx = -y_xi_f(i)/sJ;  ny = x_xi_f(i)/sJ;
            % Left face (j=1) => -sJ/w(1) * (Fnum - F_right)
            Q_left  = squeeze(Qnb(i,Np1,:));
            Q_right = squeeze(Qe(i,1,:));
            lam = max(lam_e, ws(Q_left,gamma));
            Fn_L = nflux(Q_left,nx,ny,gamma);
            Fn_R = nflux(Q_right,nx,ny,gamma);
            Fn_num = 0.5*(Fn_L+Fn_R) - 0.5*lam*(Q_right(:)-Q_left(:));
            for k=1:4
                Res(i,1,k) = Res(i,1,k) - sJ*(Fn_num(k)-Fn_R(k))/w(1);
            end
            if mu > 0
                ux_R = du_dx(i,1); uy_R = du_dy(i,1);
                vx_R = dv_dx(i,1); vy_R = dv_dy(i,1);
                du_dxi_nb = du_dxi(:,:,nb)+rxi_lift(:,:,nb);
                du_deta_nb = du_deta(:,:,nb)+reta_lift(:,:,nb);
                dv_dxi_nb = dv_dxi(:,:,nb)+ryi_lift(:,:,nb);
                dv_deta_nb = dv_deta(:,:,nb)+ryeta_lift(:,:,nb);
                xi_xn = mesh.xi_x(:,:,nb); xi_yn = mesh.xi_y(:,:,nb);
                eta_xn = mesh.eta_x(:,:,nb); eta_yn = mesh.eta_y(:,:,nb);
                ux_L = xi_xn(i,Np1)*du_dxi_nb(i,Np1)+eta_xn(i,Np1)*du_deta_nb(i,Np1);
                uy_L = xi_yn(i,Np1)*du_dxi_nb(i,Np1)+eta_yn(i,Np1)*du_deta_nb(i,Np1);
                vx_L = xi_xn(i,Np1)*dv_dxi_nb(i,Np1)+eta_xn(i,Np1)*dv_deta_nb(i,Np1);
                vy_L = xi_yn(i,Np1)*dv_dxi_nb(i,Np1)+eta_yn(i,Np1)*dv_deta_nb(i,Np1);
                ux_avg = 0.5*(ux_L+ux_R); uy_avg = 0.5*(uy_L+uy_R);
                vx_avg = 0.5*(vx_L+vx_R); vy_avg = 0.5*(vy_L+vy_R);
                Fvnx = (mu*(2*ux_avg-(2.0/3.0)*(ux_avg+vy_avg)))*nx ...
                      + mu*(uy_avg+vx_avg)*ny;
                Fvny = mu*(ux_avg+vy_avg)*nx ...
                      + (mu*(2*vy_avg-(2.0/3.0)*(ux_avg+vy_avg)))*ny;
                Res(i,1,2) = Res(i,1,2) + sJ*Fvnx/w(1);
                Res(i,1,3) = Res(i,1,3) + sJ*Fvny/w(1);
            end
        end
    else
        % Wall boundary (ej=1): no-slip via ghost cell
        for i = 1:Np1
            sJ = sqrt(y_xi_f(i)^2 + x_xi_f(i)^2);
            if sJ < 1e-30, continue; end
            nx = -y_xi_f(i)/sJ;  ny = x_xi_f(i)/sJ;
            Q_int4   = squeeze(Qe(i,1,:));
            Q_ghost4 = wall_bc(Q_int4, gamma);
            lam = max(lam_e, ws(Q_ghost4,gamma));
            Fn_ghost = nflux(Q_ghost4,nx,ny,gamma);
            Fn_int   = nflux(Q_int4,nx,ny,gamma);
            Fn_num   = 0.5*(Fn_ghost+Fn_int) - 0.5*lam*(Q_int4(:)-Q_ghost4(:));
            % Left face (j=1) => -sJ/w(1) * (Fnum - F_right)
            % F_right = Fn_int (flux with interior state = the "right" state)
            for k=1:4
                Res(i,1,k) = Res(i,1,k) - sJ*(Fn_num(k)-Fn_int(k))/w(1);
            end
            % Viscous wall: no-slip enforced — viscous flux uses
            % u_wall=0 at face. Surface flux with averaged gradient
            % where ghost side gradient trace = 0 (for no-slip)
            if mu > 0
                % Interior gradient at wall
                ux_int = du_dx(i,1); uy_int = du_dy(i,1);
                vx_int = dv_dx(i,1); vy_int = dv_dy(i,1);

                % At no-slip wall, the common velocity gradient at the face
                % is computed by averaging interior (with lifting) and ghost.
                % Ghost velocity = 0 => ghost gradient = 0 at wall
                % Common: 0.5 * interior
                ux_avg = 0.5*ux_int; uy_avg = 0.5*uy_int;
                vx_avg = 0.5*vx_int; vy_avg = 0.5*vy_int;

                % Wall shear stress = mu * du/dn (physical normal derivative)
                % tau_wall . n
                Fvnx = (mu*(2*ux_avg-(2.0/3.0)*(ux_avg+vy_avg)))*nx ...
                      + mu*(uy_avg+vx_avg)*ny;
                Fvny = mu*(ux_avg+vy_avg)*nx ...
                      + (mu*(2*vy_avg-(2.0/3.0)*(ux_avg+vy_avg)))*ny;
                % Left face: +sJ/w(1) (viscous flux leaving domain at wall)
                Res(i,1,2) = Res(i,1,2) + sJ*Fvnx/w(1);
                Res(i,1,3) = Res(i,1,3) + sJ*Fvny/w(1);

                % Viscous energy: adiabatic wall (dT/dn = 0)
                % theta_n at wall = 0 (no heat flux)
                % No energy surface term needed
            end
        end
    end

    % ── dQ/dt = -Res / J ────────────────────────────────────────────────
    for jj = 1:Np1
        for ii = 1:Np1
            dUdt(e,ii,jj,:) = -squeeze(Res(ii,jj,:)) / Je(ii,jj);
        end
    end
end
end

% ── Helpers ──────────────────────────────────────────────────────────────

function lam = ws(Q4, gamma)
    Q4  = Q4(:);
    rho = max(Q4(1), 1e-12);
    u   = Q4(2)/rho;  v = Q4(3)/rho;  E = Q4(4)/rho;
    p   = max((gamma-1)*rho*(E - 0.5*(u^2+v^2)), 1e-12);
    lam = sqrt(u^2+v^2) + sqrt(gamma*p/rho);
end

function Fn = nflux(Q4, nx, ny, gamma)
    Q4  = Q4(:);
    rho = max(Q4(1), 1e-12);
    u = Q4(2)/rho;  v = Q4(3)/rho;
    E = Q4(4);
    p = max((gamma-1)*(E - 0.5*rho*(u^2+v^2)), 1e-12);
    Vn = u*nx + v*ny;
    Fn = [rho*Vn; rho*u*Vn + p*nx; rho*v*Vn + p*ny; (E+p)*Vn];
end
