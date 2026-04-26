% =========================================================================
%  rhs_DGSEM.m
%  Spatial residual for 2D compressible Euler equations via DGSEM
%  (strong-form SBP-SAT on curvilinear C-grid)
%
%  FIX Applied:
%    - Used mesh.nelem — corrected to mesh.n_elem.
%    - Used mesh.airfoil_ids / mesh.farfield_ids — corrected to
%      mesh.airfoil_ei / mesh.ff_elems.
%    - mesh.x was treated as (n_radial, n_circ, 4_corners).  mesh.x is
%      (Np1, Np1, n_elem).  Neighbor map rewritten using nxi_c/neta.
%    - xi-direction = circumferential, eta-direction = radial (matching
%      generate_mesh_NACA0012.m).  Wall BC is on eta=-1 (ej=1), farfield
%      is on eta=+1 (ej=neta).  Previous code had these on xi faces.
%    - Volume integral now uses curvilinear metric terms (J, xi_x, xi_y,
%      eta_x, eta_y) for contravariant fluxes — fixes GCL / free-stream
%      preservation on non-Cartesian meshes.
%    - Face fluxes use Lax-Friedrichs with proper face-normal computation
%      from element geometry.
%    - Calls wall_bc.m and farfield_bc.m for boundary ghost states.
%
%  Q layout: (n_elem, Np1, Np1, 4)
%    Q(e,i,j,1:4) = [rho, rho*u, rho*v, E]
%    i = xi-node (circumferential), j = eta-node (radial)
%
%  References:
%    [1] Kopriva (2009) Ch. 8 — DGSEM with metric terms
%    [2] Hesthaven & Warburton (2008) — SBP-SAT framework
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

dUdt = zeros(size(U));

for e = 1:n_elem
    ei = mod(e-1, nxi_c) + 1;
    ej = floor((e-1) / nxi_c) + 1;

    Qe = squeeze(U(e,:,:,:));   % (Np1, Np1, 4)
    xe = mesh.x(:,:,e);
    ye = mesh.y(:,:,e);
    Je = mesh.J(:,:,e);

    % ── Primitive variables ──────────────────────────────────────────
    rho = max(Qe(:,:,1), 1e-12);
    uu  = Qe(:,:,2) ./ rho;
    vv  = Qe(:,:,3) ./ rho;
    EE  = Qe(:,:,4) ./ rho;
    pp  = max((gamma-1)*rho.*(EE - 0.5*(uu.^2 + vv.^2)), 1e-12);
    HH  = EE + pp./rho;

    % ── Physical fluxes ──────────────────────────────────────────────
    F = zeros(Np1, Np1, 4);  G = zeros(Np1, Np1, 4);
    F(:,:,1) = rho.*uu;            G(:,:,1) = rho.*vv;
    F(:,:,2) = rho.*uu.*uu + pp;   G(:,:,2) = rho.*uu.*vv;
    F(:,:,3) = rho.*uu.*vv;        G(:,:,3) = rho.*vv.*vv + pp;
    F(:,:,4) = rho.*uu.*HH;        G(:,:,4) = rho.*vv.*HH;

    % ── Contravariant fluxes (metric terms) ──────────────────────────
    xi_xe  = mesh.xi_x(:,:,e);    xi_ye  = mesh.xi_y(:,:,e);
    eta_xe = mesh.eta_x(:,:,e);   eta_ye = mesh.eta_y(:,:,e);

    Fh = zeros(Np1, Np1, 4);
    Gh = zeros(Np1, Np1, 4);
    for k = 1:4
        Fh(:,:,k) = Je .* (xi_xe.*F(:,:,k)  + xi_ye.*G(:,:,k));
        Gh(:,:,k) = Je .* (eta_xe.*F(:,:,k) + eta_ye.*G(:,:,k));
    end

    % ── Volume integral: D_xi(Fh) + D_eta(Gh) ───────────────────────
    Res = zeros(Np1, Np1, 4);
    for k = 1:4
        Res(:,:,k) = D * Fh(:,:,k) + (D * Gh(:,:,k)')';
    end

    % Max wave speed in element
    cc    = sqrt(gamma * pp ./ rho);
    lam_e = max(sqrt(uu(:).^2 + vv(:).^2) + cc(:));

    % ── FACE 1: xi = +1 (i=Np1), circumferential +1 ─────────────────
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
        for k=1:4
            Res(Np1,j,k) = Res(Np1,j,k) + sJ*(Fn_num(k)-Fn_L(k))/w(Np1);
        end
    end

    % ── FACE 2: xi = -1 (i=1), circumferential -1 ───────────────────
    ei_m = mod(ei-2, nxi_c) + 1;
    nb   = ei_m + (ej-1)*nxi_c;
    Qnb  = squeeze(U(nb,:,:,:));
    y_eta_f = (ye(1,:) * D')';   x_eta_f = (xe(1,:) * D')';
    for j = 1:Np1
        sJ = sqrt(y_eta_f(j)^2 + x_eta_f(j)^2);
        if sJ < 1e-30, continue; end
        nx = y_eta_f(j)/sJ;  ny = -x_eta_f(j)/sJ;
        Q_left  = squeeze(Qnb(Np1,j,:));
        Q_right = squeeze(Qe(1,j,:));
        lam = max(lam_e, ws(Q_left,gamma));
        Fn_L = nflux(Q_left,nx,ny,gamma);
        Fn_R = nflux(Q_right,nx,ny,gamma);
        Fn_num = 0.5*(Fn_L+Fn_R) - 0.5*lam*(Q_right(:)-Q_left(:));
        for k=1:4
            Res(1,j,k) = Res(1,j,k) - sJ*(Fn_num(k)-Fn_R(k))/w(1);
        end
    end

    % ── FACE 3: eta = +1 (j=Np1), radial outward ────────────────────
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
            for k=1:4
                Res(i,Np1,k) = Res(i,Np1,k) + sJ*(Fn_num(k)-Fn_L(k))/w(Np1);
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
        end
    end

    % ── FACE 4: eta = -1 (j=1), radial inward ───────────────────────
    y_xi_f = D * ye(:,1);   x_xi_f = D * xe(:,1);
    if ej > 1
        nb  = ei + (ej-2)*nxi_c;
        Qnb = squeeze(U(nb,:,:,:));
        for i = 1:Np1
            sJ = sqrt(y_xi_f(i)^2 + x_xi_f(i)^2);
            if sJ < 1e-30, continue; end
            nx = -y_xi_f(i)/sJ;  ny = x_xi_f(i)/sJ;
            Q_left  = squeeze(Qnb(i,Np1,:));
            Q_right = squeeze(Qe(i,1,:));
            lam = max(lam_e, ws(Q_left,gamma));
            Fn_L = nflux(Q_left,nx,ny,gamma);
            Fn_R = nflux(Q_right,nx,ny,gamma);
            Fn_num = 0.5*(Fn_L+Fn_R) - 0.5*lam*(Q_right(:)-Q_left(:));
            for k=1:4
                Res(i,1,k) = Res(i,1,k) - sJ*(Fn_num(k)-Fn_R(k))/w(1);
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
            for k=1:4
                Res(i,1,k) = Res(i,1,k) - sJ*(Fn_num(k)-Fn_int(k))/w(1);
            end
        end
    end

    % ── dQ/dt = -Res / J ─────────────────────────────────────────────
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
