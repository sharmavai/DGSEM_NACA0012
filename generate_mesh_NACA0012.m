% =========================================================================
%  generate_mesh_NACA0012.m
%  C-grid mesh for NACA 0012 with full curvilinear metric terms
%
%  PREVIOUS BUG:
%    N_surf was sized for nxi (airfoil only) but s_elem covered nxi_c
%    (airfoil + 2*nwake).  The wake portion of s_elem extrapolated beyond
%    the airfoil data => wrong TFI coordinates => folded elements.
%
%  FIX:
%    Inner boundary (eta=0) is built in three explicit sections:
%      Section A: lower wake  — straight line from downstream to TE-lower
%      Section B: airfoil     — NACA 0012 surface (lower TE -> upper TE)
%      Section C: upper wake  — straight line from TE-upper to downstream
%    Outer boundary (eta=1) is parametrized with the SAME section layout.
%    TFI blending is then consistent across all nxi_c elements.
%
%  References:
%    [1] Kopriva (2009) Ch. 8 — metric terms, GCL
%    [2] Blazek (2015) Ch. 9 — C-grid topology
% =========================================================================

function mesh = generate_mesh_NACA0012(nxi, nwake, neta, R_far, y1, N)

    Np1   = N + 1;
    nxi_c = nxi + 2*nwake;      % total wrap elements
    n_elem = nxi_c * neta;

    fprintf('[mesh] nxi=%d  nwake=%d  neta=%d  nelem=%d\n', ...
            nxi, nwake, neta, n_elem);

    %% ── 1. LGL nodes and derivative matrix ──────────────────────────────
    [xi_lgl, w_lgl] = LGL_quadrature(N);
    [D, ~, ~]       = derivative_matrix(N);

    %% ── 2. NACA 0012 surface points (airfoil section only) ──────────────
    % Airfoil parametrized with nxi*Np1+1 points using COSINE CLUSTERING.
    % FIX: Previous code used uniform linspace for s_foil, producing
    %      uniform x-spacing — no LE/TE clustering despite the README
    %      promising cosine clustering.  This caused poor LE resolution
    %      (the NACA 0012 thickness has infinite slope at x=0).
    % CORRECTED: Use x = 0.5*(1 ± cos(theta)) for proper cosine
    %            distribution that clusters points near LE and TE.
    N_foil = nxi * Np1 + 1;
    t  = 0.12;
    yt = @(x) (t/0.2)*(0.2969*sqrt(x) - 0.1260*x ...
                       - 0.3516*x.^2   + 0.2843*x.^3 - 0.1036*x.^4);

    N_half = ceil(N_foil / 2);          % points on lower surface (incl. LE)
    N_up   = N_foil - N_half + 1;       % points on upper surface (incl. shared LE)

    % Lower surface: TE (x=1) -> LE (x=0) with cosine clustering
    theta_lo = linspace(0, pi, N_half)';
    x_lo = 0.5*(1 + cos(theta_lo));     % x: 1 -> 0, clustered at both ends

    % Upper surface: LE (x=0) -> TE (x=1) with cosine clustering
    theta_up = linspace(0, pi, N_up)';
    x_up = 0.5*(1 - cos(theta_up));     % x: 0 -> 1, clustered at both ends

    x_foil = zeros(N_foil, 1);
    y_foil = zeros(N_foil, 1);

    x_foil(1:N_half)     = x_lo;
    y_foil(1:N_half)     = -yt(x_lo);   % lower surface (y < 0)

    x_foil(N_half:N_foil) = x_up;       % LE shared at index N_half
    y_foil(N_half:N_foil) =  yt(x_up);  % upper surface (y > 0)

    %% ── 3. Full C-grid inner boundary (nxi_c*Np1 + 1 points) ─────────────
    % Layout (going around the C):
    %   A: lower wake  (nwake elements): x from x_ds to 1,  y = 0  (below TE)
    %   B: airfoil     (nxi   elements): NACA0012 surface
    %   C: upper wake  (nwake elements): x from 1 to x_ds,  y = 0  (above TE)
    % x_ds = downstream x-intercept of far-field circle
    x_ds = 1.0 + R_far;    % far downstream point on x-axis

    N_wake  = nwake * Np1 + 1;   % points per wake section (including endpoint)
    N_inner = nxi_c * Np1 + 1;   % total inner boundary points

    x_inner = zeros(N_inner, 1);
    y_inner = zeros(N_inner, 1);

    % Section A: lower wake (x_ds -> x=1, y=0)
    idxA = 1 : nwake*Np1 + 1;
    x_inner(idxA) = linspace(x_ds, 1.0, N_wake);
    y_inner(idxA) = 0.0;   % y=0 (wake cut, below chord)

    % Section B: airfoil — share endpoint with section A
    idxB = nwake*Np1 + 1 : (nwake+nxi)*Np1 + 1;
    x_inner(idxB) = x_foil;
    y_inner(idxB) = y_foil;

    % Section C: upper wake (x=1 -> x_ds, y=0) — share endpoint with section B
    idxC = (nwake+nxi)*Np1 + 1 : nxi_c*Np1 + 1;
    x_inner(idxC) = linspace(1.0, x_ds, N_wake);
    y_inner(idxC) = 0.0;   % y=0 (wake cut, above chord)

    %% ── 4. Outer boundary (far-field, matching inner layout) ─────────────
    % Parametrize far-field circle with same nxi_c*Np1+1 points.
    % theta=0: (1+R_far, 0) = downstream, matching x_ds.
    % Go counter-clockwise to match C-grid orientation.
    x_outer = zeros(N_inner, 1);
    y_outer = zeros(N_inner, 1);
    theta_ff = linspace(0, -2*pi, N_inner)';
    x_outer  = 0.5 + R_far*cos(theta_ff);   % centred at mid-chord (0.5,0)
    y_outer  = R_far*sin(theta_ff);

    %% ── 5. LGL node positions along inner/outer boundaries ───────────────
    % For each of the nxi_c elements, place Np1 LGL nodes using build_s
    s_all = build_element_s(nxi_c, Np1, xi_lgl);   % [0,1] Np1*nxi_c points

    s_param = linspace(0, 1, N_inner);
    xi_s  = interp1(s_param, x_inner, s_all, 'pchip');   % inner boundary LGL x
    yi_s  = interp1(s_param, y_inner, s_all, 'pchip');   % inner boundary LGL y
    xo_s  = interp1(s_param, x_outer, s_all, 'pchip');   % outer boundary LGL x
    yo_s  = interp1(s_param, y_outer, s_all, 'pchip');   % outer boundary LGL y

    %% ── 6. Radial stretching ──────────────────────────────────────────────
    r = compute_growth_ratio(y1, R_far, neta);
    if r > 1.30
        r = 1.25;
        fprintf('[mesh] WARNING: growth ratio capped at %.2f\n', r);
    end
    fprintf('[mesh] Radial growth ratio r = %.4f  (y1=%.2e)\n', r, y1);

    % Physical distances at element boundaries
    eta_phys = zeros(neta+1, 1);
    for j = 2:neta+1
        eta_phys(j) = y1 * (r^(j-1) - 1) / (r - 1);
    end
    eta_dist = eta_phys / eta_phys(end);   % normalise [0,1]

    % Stretched LGL positions within each radial element
    N_wrap = nxi_c * Np1;
    N_rad  = neta  * Np1;

    eta_lgl_global = zeros(N_rad, 1);
    for ej = 1:neta
        eta_lo = eta_dist(ej);
        eta_hi = eta_dist(ej+1);
        idx    = (ej-1)*Np1 + (1:Np1);
        eta_lgl_global(idx) = eta_lo + 0.5*(xi_lgl+1)*(eta_hi - eta_lo);
    end

    %% ── 7. TFI: build all node coordinates ───────────────────────────────
    Xn = zeros(N_wrap, N_rad);
    Yn = zeros(N_wrap, N_rad);
    for ii = 1:N_wrap
        for jj = 1:N_rad
            t = eta_lgl_global(jj);
            Xn(ii,jj) = (1-t)*xi_s(ii) + t*xo_s(ii);
            Yn(ii,jj) = (1-t)*yi_s(ii) + t*yo_s(ii);
        end
    end

    %% ── 8. Metric terms via DGSEM D-matrix ──────────────────────────────
    J     = zeros(Np1, Np1, n_elem);
    xi_x  = zeros(Np1, Np1, n_elem);
    xi_y  = zeros(Np1, Np1, n_elem);
    eta_x = zeros(Np1, Np1, n_elem);
    eta_y = zeros(Np1, Np1, n_elem);
    x_elem = zeros(Np1, Np1, n_elem);
    y_elem = zeros(Np1, Np1, n_elem);

    elem = 0;
    for ej = 1:neta
        for ei = 1:nxi_c
            elem = elem + 1;
            i0 = (ei-1)*Np1+1;  i1 = ei*Np1;
            j0 = (ej-1)*Np1+1;  j1 = ej*Np1;

            xe = Xn(i0:i1, j0:j1);
            ye = Yn(i0:i1, j0:j1);
            x_elem(:,:,elem) = xe;
            y_elem(:,:,elem) = ye;

            x_xi  = D*xe;      y_xi  = D*ye;
            x_eta = xe*D';     y_eta = ye*D';
            Je    = x_xi.*y_eta - x_eta.*y_xi;

            if any(Je(:) <= 0)
                fprintf('[mesh] Element (%d,%d): x range [%.4f %.4f] y range [%.4f %.4f]\n', ...
                        ei, ej, min(xe(:)), max(xe(:)), min(ye(:)), max(ye(:)));
                error('[mesh] Negative Jacobian in element (%d,%d)! min J = %e', ...
                      ei, ej, min(Je(:)));
            end

            J(:,:,elem)     = Je;
            xi_x(:,:,elem)  =  y_eta./Je;
            xi_y(:,:,elem)  = -x_eta./Je;
            eta_x(:,:,elem) = -y_xi ./Je;
            eta_y(:,:,elem) =  x_xi ./Je;
        end
    end

    fprintf('[mesh] Metric computation done. min(J) = %.4e\n', min(J(:)));

    %% ── 9. Face connectivity ─────────────────────────────────────────────
    elem_id = reshape(1:n_elem, nxi_c, neta);

    n_faces_est = 2*n_elem + nxi_c + neta;
    face_conn = zeros(n_faces_est, 4);
    face_nx   = zeros(Np1, n_faces_est);
    face_ny   = zeros(Np1, n_faces_est);
    face_ds   = zeros(Np1, n_faces_est);
    nf = 0;

    for ej = 1:neta
        for ei = 1:nxi_c-1
            nf = nf+1;
            eL = elem_id(ei,ej);  eR = elem_id(ei+1,ej);
            face_conn(nf,:) = [eL,2,eR,1];
            [nx,ny,ds] = face_normal_xi(x_elem(:,:,eL),y_elem(:,:,eL),D,w_lgl,+1);
            face_nx(:,nf)=nx; face_ny(:,nf)=ny; face_ds(:,nf)=ds;
        end
        nf = nf+1;
        eL = elem_id(nxi_c,ej);  eR = elem_id(1,ej);
        face_conn(nf,:) = [eL,2,eR,1];
        [nx,ny,ds] = face_normal_xi(x_elem(:,:,eL),y_elem(:,:,eL),D,w_lgl,+1);
        face_nx(:,nf)=nx; face_ny(:,nf)=ny; face_ds(:,nf)=ds;
    end

    for ej = 1:neta-1
        for ei = 1:nxi_c
            nf = nf+1;
            eL = elem_id(ei,ej);  eR = elem_id(ei,ej+1);
            face_conn(nf,:) = [eL,4,eR,3];
            [nx,ny,ds] = face_normal_eta(x_elem(:,:,eL),y_elem(:,:,eL),D,w_lgl,+1);
            face_nx(:,nf)=nx; face_ny(:,nf)=ny; face_ds(:,nf)=ds;
        end
    end

    n_faces   = nf;
    face_conn = face_conn(1:n_faces,:);
    face_nx   = face_nx(:,1:n_faces);
    face_ny   = face_ny(:,1:n_faces);
    face_ds   = face_ds(:,1:n_faces);

    %% ── 10. Boundary lists ───────────────────────────────────────────────
    wall_elems = zeros(nxi_c,1);
    ff_elems   = zeros(nxi_c,1);
    wall_nx = zeros(Np1,nxi_c);
    wall_ny = zeros(Np1,nxi_c);
    wall_ds = zeros(Np1,nxi_c);
    wall_x  = zeros(Np1,nxi_c);
    wall_y  = zeros(Np1,nxi_c);

    for ei = 1:nxi_c
        wall_elems(ei) = elem_id(ei,1);
        ff_elems(ei)   = elem_id(ei,neta);
        ew = elem_id(ei,1);
        xe = x_elem(:,:,ew);  ye = y_elem(:,:,ew);
        [nx,ny,ds] = face_normal_eta(xe,ye,D,w_lgl,-1);
        wall_nx(:,ei) = -nx;
        wall_ny(:,ei) = -ny;
        wall_ds(:,ei) =  ds;
        wall_x(:,ei)  = xe(:,1);
        wall_y(:,ei)  = ye(:,1);
    end

    airfoil_ei = nwake+1 : nwake+nxi;   % only true airfoil elements

    h_min = zeros(n_elem,1);
    for k = 1:n_elem
        h_min(k) = min(J(:,:,k),[],'all') * 2;
    end

    %% ── 11. Mesh struct ──────────────────────────────────────────────────
    mesh.n_elem      = n_elem;
    mesh.n_faces     = n_faces;
    mesh.nxi_c       = nxi_c;
    mesh.neta        = neta;
    mesh.nxi_airfoil = nxi;
    mesh.nwake       = nwake;
    mesh.N           = N;
    mesh.x           = x_elem;
    mesh.y           = y_elem;
    mesh.J           = J;
    mesh.xi_x        = xi_x;
    mesh.xi_y        = xi_y;
    mesh.eta_x       = eta_x;
    mesh.eta_y       = eta_y;
    mesh.face_conn   = face_conn;
    mesh.face_nx     = face_nx;
    mesh.face_ny     = face_ny;
    mesh.face_ds     = face_ds;
    mesh.wall_elems  = wall_elems;
    mesh.ff_elems    = ff_elems;
    mesh.airfoil_ei  = airfoil_ei;
    mesh.wall_nx     = wall_nx;
    mesh.wall_ny     = wall_ny;
    mesh.wall_ds     = wall_ds;
    mesh.wall_x      = wall_x;
    mesh.wall_y      = wall_y;
    mesh.h_min       = h_min;
    mesh.xi_lgl      = xi_lgl;
    mesh.w_lgl       = w_lgl;
    mesh.D           = D;

    fprintf('[mesh] Complete. n_elem=%d  min(h)=%.3e  min(J)=%.3e\n', ...
            n_elem, min(h_min), min(J(:)));
end

% =========================================================================
%  HELPERS
% =========================================================================

function r = compute_growth_ratio(y1, L, n)
    f = @(rv) y1*(rv.^n - 1)./(rv - 1) - L;
    if f(1.0001) >= 0, r = 1.0001; return; end
    r_lo = 1.0001;  r_hi = 5.0;
    for k = 1:200
        rm = 0.5*(r_lo+r_hi);
        if f(rm) < 0, r_lo = rm; else, r_hi = rm; end
        if r_hi-r_lo < 1e-10, break; end
    end
    r = 0.5*(r_lo+r_hi);
end

function s_nodes = build_element_s(n_elem, Np1, xi_lgl)
    ds = 1.0/n_elem;
    s_nodes = zeros(n_elem*Np1, 1);
    for e = 1:n_elem
        idx = (e-1)*Np1 + (1:Np1);
        s_nodes(idx) = (e-1)*ds + 0.5*(xi_lgl+1)*ds;
    end
end

function [nx,ny,ds] = face_normal_xi(xe, ye, D, w, side)
    if side > 0
        x_eta_f = (xe(end,:)*D')';  y_eta_f = (ye(end,:)*D')';
    else
        x_eta_f = (xe(1,:)*D')';    y_eta_f = (ye(1,:)*D')';
    end
    len = sqrt(x_eta_f.^2 + y_eta_f.^2) + 1e-300;
    nx = sign(side) *  y_eta_f./len;
    ny = sign(side) * -x_eta_f./len;
    ds = len .* w(:);
end

function [nx,ny,ds] = face_normal_eta(xe, ye, D, w, side)
    if side > 0
        x_xi_f = D*xe(:,end);  y_xi_f = D*ye(:,end);
    else
        x_xi_f = D*xe(:,1);    y_xi_f = D*ye(:,1);
    end
    len = sqrt(x_xi_f.^2 + y_xi_f.^2) + 1e-300;
    nx = sign(side) * -y_xi_f./len;
    ny = sign(side) *  x_xi_f./len;
    ds = len .* w(:);
end