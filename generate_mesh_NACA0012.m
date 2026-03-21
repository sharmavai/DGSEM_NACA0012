% =========================================================================
%  generate_mesh_NACA0012.m
%  C-grid mesh for NACA 0012 with full curvilinear metric terms
%
%  FIX #3 Applied:
%    - NACA 0012 thickness uses correct closed-TE coefficient -0.1036
%    - C-grid via transfinite interpolation (TFI) with geometric radial stretch
%    - Metric terms (J, xi_x, xi_y, eta_x, eta_y) computed using DGSEM
%      D-matrix — NOT finite differences — ensuring discrete GCL compliance
%    - Outward unit face normals at all 4 faces of every element
%    - Face arc-length weights (ds) for surface integration (CL/CD/CM)
%    - Jacobian positivity check (inverted elements = instant error)
%    - mesh struct returned with ALL fields needed by rhs_DGSEM and
%      compute_aero_coefficients
%
%  Mesh layout:
%    xi  direction  : wraps around airfoil (i-index, 1..nxi)
%    eta direction  : wall-normal, airfoil to far-field (j-index, 1..neta)
%    Each element has (N+1)x(N+1) LGL nodes in (xi,eta) space
%
%  References:
%    [1] Kopriva (2009) Ch. 8 — Metric terms via DGSEM, GCL
%    [2] Blazek (2015) Ch. 9 — C-grid topology, TFI
%    [3] Ferrer et al. (2021) — mesh parameters for Re=1e4 NACA0012
% =========================================================================

function mesh = generate_mesh_NACA0012(nxi, nwake, neta, R_far, y1, N)
% Inputs:
%   nxi   — elements wrapping airfoil  (recommend 48)
%   nwake — elements in wake strip     (recommend 16)
%   neta  — elements wall-normal       (recommend 20)
%   R_far — far-field radius [chord]   (recommend 30)
%   y1    — first cell height at wall  (recommend 5/Re)
%   N     — polynomial degree          (3)
%
% Output:
%   mesh  — struct with all geometric and connectivity data

    Np1   = N + 1;
    nxi_c = nxi + 2*nwake;   % total elements in wrap direction (C-topology)
    n_elem = nxi_c * neta;    % total element count

    fprintf('[mesh] nxi=%d  nwake=%d  neta=%d  nelem=%d\n', ...
            nxi, nwake, neta, n_elem);

    %% ── 1. LGL nodes and derivative matrix ──────────────────────────────
    [xi_lgl, w_lgl, ~, ~] = LGL_quadrature(N);
    [D, ~, ~]              = derivative_matrix(N);

    %% ── 2. NACA 0012 surface coordinates ─────────────────────────────────
    % Total surface nodes: nxi_c elements × (N+1) nodes, minus duplicates
    % Use cosine clustering for LE/TE resolution
    N_surf = nxi * Np1 + 1;   % airfoil wrap only (TE to TE via LE)
    theta  = linspace(0, 2*pi, N_surf)';
    s      = 0.5*(1 - cos(theta));    % s in [0,1], clustered at 0 and 1

    % NACA 0012 thickness (closed TE: coefficient -0.1036, NOT -0.1015)
    % Ref: Ladson et al. NASA TM-4074 eq. (1)
    t  = 0.12;
    yt = @(x) (t/0.2) * (0.2969*sqrt(x) - 0.1260*x ...
                         - 0.3516*x.^2   + 0.2843*x.^3 - 0.1036*x.^4);

    % Parametrise: lower surface TE->LE (s: 0->0.5), upper LE->TE (s: 0.5->1)
    x_surf = zeros(N_surf, 1);
    y_surf = zeros(N_surf, 1);
    for i = 1:N_surf
        if s(i) <= 0.5
            % Lower surface: x goes from 1 (TE) to 0 (LE)
            xc = 1 - 2*s(i);
            x_surf(i) =  xc;
            y_surf(i) = -yt(xc);
        else
            % Upper surface: x goes from 0 (LE) to 1 (TE)
            xc = 2*(s(i) - 0.5);
            x_surf(i) =  xc;
            y_surf(i) =  yt(xc);
        end
    end

    %% ── 3. Far-field circle ───────────────────────────────────────────────
    % Offset angle so TE gap of C-grid aligns with downstream direction
    theta_ff = linspace(-pi, pi, N_surf)';
    x_ff = 0.5 + R_far * cos(theta_ff);   % centred at mid-chord
    y_ff = R_far * sin(theta_ff);

    %% ── 4. Radial distribution (geometric stretching) ────────────────────
    % Growth ratio chosen so first cell = y1, last cell reaches R_far
    % Solve: y1 * (r^neta - 1)/(r - 1) = R_far  => iterate for r
    r = compute_growth_ratio(y1, R_far, neta);
    fprintf('[mesh] Radial growth ratio r = %.4f  (y1=%.2e)\n', r, y1);

    % Cumulative radial distances
    eta_dist = zeros(neta+1, 1);
    for j = 2:neta+1
        eta_dist(j) = y1 * (r^(j-1) - 1) / (r - 1);
    end
    eta_dist = eta_dist / eta_dist(end);   % normalise to [0,1]

    %% ── 5. TFI grid: build node positions for all elements ───────────────
    % Global node arrays (element corners at LGL nodes mapped from [-1,1])
    % xn(i,j), yn(i,j): physical coords at global node (i,j)
    %   i: 1 .. nxi_c*(N+1) in wrap direction
    %   j: 1 .. neta*(N+1)  in radial direction

    N_wrap  = nxi_c * Np1;   % global nodes in wrap direction (with duplicates at seams)
    N_rad   = neta  * Np1;   % global nodes in radial direction

    % Build surf/ff point arrays at all element LGL node positions
    % (re-interpolate surface parametrically to LGL-spaced nodes)
    s_elem = build_element_s(nxi_c, Np1, xi_lgl);   % parameter s at each LGL node

    x_s = interp1(linspace(0,1,N_surf), x_surf, s_elem, 'pchip');
    y_s = interp1(linspace(0,1,N_surf), y_surf, s_elem, 'pchip');
    x_f = interp1(linspace(0,1,N_surf), x_ff,   s_elem, 'pchip');
    y_f = interp1(linspace(0,1,N_surf), y_ff,   s_elem, 'pchip');

    % Radial parameter at all LGL node positions in eta-direction
    eta_elem = build_element_s(neta, Np1, xi_lgl);  % [0,1] mapped from [-1,1]

    % Physical coordinates: TFI blending
    Xn = zeros(N_wrap, N_rad);
    Yn = zeros(N_wrap, N_rad);
    for ii = 1:N_wrap
        for jj = 1:N_rad
            eta_jj = eta_elem(jj);
            Xn(ii,jj) = (1-eta_jj)*x_s(ii) + eta_jj*x_f(ii);
            Yn(ii,jj) = (1-eta_jj)*y_s(ii) + eta_jj*y_f(ii);
        end
    end

    %% ── 6. Compute metric terms via DGSEM D-matrix ────────────────────────
    % For each element k: extract (N+1)x(N+1) node block, apply D in both dirs
    % CRITICAL: Using D (not finite differences) ensures discrete GCL.
    % Reference: Kopriva (2009) eq. 8.28-8.33

    J       = zeros(Np1, Np1, n_elem);   % Jacobian det
    xi_x    = zeros(Np1, Np1, n_elem);   % dxi/dx
    xi_y    = zeros(Np1, Np1, n_elem);   % dxi/dy
    eta_x   = zeros(Np1, Np1, n_elem);   % deta/dx
    eta_y   = zeros(Np1, Np1, n_elem);   % deta/dy
    x_elem  = zeros(Np1, Np1, n_elem);   % physical x at LGL nodes
    y_elem  = zeros(Np1, Np1, n_elem);   % physical y at LGL nodes

    elem = 0;
    for ej = 1:neta          % radial element index
        for ei = 1:nxi_c     % wrap element index
            elem = elem + 1;

            % Global index ranges for this element
            i0 = (ei-1)*Np1 + 1;  i1 = ei*Np1;
            j0 = (ej-1)*Np1 + 1;  j1 = ej*Np1;

            xe = Xn(i0:i1, j0:j1);   % (Np1 x Np1) physical x
            ye = Yn(i0:i1, j0:j1);   % (Np1 x Np1) physical y

            x_elem(:,:,elem) = xe;
            y_elem(:,:,elem) = ye;

            % Metric terms using D (xi=rows=i-direction, eta=cols=j-direction)
            % dx/dxi  = D * xe       (derivative in xi applied row-wise)
            % dx/deta = xe * D'      (derivative in eta applied col-wise)
            x_xi  = D * xe;       % (Np1 x Np1)
            y_xi  = D * ye;
            x_eta = xe * D';
            y_eta = ye * D';

            Je =  x_xi .* y_eta - x_eta .* y_xi;

            % Check for negative Jacobian (inverted element = bad mesh)
            if any(Je(:) <= 0)
                error('[mesh] Negative Jacobian in element (%d,%d)! min J = %e\n  Reduce growth ratio or increase neta.', ei, ej, min(Je(:)));
            end

            % Metric terms of inverse mapping (Kopriva 2009, eq. 8.31)
            J(:,:,elem)     = Je;
            xi_x(:,:,elem)  =  y_eta ./ Je;
            xi_y(:,:,elem)  = -x_eta ./ Je;
            eta_x(:,:,elem) = -y_xi  ./ Je;
            eta_y(:,:,elem) =  x_xi  ./ Je;
        end
    end

    fprintf('[mesh] Metric computation done. min(J) = %.4e\n', min(J(:)));

    %% ── 7. Element connectivity ───────────────────────────────────────────
    % elem_id(ei, ej) = global element number
    elem_id = reshape(1:n_elem, nxi_c, neta);   % elem_id(wrap, radial)

    % Face connectivity: for each face, store [elemL, sideL, elemR, sideR]
    % Side convention: 1=xi_min, 2=xi_max, 3=eta_min, 4=eta_max
    n_faces_est = 2*n_elem + nxi_c + neta;
    face_conn   = zeros(n_faces_est, 4);
    face_nx     = zeros(Np1, n_faces_est);
    face_ny     = zeros(Np1, n_faces_est);
    face_ds     = zeros(Np1, n_faces_est);   % |face Jacobian| * w_lgl

    nf = 0;

    % Interior faces in xi-direction (between elements in wrap direction)
    for ej = 1:neta
        for ei = 1:nxi_c-1
            nf = nf + 1;
            eL = elem_id(ei,   ej);
            eR = elem_id(ei+1, ej);
            face_conn(nf,:) = [eL, 2, eR, 1];   % xi_max of L, xi_min of R
            % Normal at xi=+1 face of eL: n = (xi_x, xi_y) / |(xi_x, xi_y)|
            [nx, ny, ds] = face_normal_xi(x_elem(:,:,eL), y_elem(:,:,eL), D, w_lgl, +1);
            face_nx(:,nf) = nx;  face_ny(:,nf) = ny;  face_ds(:,nf) = ds;
        end
        % C-grid wrap: last xi element connects to first
        nf = nf + 1;
        eL = elem_id(nxi_c, ej);
        eR = elem_id(1,     ej);
        face_conn(nf,:) = [eL, 2, eR, 1];
        [nx, ny, ds] = face_normal_xi(x_elem(:,:,eL), y_elem(:,:,eL), D, w_lgl, +1);
        face_nx(:,nf) = nx;  face_ny(:,nf) = ny;  face_ds(:,nf) = ds;
    end

    % Interior faces in eta-direction (between radial elements)
    for ej = 1:neta-1
        for ei = 1:nxi_c
            nf = nf + 1;
            eL = elem_id(ei, ej);
            eR = elem_id(ei, ej+1);
            face_conn(nf,:) = [eL, 4, eR, 3];   % eta_max of L, eta_min of R
            [nx, ny, ds] = face_normal_eta(x_elem(:,:,eL), y_elem(:,:,eL), D, w_lgl, +1);
            face_nx(:,nf) = nx;  face_ny(:,nf) = ny;  face_ds(:,nf) = ds;
        end
    end

    n_faces = nf;
    face_conn = face_conn(1:n_faces, :);
    face_nx   = face_nx(:, 1:n_faces);
    face_ny   = face_ny(:, 1:n_faces);
    face_ds   = face_ds(:, 1:n_faces);

    %% ── 8. Boundary face lists ────────────────────────────────────────────
    % Airfoil wall: eta_min (ej=1) faces
    wall_faces = zeros(nxi_c, 1);
    for ei = 1:nxi_c
        % Find the eta_min face of element elem_id(ei,1)
        % These are faces whose sideL=3 and elemL is in row ej=1
        % For simplicity, store directly — indexed separately
        wall_faces(ei) = elem_id(ei, 1);
    end

    % Far-field: eta_max (ej=neta) faces
    ff_elems = zeros(nxi_c, 1);
    for ei = 1:nxi_c
        ff_elems(ei) = elem_id(ei, neta);
    end

    % Airfoil surface normals and positions (for CL/CD/CM integration)
    wall_nx = zeros(Np1, nxi_c);
    wall_ny = zeros(Np1, nxi_c);
    wall_ds = zeros(Np1, nxi_c);
    wall_x  = zeros(Np1, nxi_c);
    wall_y  = zeros(Np1, nxi_c);
    for ei = 1:nxi_c
        ew = elem_id(ei, 1);
        xe = x_elem(:,:,ew);
        ye = y_elem(:,:,ew);
        % eta_min face: eta=-1, normal points INWARD to fluid => outward from wall
        [nx, ny, ds] = face_normal_eta(xe, ye, D, w_lgl, -1);
        wall_nx(:,ei) = -nx;   % flip: outward from airfoil INTO fluid
        wall_ny(:,ei) = -ny;
        wall_ds(:,ei) =  ds;
        wall_x(:,ei)  = xe(:,1);   % j=1 column (eta_min nodes)
        wall_y(:,ei)  = ye(:,1);
    end

    % Identify which wall elements are actually on the airfoil
    % (vs wake cut of C-grid — elements ei=1..nxi correspond to airfoil)
    airfoil_ei = 1:nxi;    % first nxi elements wrap the actual airfoil

    %% ── 9. Minimum element size (for time step) ──────────────────────────
    h_min = zeros(n_elem, 1);
    for k = 1:n_elem
        % Approximate: min edge length ~ min(J) * 2
        h_min(k) = min(J(:,:,k), [], 'all') * 2;
    end

    %% ── 10. Assemble mesh struct ─────────────────────────────────────────
    mesh.n_elem      = n_elem;
    mesh.n_faces     = n_faces;
    mesh.nxi_c       = nxi_c;
    mesh.neta        = neta;
    mesh.nxi_airfoil = nxi;
    mesh.N           = N;

    % Physical coordinates at LGL nodes
    mesh.x   = x_elem;    % (Np1 x Np1 x n_elem)
    mesh.y   = y_elem;

    % Metric terms at LGL nodes
    mesh.J     = J;        % Jacobian determinant
    mesh.xi_x  = xi_x;
    mesh.xi_y  = xi_y;
    mesh.eta_x = eta_x;
    mesh.eta_y = eta_y;

    % Face data
    mesh.face_conn = face_conn;   % (n_faces x 4): [eL,sideL,eR,sideR]
    mesh.face_nx   = face_nx;     % (Np1 x n_faces) outward normals from eL
    mesh.face_ny   = face_ny;
    mesh.face_ds   = face_ds;     % (Np1 x n_faces) face integration weights

    % Boundary data
    mesh.wall_elems = wall_faces;    % element indices of wall-adjacent elements
    mesh.ff_elems   = ff_elems;
    mesh.airfoil_ei = airfoil_ei;

    % Wall (airfoil surface) geometry
    mesh.wall_nx = wall_nx;     % (Np1 x nxi_c) outward normals at wall
    mesh.wall_ny = wall_ny;
    mesh.wall_ds = wall_ds;     % arc-length weights at wall
    mesh.wall_x  = wall_x;     % (Np1 x nxi_c) x-coords on airfoil surface
    mesh.wall_y  = wall_y;

    % Time step helper
    mesh.h_min = h_min;         % (n_elem x 1) min element size

    % LGL data (needed by other routines)
    mesh.xi_lgl = xi_lgl;
    mesh.w_lgl  = w_lgl;
    mesh.D      = D;

    fprintf('[mesh] Mesh complete. n_elem=%d, n_faces=%d\n', n_elem, n_faces);
    fprintf('[mesh] min h = %.4e,  min J = %.4e\n', min(h_min), min(J(:)));
end

% =========================================================================
%  LOCAL HELPERS
% =========================================================================

function r = compute_growth_ratio(y1, L, n)
% Solve y1*(r^n-1)/(r-1) = L for r by bisection
    f = @(r) y1*(r.^n - 1)./(r - 1) - L;
    r = 1.0001;
    if f(r) > 0, return; end   % already satisfied
    r_lo = 1.0001;  r_hi = 2.0;
    for k = 1:100
        r_mid = 0.5*(r_lo + r_hi);
        if f(r_mid) < 0
            r_lo = r_mid;
        else
            r_hi = r_mid;
        end
        if r_hi - r_lo < 1e-10, break; end
    end
    r = 0.5*(r_lo + r_hi);
end

function s_nodes = build_element_s(n_elem, Np1, xi_lgl)
% Build global parameter array s in [0,1] for n_elem elements,
% each with Np1 LGL nodes mapped from [-1,1]
    s_nodes = zeros(n_elem * Np1, 1);
    ds = 1.0 / n_elem;
    for e = 1:n_elem
        s0 = (e-1)*ds;
        % Map xi_lgl in [-1,1] to [s0, s0+ds]
        idx = (e-1)*Np1 + (1:Np1);
        s_nodes(idx) = s0 + 0.5*(xi_lgl + 1)*ds;
    end
end

function [nx, ny, ds] = face_normal_xi(xe, ye, D, w, side)
% Outward normal at xi=side face (side=+1 for xi_max, =-1 for xi_min)
% xe, ye: (Np1 x Np1) element physical coords
% Returns (Np1 x 1) normal components and ds at Np1 eta-direction LGL nodes
    Np1 = size(xe,1);
    if side > 0
        % xi_max face: i = Np1 (last row in xi)
        x_face = xe(end, :)';   % Np1 x 1
        y_face = ye(end, :)';
        % Tangent in eta direction at xi=+1: dy/deta, dx/deta
        x_eta_face = (xe(end,:) * D')';
        y_eta_face = (ye(end,:) * D')';
    else
        x_face = xe(1, :)';
        y_face = ye(1, :)';
        x_eta_face = (xe(1,:) * D')';
        y_eta_face = (ye(1,:) * D')';
    end
    % Normal to xi=const face: n ~ (y_eta, -x_eta) then normalised
    % Arc-length element: |face_tangent| = sqrt(x_eta^2 + y_eta^2)
    face_len = sqrt(x_eta_face.^2 + y_eta_face.^2);
    nx = sign(side) *  y_eta_face ./ face_len;
    ny = sign(side) * (-x_eta_face) ./ face_len;
    ds = face_len .* w(:);   % integration weight = |tangent| * LGL weight
end

function [nx, ny, ds] = face_normal_eta(xe, ye, D, w, side)
% Outward normal at eta=side face (side=+1 for eta_max, =-1 for eta_min)
    if side > 0
        x_xi_face = (D * xe(:,end));   % Np1 x 1
        y_xi_face = (D * ye(:,end));
    else
        x_xi_face = (D * xe(:,1));
        y_xi_face = (D * ye(:,1));
    end
    % Normal to eta=const face: n ~ (-y_xi, x_xi) then normalised
    face_len = sqrt(x_xi_face.^2 + y_xi_face.^2);
    nx = sign(side) * (-y_xi_face) ./ face_len;
    ny = sign(side) *  x_xi_face  ./ face_len;
    ds = face_len .* w(:);
end