function mesh = generate_mesh_NACA0012(nxi, nwake, neta, R_far, y1, N)
% C-grid mesh for NACA 0012

    % Validate inputs
    if ~isscalar(nxi) || nxi < 8 || mod(nxi,2) ~= 0
        error('nxi must be an even integer >= 8');
    end
    if ~isscalar(nwake) || nwake < 0
        error('nwake must be a non-negative integer');
    end
    if ~isscalar(neta) || neta < 4
        error('neta must be an integer >= 4');
    end
    if ~isscalar(N) || N < 1
        error('N must be a positive integer');
    end
    if ~isscalar(R_far) || R_far <= 2
        error('R_far must be > 2');
    end
    if ~isscalar(y1) || y1 <= 0
        error('y1 must be > 0');
    end
    
    Np1 = N + 1;
    nxi_c = nxi + 2*nwake;
    n_elem = nxi_c * neta;
    
    fprintf('[mesh] nxi=%d  nwake=%d  neta=%d  nelem=%d  N=%d\n', ...
            nxi, nwake, neta, n_elem, N);
    
    [xi_lgl, w_lgl] = LGL_quadrature(N);
    [D, ~, ~] = derivative_matrix(N);
    
    % NACA 0012
    t = 0.12;
    yt = @(x) (t/0.2)*(0.2969*sqrt(x) - 0.1260*x - 0.3516*x.^2 + 0.2843*x.^3 - 0.1036*x.^4);
    
    % Airfoil points
    N_foil = nxi * Np1 + 1;
    N_half = ceil(N_foil / 2);
    N_up = N_foil - N_half + 1;
    
    theta_lo = linspace(0, pi, N_half)';
    x_lo = 0.5*(1 + cos(theta_lo));
    
    theta_up = linspace(0, pi, N_up)';
    x_up = 0.5*(1 - cos(theta_up));
    
    x_foil = zeros(N_foil, 1);
    y_foil = zeros(N_foil, 1);
    x_foil(1:N_half) = x_lo;
    y_foil(1:N_half) = -yt(x_lo);
    x_foil(N_half:N_foil) = x_up;
    y_foil(N_half:N_foil) = yt(x_up);
    
    % Inner boundary
    N_wake = nwake * Np1 + 1;
    x_ds = 0.5 + R_far;
    
    N_inner = nxi_c * Np1 + 1;
    x_inner = zeros(N_inner, 1);
    y_inner = zeros(N_inner, 1);
    
    % Lower wake
    idx = 1:(nwake*Np1+1);
    x_inner(idx) = linspace(x_ds, 1.0, N_wake)';
    y_inner(idx) = 0;
    
    % Airfoil
    idx = (nwake*Np1+1):((nwake+nxi)*Np1+1);
    x_inner(idx) = x_foil;
    y_inner(idx) = y_foil;
    
    % Upper wake
    idx = ((nwake+nxi)*Np1+1):N_inner;
    x_inner(idx) = linspace(1.0, x_ds, N_wake)';
    y_inner(idx) = 0;
    
    % Outer boundary
    theta_ff = linspace(0, -2*pi, N_inner)';
    x_outer = 0.5 + R_far*cos(theta_ff);
    y_outer = R_far*sin(theta_ff);
    
    % TFI
    s_param = linspace(0, 1, N_inner);
    s_all = build_element_s(nxi_c, Np1, xi_lgl);
    
    xi_s = interp1(s_param, x_inner, s_all, 'pchip');
    yi_s = interp1(s_param, y_inner, s_all, 'pchip');
    xo_s = interp1(s_param, x_outer, s_all, 'pchip');
    yo_s = interp1(s_param, y_outer, s_all, 'pchip');
    
    % Radial stretching
    r = compute_growth_ratio(y1, R_far, neta);
    fprintf('[mesh] Radial growth ratio r = %.4f\n', r);
    
    eta_phys = zeros(neta+1, 1);
    for j = 2:neta+1
        eta_phys(j) = y1 * (r^(j-1) - 1) / (r - 1);
    end
    eta_dist = eta_phys / eta_phys(end);
    
    N_wrap = nxi_c * Np1;
    N_rad = neta * Np1;
    eta_lgl_global = zeros(N_rad, 1);
    for ej = 1:neta
        eta_lo = eta_dist(ej);
        eta_hi = eta_dist(ej+1);
        idx = (ej-1)*Np1 + (1:Np1);
        eta_lgl_global(idx) = eta_lo + 0.5*(xi_lgl+1)*(eta_hi - eta_lo);
    end
    
    % Build coordinates
    Xn = zeros(N_wrap, N_rad);
    Yn = zeros(N_wrap, N_rad);
    for ii = 1:N_wrap
        for jj = 1:N_rad
            t = eta_lgl_global(jj);
            Xn(ii,jj) = (1-t)*xi_s(ii) + t*xo_s(ii);
            Yn(ii,jj) = (1-t)*yi_s(ii) + t*yo_s(ii);
        end
    end
    
    % Metrics
    J = zeros(Np1, Np1, n_elem);
    xi_x = zeros(Np1, Np1, n_elem);
    xi_y = zeros(Np1, Np1, n_elem);
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
            
            x_xi = D*xe;      y_xi = D*ye;
            x_eta = xe*D';    y_eta = ye*D';
            Je = x_xi.*y_eta - x_eta.*y_xi;
            
            if any(Je(:) <= 0)
                fprintf('[mesh] Element (%d,%d): min J = %e\n', ei, ej, min(Je(:)));
                error('[mesh] Negative Jacobian in element (%d,%d)!', ei, ej);
            end
            
            J(:,:,elem) = Je;
            xi_x(:,:,elem) = y_eta./Je;
            xi_y(:,:,elem) = -x_eta./Je;
            eta_x(:,:,elem) = -y_xi./Je;
            eta_y(:,:,elem) = x_xi./Je;
        end
    end
    
    fprintf('[mesh] Metric computation done. min(J) = %.4e\n', min(J(:)));
    
    % Face connectivity
    elem_id = reshape(1:n_elem, nxi_c, neta);
    
    n_faces_est = 2*n_elem + nxi_c + neta;
    face_conn = zeros(n_faces_est, 4);
    face_nx = zeros(Np1, n_faces_est);
    face_ny = zeros(Np1, n_faces_est);
    face_ds = zeros(Np1, n_faces_est);
    nf = 0;
    
    for ej = 1:neta
        for ei = 1:nxi_c-1
            nf = nf+1;
            eL = elem_id(ei,ej);
            eR = elem_id(ei+1,ej);
            face_conn(nf,:) = [eL, 2, eR, 1];
            [nx,ny,ds] = face_normal_xi(x_elem(:,:,eL), y_elem(:,:,eL), D, w_lgl, +1);
            face_nx(:,nf) = nx; face_ny(:,nf) = ny; face_ds(:,nf) = ds;
        end
    end
    
    for ej = 1:neta-1
        for ei = 1:nxi_c
            nf = nf+1;
            eL = elem_id(ei,ej);
            eR = elem_id(ei,ej+1);
            face_conn(nf,:) = [eL, 4, eR, 3];
            [nx,ny,ds] = face_normal_eta(x_elem(:,:,eL), y_elem(:,:,eL), D, w_lgl, +1);
            face_nx(:,nf) = nx; face_ny(:,nf) = ny; face_ds(:,nf) = ds;
        end
    end
    
    n_faces = nf;
    face_conn = face_conn(1:n_faces, :);
    face_nx = face_nx(:, 1:n_faces);
    face_ny = face_ny(:, 1:n_faces);
    face_ds = face_ds(:, 1:n_faces);
    
    % Boundaries
    airfoil_ei = nwake+1 : nwake+nxi;
    n_wall = length(airfoil_ei);
    
    wall_elems = zeros(n_wall, 1);
    wall_nx = zeros(Np1, n_wall);
    wall_ny = zeros(Np1, n_wall);
    wall_ds = zeros(Np1, n_wall);
    wall_x = zeros(Np1, n_wall);
    wall_y = zeros(Np1, n_wall);
    
    for idx = 1:n_wall
        ei = airfoil_ei(idx);
        e = elem_id(ei, 1);
        wall_elems(idx) = e;
        xe = x_elem(:,:,e);
        ye = y_elem(:,:,e);
        [nx,ny,ds] = face_normal_eta(xe, ye, D, w_lgl, -1);
        wall_nx(:,idx) = -nx;
        wall_ny(:,idx) = -ny;
        wall_ds(:,idx) = ds;
        wall_x(:,idx) = xe(:,1);
        wall_y(:,idx) = ye(:,1);
    end
    
    ff_elems = elem_id(:, neta)';
    
    % h_min
    h_min = zeros(n_elem, 1);
    for k = 1:n_elem
        xe = x_elem(:,:,k);
        ye = y_elem(:,:,k);
        
        dx_dxi = D * xe;
        dy_dxi = D * ye;
        dx_deta = xe * D';
        dy_deta = ye * D';
        
        edge_xi = sqrt(dx_dxi.^2 + dy_dxi.^2);
        edge_eta = sqrt(dx_deta.^2 + dy_deta.^2);
        
        h_min(k) = min([min(edge_xi(:)), min(edge_eta(:))]);
    end
    
    fprintf('[mesh] Complete. n_elem=%d  min(h)=%.3e  min(J)=%.3e\n', ...
            n_elem, min(h_min), min(J(:)));
    
    % Mesh struct
    mesh.n_elem = n_elem;
    mesh.n_faces = n_faces;
    mesh.nxi_c = nxi_c;
    mesh.neta = neta;
    mesh.nxi_airfoil = nxi;
    mesh.nwake = nwake;
    mesh.N = N;
    mesh.x = x_elem;
    mesh.y = y_elem;
    mesh.J = J;
    mesh.xi_x = xi_x;
    mesh.xi_y = xi_y;
    mesh.eta_x = eta_x;
    mesh.eta_y = eta_y;
    mesh.face_conn = face_conn;
    mesh.face_nx = face_nx;
    mesh.face_ny = face_ny;
    mesh.face_ds = face_ds;
    mesh.wall_elems = wall_elems;
    mesh.ff_elems = ff_elems;
    mesh.airfoil_ei = airfoil_ei;
    mesh.wall_nx = wall_nx;
    mesh.wall_ny = wall_ny;
    mesh.wall_ds = wall_ds;
    mesh.wall_x = wall_x;
    mesh.wall_y = wall_y;
    mesh.h_min = h_min;
    mesh.xi_lgl = xi_lgl;
    mesh.w_lgl = w_lgl;
    mesh.D = D;
end

function s_nodes = build_element_s(n_elem, Np1, xi_lgl)
    ds = 1.0 / n_elem;
    s_nodes = zeros(n_elem*Np1, 1);
    for e = 1:n_elem
        idx = (e-1)*Np1 + (1:Np1);
        s_nodes(idx) = (e-1)*ds + 0.5*(xi_lgl+1)*ds;
    end
end

function r = compute_growth_ratio(y1, L, n)
    f = @(rv) y1*(rv.^n - 1)./(rv - 1) - L;
    if f(1.0001) >= 0
        r = 1.0001;
        return;
    end
    r_lo = 1.0001;
    r_hi = 5.0;
    for k = 1:200
        rm = 0.5*(r_lo + r_hi);
        if f(rm) < 0
            r_lo = rm;
        else
            r_hi = rm;
        end
        if r_hi - r_lo < 1e-10
            break;
        end
    end
    r = 0.5*(r_lo + r_hi);
end

function [nx,ny,ds] = face_normal_xi(xe, ye, D, w, side)
    if side > 0
        x_eta_f = (xe(end,:)*D')';
        y_eta_f = (ye(end,:)*D')';
    else
        x_eta_f = (xe(1,:)*D')';
        y_eta_f = (ye(1,:)*D')';
    end
    len = sqrt(x_eta_f.^2 + y_eta_f.^2) + 1e-300;
    nx = sign(side) * y_eta_f./len;
    ny = sign(side) * -x_eta_f./len;
    ds = len .* w(:);
end

function [nx,ny,ds] = face_normal_eta(xe, ye, D, w, side)
    if side > 0
        x_xi_f = D*xe(:,end);
        y_xi_f = D*ye(:,end);
    else
        x_xi_f = D*xe(:,1);
        y_xi_f = D*ye(:,1);
    end
    len = sqrt(x_xi_f.^2 + y_xi_f.^2) + 1e-300;
    nx = sign(side) * -y_xi_f./len;
    ny = sign(side) * x_xi_f./len;
    ds = len .* w(:);
end