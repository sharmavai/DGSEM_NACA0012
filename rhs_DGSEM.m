function dUdt = rhs_DGSEM(U, mesh, N, D, w, gamma, mu, Pr, R, U_inf, alpha, p_inf, T_inf)
% Non-dimensional: rho_inf=1, U_inf=1, chord=1, p_inf=1/gamma

nelem = mesh.nelem;
dUdt  = zeros(size(U));
[n_radial, n_circ, ~] = size(mesh.x);

alpha_rad = alpha * pi / 180;
u_inf  = U_inf * cos(alpha_rad);
v_inf  = U_inf * sin(alpha_rad);
rho_ff = 1.0;
p_ff   = p_inf;
E_ff   = p_ff/((gamma-1)*rho_ff) + 0.5*(u_inf^2+v_inf^2);
Uff    = [rho_ff; rho_ff*u_inf; rho_ff*v_inf; rho_ff*E_ff];

is_wall = false(nelem,1);  is_wall(mesh.airfoil_ids)  = true;
is_far  = false(nelem,1);  is_far(mesh.farfield_ids)   = true;

% neighbour index map
nbr = zeros(nelem,4);  % 1=i+, 2=i-, 3=j+, 4=j-
for j = 1:n_circ
    for i = 1:n_radial
        e = i + (j-1)*n_radial;
        if i < n_radial, nbr(e,1) = (i+1)+(j-1)*n_radial; end
        if i > 1,        nbr(e,2) = (i-1)+(j-1)*n_radial; end
        jp = mod(j, n_circ)+1;  jm = mod(j-2,n_circ)+1;
        nbr(e,3) = i+(jp-1)*n_radial;
        nbr(e,4) = i+(jm-1)*n_radial;
    end
end

for e = 1:nelem
    Ue = squeeze(U(e,:,:,:));
    [rho,u,v,p,H,ok] = get_prim(Ue, gamma);
    if ~ok, dUdt(e,:,:,:) = 0; continue; end

    F = Fx(rho,u,v,p,H);
    G = Gy(rho,u,v,p,H);
    dFdx = zeros(N+1,N+1,4);
    dGdy = zeros(N+1,N+1,4);
    for k=1:4
        dFdx(:,:,k) = D*F(:,:,k);
        dGdy(:,:,k) = (D*G(:,:,k)')';
    end
    Res = dFdx + dGdy;

    c_el  = sqrt(gamma*p./rho);
    lam_e = max(max(sqrt(u.^2+v.^2) + c_el));

    % --- face 1: xi=+1 (radial outward, i=N+1 row) ---
    nb = nbr(e,1);
    if nb > 0 && ~is_far(e)
        Une = squeeze(U(nb,:,:,:));
        [rn,un,vn,pn,Hn,okn] = get_prim(Une,gamma);
        if okn
            lf = max(lam_e, max(max(sqrt(un.^2+vn.^2)+sqrt(gamma*pn./rn))));
            Fn = Fx(rn,un,vn,pn,Hn);
            for k=1:4
                num = 0.5*(F(N+1,:,k)+Fn(1,:,k)) - 0.5*lf*(Une(1,:,k)-Ue(N+1,:,k));
                Res(N+1,:,k) = Res(N+1,:,k) + (num-F(N+1,:,k))/w(N+1);
            end
        end
    elseif is_far(e)
        lf = abs(u_inf)+sqrt(gamma*p_ff/rho_ff);
        for j=1:N+1
            fin = Fxv(squeeze(Ue(N+1,j,:)), gamma);
            ff  = Fxv(Uff, gamma);
            num = 0.5*(fin+ff) - 0.5*lf*(Uff - squeeze(Ue(N+1,j,:)));
            for k=1:4, Res(N+1,j,k) = Res(N+1,j,k)+(num(k)-fin(k))/w(N+1); end
        end
    end

    % --- face 2: xi=-1 (radial inward, i=1 row) ---
    nb = nbr(e,2);
    if nb > 0
        Une = squeeze(U(nb,:,:,:));
        [rn,un,vn,pn,Hn,okn] = get_prim(Une,gamma);
        if okn
            lf = max(lam_e, max(max(sqrt(un.^2+vn.^2)+sqrt(gamma*pn./rn))));
            Fn = Fx(rn,un,vn,pn,Hn);
            for k=1:4
                num = 0.5*(Fn(N+1,:,k)+F(1,:,k)) - 0.5*lf*(Ue(1,:,k)-Une(N+1,:,k));
                Res(1,:,k) = Res(1,:,k) - (num-F(1,:,k))/w(1);
            end
        end
    elseif is_wall(e)
        for j=1:N+1
            Ui = squeeze(Ue(1,j,:));
            ri = max(Ui(1),1e-12);
            Ug = [ri; -Ui(2); -Ui(3); Ui(4)];  % mirror velocity → no-slip
            fin = Fxv(Ui, gamma);
            num = LF(Ui, Ug, gamma, lam_e);
            for k=1:4, Res(1,j,k) = Res(1,j,k)-(num(k)-fin(k))/w(1); end
        end
    end

    % --- face 3: eta=+1 (j=N+1 col) ---
    nb = nbr(e,3);
    if nb > 0
        Une = squeeze(U(nb,:,:,:));
        [rn,un,vn,pn,Hn,okn] = get_prim(Une,gamma);
        if okn
            lf = max(lam_e, max(max(sqrt(un.^2+vn.^2)+sqrt(gamma*pn./rn))));
            Gn = Gy(rn,un,vn,pn,Hn);
            for k=1:4
                num = 0.5*(G(:,N+1,k)+Gn(:,1,k)) - 0.5*lf*(Une(:,1,k)-Ue(:,N+1,k));
                Res(:,N+1,k) = Res(:,N+1,k) + (num-G(:,N+1,k))/w(N+1);
            end
        end
    end

    % --- face 4: eta=-1 (j=1 col) ---
    nb = nbr(e,4);
    if nb > 0
        Une = squeeze(U(nb,:,:,:));
        [rn,un,vn,pn,Hn,okn] = get_prim(Une,gamma);
        if okn
            lf = max(lam_e, max(max(sqrt(un.^2+vn.^2)+sqrt(gamma*pn./rn))));
            Gn = Gy(rn,un,vn,pn,Hn);
            for k=1:4
                num = 0.5*(Gn(:,N+1,k)+G(:,1,k)) - 0.5*lf*(Ue(:,1,k)-Une(:,N+1,k));
                Res(:,1,k) = Res(:,1,k) - (num-G(:,1,k))/w(1);
            end
        end
    end

    for jj=1:N+1
        for ii=1:N+1
            dUdt(e,ii,jj,:) = -squeeze(Res(ii,jj,:)) / (w(ii)*w(jj));
        end
    end
end
end

function [rho,u,v,p,H,ok] = get_prim(Ue, gamma)
    rho = max(Ue(:,:,1), 1e-12);
    u   = Ue(:,:,2)./rho;
    v   = Ue(:,:,3)./rho;
    E   = Ue(:,:,4)./rho;
    p   = (gamma-1)*rho.*(E - 0.5*(u.^2+v.^2));
    p   = max(p, 1e-12);
    H   = E + p./rho;
    ok  = all(isfinite(rho(:))) && all(isfinite(p(:)));
end

function F = Fx(rho,u,v,p,H)
    F = cat(3, rho.*u, rho.*u.*u+p, rho.*u.*v, rho.*u.*H);
end

function G = Gy(rho,u,v,p,H)
    G = cat(3, rho.*v, rho.*u.*v, rho.*v.*v+p, rho.*v.*H);
end

function f = Fxv(U4, gamma)
    U4 = U4(:);
    r = max(U4(1),1e-12); u=U4(2)/r; v=U4(3)/r; E=U4(4)/r;
    p = max((gamma-1)*r*(E-0.5*(u^2+v^2)), 1e-12);
    H = E+p/r;
    f = [r*u; r*u*u+p; r*u*v; r*u*H];
end

function fn = LF(UL, UR, gamma, lam)
    fL = Fxv(UL,gamma); fR = Fxv(UR,gamma);
    fn = 0.5*(fL+fR) - 0.5*lam*(UR(:)-UL(:));
end