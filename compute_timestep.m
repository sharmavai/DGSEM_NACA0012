function dt = compute_timestep(U, mesh, N, CFL, gamma)

[n_radial, n_circ, ~] = size(mesh.x);
dt = inf;

for j = 1:n_circ
    for i = 1:n_radial
        e = i + (j-1)*n_radial;

        lam_elem = 0;
        for jj = 1:N+1
            for ii = 1:N+1
                rho  = U(e, ii, jj, 1);
                rhou = U(e, ii, jj, 2);
                rhov = U(e, ii, jj, 3);
                rhoE = U(e, ii, jj, 4);

                rho  = max(rho,  1e-10);
                u    = rhou / rho;
                v    = rhov / rho;
                E    = rhoE / rho;
                p    = max((gamma-1)*rho*(E - 0.5*(u^2+v^2)), 1e-10);
                c    = sqrt(gamma * p / rho);
                lam_elem = max(lam_elem, sqrt(u^2 + v^2) + c);
            end
        end

        lam_elem = max(lam_elem, 1e-10);

        x1=mesh.x(i,j,1); y1=mesh.y(i,j,1);
        x2=mesh.x(i,j,2); y2=mesh.y(i,j,2);
        x3=mesh.x(i,j,3); y3=mesh.y(i,j,3);
        x4=mesh.x(i,j,4); y4=mesh.y(i,j,4);

        h = min([ sqrt((x2-x1)^2+(y2-y1)^2), sqrt((x3-x2)^2+(y3-y2)^2), ...
                  sqrt((x4-x3)^2+(y4-y3)^2), sqrt((x1-x4)^2+(y1-y4)^2) ]);
        h = max(h, 1e-12);

        dt = min(dt, CFL * h / (lam_elem * (2*N+1)));
    end
end

if ~isfinite(dt) || dt <= 0
    error('Invalid time step: dt = %e', dt);
end
end