function [CL, CD, CM] = compute_aero_coefficients(U, mesh, N, gamma, rho_inf, U_inf, chord)
% Surface pressure integration on airfoil wall elements.
%
% Airfoil path winds CLOCKWISE (TE -> lower -> LE -> upper -> TE).
% Outward normal into fluid = rotate tangent (dx,dy) by +90 deg CCW:
%   nx = -dy/seg,  ny = dx/seg
%
% Force ON airfoil:  dF = -p * n_outward * ds
%   Fx += -p * nx * ds
%   Fy += -p * ny * ds
%
% CL = Fy / (q_inf * chord)      [lift perp to freestream, alpha=0 => y]
% CD = -Fx / (q_inf * chord)     [drag opposing freestream direction]
% CM = moment / (q_inf * chord^2)

q_inf = 0.5 * rho_inf * U_inf^2;
[n_radial, ~, ~] = size(mesh.x);

Fx = 0;  Fy = 0;  Mt = 0;

for eid = 1:length(mesh.airfoil_ids)
    e = mesh.airfoil_ids(eid);
    j = ceil(e / n_radial);

    % inner face edge of wall element: mesh corners 1 and 2 at radial layer i=1
    x1 = mesh.x(1, j, 1);  y1 = mesh.y(1, j, 1);
    x2 = mesh.x(1, j, 2);  y2 = mesh.y(1, j, 2);

    dx  = x2 - x1;
    dy  = y2 - y1;
    seg = sqrt(dx^2 + dy^2);
    if seg < 1e-14, continue; end

    % outward normal FROM airfoil INTO fluid
    % (clockwise path => CCW rotation of tangent = outward)
    nx = -dy / seg;
    ny =  dx / seg;

    ds = seg / (N+1);   % arc length per LGL point

    for jj = 1:N+1
        rho  = U(e, 1, jj, 1);
        if ~isfinite(rho) || rho <= 0, continue; end

        u_v  = U(e, 1, jj, 2) / rho;
        v_v  = U(e, 1, jj, 3) / rho;
        E_v  = U(e, 1, jj, 4) / rho;
        p    = (gamma-1) * rho * (E_v - 0.5*(u_v^2 + v_v^2));
        if ~isfinite(p) || p < 0, continue; end

        % position along face
        t   = (jj - 0.5) / (N+1);
        x_s = x1 + t*dx;
        y_s = y1 + t*dy;

        % force on airfoil = -p * n_outward * ds
        Fx = Fx - p * nx * ds;
        Fy = Fy - p * ny * ds;

        % moment about quarter chord  (x=0.25c, y=0)
        Mt = Mt - p * ((x_s - 0.25*chord)*ny - y_s*nx) * ds;
    end
end

q_chord = q_inf * chord;
CL =  Fy / q_chord;
CD = -Fx / q_chord;
CM =  Mt / (q_chord * chord);

if ~isfinite(CL), CL = 0; end
if ~isfinite(CD), CD = 0; end
if ~isfinite(CM), CM = 0; end
end