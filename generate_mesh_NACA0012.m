function mesh = generate_mesh_NACA0012(nelem_airfoil, nelem_wake, nelem_radial, chord)
% GENERATE_MESH_NACA0012  C-grid mesh around NACA 0012 airfoil.

farfield_radius  = 20 * chord;
min_element_size = 0.02 * chord;
t = 0.12;

n_points  = nelem_airfoil + 1;
th_lower  = linspace(0, pi, ceil(n_points/2)+1);
th_upper  = linspace(0, pi, floor(n_points/2)+1);
x_lo = 0.5*(1 - cos(th_lower));
x_up = 0.5*(1 + cos(th_upper));

naca = @(x) (t/0.2)*(0.2969*sqrt(x) - 0.1260*x - 0.3516*x.^2 + 0.2843*x.^3 - 0.1036*x.^4);

x_airfoil = [x_lo(end:-1:1), x_up(2:end)];
y_airfoil = [-naca(x_lo(end:-1:1)), naca(x_up(2:end))];
x_airfoil(1) = chord; y_airfoil(1) = 0;
x_airfoil(end) = chord; y_airfoil(end) = 0;

n_wake = nelem_wake + 1;
x_wake = linspace(chord, chord + 10*chord, n_wake);
y_wake = zeros(size(x_wake));

x_path = [x_airfoil(1:end-1), x_wake(2:end)];
y_path = [y_airfoil(1:end-1), y_wake(2:end)];

% Enforce minimum spacing
s = [0, cumsum(sqrt(diff(x_path).^2 + diff(y_path).^2))];
total_len = s(end);
small = sum(diff(s) < min_element_size*0.5);
nelem_circ = length(x_path) - 1;

if small > nelem_circ*0.1
    n_new = max(nelem_circ, floor(total_len/min_element_size));
    s_new = linspace(0, total_len, n_new+1);
    x_path = interp1(s, x_path, s_new, 'linear');
    y_path = interp1(s, y_path, s_new, 'linear');
    nelem_circ = length(x_path) - 1;
end

mesh.nelem        = nelem_circ * nelem_radial;
mesh.x            = zeros(nelem_radial, nelem_circ, 4);
mesh.y            = zeros(nelem_radial, nelem_circ, 4);
mesh.airfoil_ids  = [];
mesh.farfield_ids = [];
mesh.type         = zeros(mesh.nelem, 1);

eid = 0;
for j = 1:nelem_circ
    xi1 = x_path(j);   yi1 = y_path(j);
    xi2 = x_path(j+1); yi2 = y_path(j+1);

    dx = xi2 - xi1; dy = yi2 - yi1;
    nm = sqrt(dx^2 + dy^2);
    if nm > 1e-10
        nx = -dy/nm; ny = dx/nm;
    else
        nx = 0; ny = 1;
    end

    for i = 1:nelem_radial
        eid = eid + 1;
        r1 = ((i-1)/nelem_radial)^1.5 * farfield_radius;
        r2 = (i    /nelem_radial)^1.5 * farfield_radius;

        mesh.x(i,j,1) = xi1 + r1*nx;  mesh.y(i,j,1) = yi1 + r1*ny;
        mesh.x(i,j,2) = xi2 + r1*nx;  mesh.y(i,j,2) = yi2 + r1*ny;
        mesh.x(i,j,3) = xi2 + r2*nx;  mesh.y(i,j,3) = yi2 + r2*ny;
        mesh.x(i,j,4) = xi1 + r2*nx;  mesh.y(i,j,4) = yi1 + r2*ny;

        if i == 1 && j <= nelem_airfoil
            mesh.airfoil_ids(end+1) = eid;
            mesh.type(eid) = 2;
        elseif i == nelem_radial
            mesh.farfield_ids(end+1) = eid;
            mesh.type(eid) = 3;
        else
            mesh.type(eid) = 1;
        end
    end
end

fprintf('Mesh: %d elements (%d airfoil, %d farfield)\n', ...
        mesh.nelem, length(mesh.airfoil_ids), length(mesh.farfield_ids));
end
