# DGSEM Solver — NACA 0012 Viscous Flow (MATLAB)

MATLAB implementation of the Discontinuous Galerkin Spectral Element Method (DGSEM) for 2-D compressible viscous flow around the NACA 0012 airfoil.

## Method

- **Equations:** 2-D compressible Navier–Stokes (conservative form)
- **Spatial discretisation:** DGSEM, polynomial degree N = 3, Legendre–Gauss–Lobatto quadrature
- **Inviscid flux:** Roe approximate Riemann solver with Harten–Hyman entropy fix
- **Viscous flux:** BR2 scheme
- **Time integration:** Classical RK4 (explicit)
- **Mesh:** C-grid with cosine clustering on the airfoil surface and geometric radial stretching to the farfield
- **Boundary conditions:** Characteristic (Riemann-invariant) far-field BC; no-slip adiabatic wall BC via ghost cells

## Files

| File | Description |
|---|---|
| `main_NACA0012.m` | Main simulation driver |
| `generate_mesh_NACA0012.m` | C-grid mesh generation with cosine-clustered airfoil nodes, radial stretching, and curvilinear metric computation |
| `LGL_quadrature.m` | LGL nodes and weights |
| `derivative_matrix.m` | Spectral derivative matrix (SBP property verified) |
| `rhs_DGSEM.m` | Spatial residual — contravariant fluxes with curvilinear metric terms, Lax-Friedrichs face fluxes |
| `rk4_step.m` | Classical four-stage Runge–Kutta time integrator |
| `roe_flux.m` | Roe approximate Riemann solver with entropy fix |
| `initialize_solution.m` | Free-stream initial condition (non-dimensional, with angle of attack) |
| `compute_timestep.m` | CFL-based time step (inviscid + viscous Fourier limit) |
| `compute_aero_coefficients.m` | CL, CD, CM via surface pressure integration with LGL quadrature |
| `compute_pressure_coefficient.m` | Surface Cp distribution |
| `compute_total_energy.m` | Jacobian-weighted energy conservation check |
| `wall_bc.m` | No-slip adiabatic wall boundary condition (ghost cell) |
| `farfield_bc.m` | Characteristic far-field boundary condition (Riemann invariants) |
| `diagnostic_check.m` | Initialization and timestep diagnostics |

### Test / Verification Files

| File | Description |
|---|---|
| `test_LGL_and_D.m` | Validates LGL nodes/weights, SBP property, polynomial exactness |
| `test_roe_flux.m` | Validates Roe solver: consistency, conservation, entropy fix, rotational invariance |
| `test_boundary_conditions.m` | Validates wall BC (no-slip, adiabatic) and far-field BC (characteristic) |
| `test_initialize_solution.m` | Validates Q layout, AoA decomposition, energy consistency |
| `test_mesh.m` | Validates Jacobian positivity, GCL (free-stream preservation), airfoil geometry |
| `test_isentropic_vortex.m` | Verification test with exact solution (isentropic vortex) |

## Quick Start

```matlab
% 1. Run validation tests (recommended before simulation)
test_LGL_and_D
test_roe_flux
test_boundary_conditions
test_initialize_solution
test_mesh                    % also generates mesh_validation.png

% 2. NACA 0012 simulation
main_NACA0012

% 3. Diagnostics
diagnostic_check
```

## Key Parameters (`main_NACA0012.m`)

```matlab
Re    = 1e4;       % Reynolds number
Mach  = 0.3;       % Mach number
alpha = 5.0;       % Angle of attack (degrees)
N     = 3;         % Polynomial degree
CFL   = 0.10;      % Inviscid CFL number (stability limit ≈ 1/(2N+1) ≈ 0.143)
CFL_v = 0.10;      % Viscous CFL (Fourier number)
nelem_airfoil = 48;  % Elements wrapping the airfoil
nelem_wake    = 16;  % Elements in the wake extension
nelem_radial  = 20;  % Elements in the wall-normal direction
R_farfield    = 30.0; % Far-field radius (chord lengths)
```

## Non-Dimensionalisation

All quantities are non-dimensional:
- `rho_inf = 1`, `c_inf = 1`, `chord = 1`
- `p_inf = 1/gamma`
- `V_inf = Mach`
- `mu = 1/Re`

## Solution Layout

The conservative variable array `Q` is stored as:

```
Q(e, i, j, k)    % e = element, i = xi-node, j = eta-node, k = variable
                  % k: 1=rho, 2=rho*u, 3=rho*v, 4=E
```

Shape: `(n_elem, Np1, Np1, 4)` where `Np1 = N + 1`.

## SBP Property

The derivative matrix **D** and diagonal mass matrix **M** (with LGL weights) satisfy the summation-by-parts identity:

```
M·D + (M·D)ᵀ = B
```

where **B** is the boundary operator. This is verified via `test_LGL_and_D.m`.

## Verification

The isentropic vortex (`test_isentropic_vortex.m`) provides an exact solution to the Euler equations. Expected L2 errors for N = 3 on a 10×10 grid: O(10⁻⁴) – O(10⁻⁶), with 4th-order spatial convergence.

## NACA 0012 Geometry

```
yₜ = (t/0.2) [0.2969√x − 0.1260x − 0.3516x² + 0.2843x³ − 0.1036x⁴]
```

t = 0.12. Trailing edge coefficient −0.1036 for proper closure. The airfoil surface is parametrised with **cosine clustering** to concentrate nodes near the leading and trailing edges.

## Mesh Topology

The C-grid wraps around the airfoil and extends downstream:
- **ξ-direction (circumferential):** lower wake → lower surface → LE → upper surface → upper wake (periodic)
- **η-direction (radial):** airfoil wall (η = −1) → far-field (η = +1)
- Wall boundary at `ej = 1`, far-field boundary at `ej = neta`

## References

1. Kopriva, D. A., *Implementing Spectral Methods for Partial Differential Equations*, Springer, 2009.
2. Hesthaven, J. S. & Warburton, T., *Nodal Discontinuous Galerkin Methods*, Springer, 2008.
3. Bassi, F. & Rebay, S., J. Comput. Phys., 1997.
4. Roe, P. L., J. Comput. Phys., 1981.
5. Toro, E. F., *Riemann Solvers and Numerical Methods for Fluid Dynamics*, 3rd ed., Springer, 2009.
6. Ferrer, E. et al., *DGSEM iLES of NACA0012*, NNFM 143, Springer, 2021.
7. Gregory, N. & O'Reilly, C. L., R&M No. 3726, ARC, 1970.
