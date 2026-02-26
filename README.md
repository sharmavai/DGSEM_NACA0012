# DGSEM Solver — NACA 0012 Viscous Flow (MATLAB)

MATLAB implementation of the Discontinuous Galerkin Spectral Element Method (DGSEM) for 2-D compressible viscous flow around the NACA 0012 airfoil.

## Method

- **Equations:** 2-D compressible Navier–Stokes (conservative form)
- **Spatial discretisation:** DGSEM, polynomial degree N = 3, Legendre–Gauss–Lobatto quadrature
- **Inviscid flux:** Roe approximate Riemann solver with entropy fix
- **Viscous flux:** BR2 scheme
- **Time integration:** Classical RK4 (explicit)
- **Mesh:** C-grid with cosine clustering on airfoil and radial stretching to farfield

## Files

| File | Description |
|---|---|
| `main_NACA0012.m` | Main simulation driver |
| `test_isentropic_vortex.m` | Verification test (exact solution) |
| `LGL_quadrature.m` | LGL points and weights |
| `derivative_matrix.m` | Spectral derivative matrix |
| `rhs_DGSEM.m` | Spatial residual computation |
| `roe_flux.m` | Roe Riemann solver |
| `generate_mesh_NACA0012.m` | C-grid mesh generation |
| `initialize_solution.m` | Free-stream initial condition |
| `compute_timestep.m` | CFL-based time step |
| `compute_aero_coefficients.m` | CL, CD, CM integration |
| `compute_pressure_coefficient.m` | Surface Cp |
| `compute_total_energy.m` | Energy conservation check |
| `diagnostic_check.m` | Initialization / timestep diagnostics |

## Quick Start

```matlab
% Verification (isentropic vortex, exact solution)
test_isentropic_vortex

% NACA 0012 simulation
main_NACA0012

% Diagnostics
diagnostic_check
```

## Key Parameters (`main_NACA0012.m`)

```matlab
Re    = 1e4;      % Reynolds number
Mach  = 0.3;      % Mach number
alpha = 0.0;      % Angle of attack (degrees)
N     = 3;        % Polynomial degree
CFL   = 0.2;      % CFL number (stable limit ≈ 1/(2N+1) ≈ 0.14)
nelem_airfoil = 24;
nelem_wake    = 12;
nelem_radial  = 10;
```

## SBP Property

The derivative matrix **D** and diagonal mass matrix **M** (with LGL weights) satisfy the summation-by-parts identity:

```
M·D + (M·D)ᵀ = B
```

where **B** is the boundary operator. This is verified via `derivative_matrix.m` and `LGL_quadrature.m`.

## Verification

The isentropic vortex (`test_isentropic_vortex.m`) provides an exact solution to the Euler equations. Expected L2 errors for N = 3 on a 10×10 grid: O(10⁻⁴) – O(10⁻⁶), with 4th-order spatial convergence.

## NACA 0012 Geometry

```
yₜ = (t/0.2) [0.2969√x − 0.1260x − 0.3516x² + 0.2843x³ − 0.1036x⁴]
```
t = 0.12. Trailing edge coefficient −0.1036 for proper closure.

## Known Limitations

- Surface normal vectors and face integration use simplified placeholders (full curvilinear metric terms not yet implemented)
- Viscous gradient computation at boundaries is approximate
- Interface flux uses element-local dissipation rather than full face-coupled Roe flux

## References

1. Kopriva, D. A., *Implementing Spectral Methods for Partial Differential Equations*, Springer, 2009.
2. Hesthaven, J. S. & Warburton, T., *Nodal Discontinuous Galerkin Methods*, Springer, 2008.
3. Bassi, F. & Rebay, S., J. Comput. Phys., 1997.
4. Roe, P. L., J. Comput. Phys., 1981.
5. Gregory, N. & O'Reilly, C. L., R&M No. 3726, ARC, 1970.
"# DGSEM_NACA0012" 
