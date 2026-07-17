# DGSEM Solver — NACA 0012 Viscous Flow (MATLAB)

## 2-D Compressible Navier-Stokes via Discontinuous Galerkin Spectral Element Method

### Overview

This solver implements a high-order discontinuous Galerkin spectral element method (DGSEM) for the two-dimensional compressible Navier-Stokes equations around the NACA 0012 airfoil. The code is designed for academic research and is validated against published benchmarks.

### Method

| Component | Method | Reference |
|-----------|--------|-----------|
| Equations | 2-D compressible Navier-Stokes (conservative form) | — |
| Spatial discretisation | DGSEM, polynomial degree N, LGL quadrature | Kopriva (2009), Hesthaven & Warburton (2008) |
| Inviscid flux | Lax-Friedrichs (Rusanov) at element interfaces | Hesthaven & Warburton (2008) Ch. 6 |
| Viscous flux | BR2 scheme (Bassi-Rebay 2) | Bassi & Rebay (1997), Cockburn & Shu (1998) |
| Time integration | Classical four-stage Runge-Kutta (explicit) | Hesthaven & Warburton (2008) Ch. 6 |
| Mesh | C-grid with cosine clustering, geometric radial stretching | Blazek (2015) Ch. 9 |
| Wall BC | No-slip adiabatic (ghost-cell approach) | Blazek (2015) Sec. 8.4, Kopriva (2009) Sec. 9.2 |
| Far-field BC | Characteristic (Riemann invariants) | Thompson (1987), Poinsot & Lele (1992) |

### Non-Dimensionalisation

All quantities are non-dimensionalised by reference values at infinity:

| Quantity | Reference | Non-dim value |
|----------|-----------|---------------|
| Density | ρ∞ | 1 |
| Velocity | c∞ (speed of sound) | — |
| Length | chord c | 1 |
| Pressure | ρ∞c∞² | — |

Derived quantities:

- p∞ = 1/γ = 0.71429 (for γ = 1.4)
- V∞ = M∞·c∞ = M∞ (in these units)
- μ = 1/Re
- k = μγ / (Pr(γ − 1))

### Solver Equations

The compressible Navier-Stokes equations in conservative form:

∂ρ/∂t + ∇·(ρ**u**) = 0

∂(ρ**u**)/∂t + ∇·(ρ**u**⊗**u** + pI) = ∇·**τ**

∂E/∂t + ∇·((E+p)**u**) = ∇·(**τ**·**u** + k∇T)

where the viscous stress tensor (Stokes hypothesis):

τᵢⱼ = μ(∂uᵢ/∂xⱼ + ∂uⱼ/∂xᵢ − ⅔δᵢⱼ∇·**u**)

### BR2 Viscous Scheme

The Bassi-Rebay 2 (BR2) scheme for viscous fluxes:

1. **Volume integral**: Compute velocity gradients ∂u/∂x, ∂u/∂y at all LGL nodes using the spectral derivative matrix D and curvilinear metric terms
2. **Lifting**: At each element interface, compute the jump in the gradient trace and distribute it to both elements via the BR2 lifting operator with penalty parameter η_B = 1
3. **Corrected gradient**: Add the lifting correction to the reference-space gradient, then transform to physical space via metric terms
4. **Viscous flux**: Compute τ·**n** at all face nodes using the averaged corrected gradient
5. **Energy**: Compute heat flux k∇T for viscous energy contribution

### Solution Layout

The conservative variable array Q is stored as:

```
Q(e, i, j, k)    % e = element index
                  % i = ξ-node (circumferential, Np1 nodes)
                  % j = η-node (radial, Np1 nodes)
                  % k = variable index:
                  %     1 = ρ      (density)
                  %     2 = ρu     (x-momentum)
                  %     3 = ρv     (y-momentum)
                  %     4 = E      (total energy)
```

Shape: `(n_elem, Np1, Np1, 4)` where `Np1 = N + 1`.

### SBP Property

The derivative matrix D and diagonal mass matrix M (with LGL weights) satisfy the summation-by-parts identity:

**M·D + (M·D)ᵀ = B**

where B is the boundary operator diag(−1, 0, ..., 0, 1). This property ensures the stability and conservation of the DGSEM discretisation.

### Files

| File | Purpose |
|------|---------|
| `main_NACA0012.m` | Main simulation driver with grid convergence study |
| `generate_mesh_NACA0012.m` | C-grid mesh generation with cosine clustering and curvilinear metrics |
| `rhs_DGSEM.m` | Spatial residual — inviscid + viscous (BR2) contravariant fluxes with SAT |
| `rk4_step.m` | Classical four-stage Runge-Kutta time integrator |
| `roe_flux.m` | Roe approximate Riemann solver (available but not used by default) |
| `LGL_quadrature.m` | Legendre-Gauss-Lobatto nodes and weights |
| `derivative_matrix.m` | Spectral derivative matrix with SBP verification |
| `initialize_solution.m` | Freestream initial condition (non-dim, with AoA) |
| `compute_timestep.m` | CFL-based time step (inviscid + viscous Fourier limits) |
| `compute_aero_coefficients.m` | CL, CD, CM via surface pressure integration |
| `compute_pressure_coefficient.m` | Surface Cp distribution |
| `compute_total_energy.m` | Jacobian-weighted energy conservation check |
| `wall_bc.m` | No-slip adiabatic wall BC (ghost cell) |
| `farfield_bc.m` | Characteristic far-field BC (Riemann invariants) |
| `diagnostic_check.m` | Pre-flight diagnostics (mesh, init, timestep, GCL) |

### Verification Tests

| File | Purpose |
|------|---------|
| `test_LGL_and_D.m` | LGL nodes/weights, SBP property, polynomial exactness |
| `test_roe_flux.m` | Roe solver: consistency, conservation, entropy fix, rotational invariance |
| `test_boundary_conditions.m` | Wall BC (no-slip, adiabatic) and far-field BC validation |
| `test_initialize_solution.m` | Q layout, AoA decomposition, energy consistency |
| `test_mesh.m` | Jacobian positivity, GCL, unit normals, airfoil geometry |
| `test_isentropic_vortex.m` | Full DGSEM verification with exact solution (grid convergence) |

### Quick Start

```matlab
% Add to path
addpath('path/to/DGSEM_NACA0012');

% 1. Run verification tests (recommended)
test_LGL_and_D              % Spectral operator verification
test_roe_flux               % Roe solver verification
test_boundary_conditions    % BC verification
test_initialize_solution    % Initial condition verification
test_mesh                   % Mesh + GCL verification
test_isentropic_vortex       % Full DGSEM operator verification

% 2. Pre-flight diagnostics
diagnostic_check

% 3. Main simulation (single mesh level)
%    Set do_convergence_study = false in main_NACA0012.m
main_NACA0012

% 4. Grid convergence study (3 mesh levels)
%    Set do_convergence_study = true in main_NACA0012.m
main_NACA0012
```

### Default Parameters

| Parameter | Symbol | Value | Notes |
|-----------|--------|-------|-------|
| Reynolds number | Re | 1×10⁴ | Laminar vortex shedding regime |
| Mach number | M∞ | 0.3 | Low subsonic (compressibility effects moderate) |
| Angle of attack | α | 5.0° | Lift-generating case |
| Specific heats ratio | γ | 1.4 | Ideal diatomic gas |
| Prandtl number | Pr | 0.72 | Air at standard conditions |
| Polynomial degree | N | 3 | P4 elements (4×4 = 16 DOFs per element per variable) |
| Inviscid CFL | — | 0.08 | Stability limit ≈ 1/(2N+1) ≈ 0.143 |
| Viscous CFL | — | 0.08 | Fourier number |
| Airfoil elements | — | 48 | Circumferential elements around airfoil |
| Wake elements | — | 16 | Downstream wake extension |
| Radial elements | — | 20 | Wall-normal direction |
| Far-field radius | — | 30c | ≥30c recommended |
| Simulation time | — | 30 conv. times | 20 transient + 10 averaging |

### NACA 0012 Geometry

```
yₜ = (t/0.2) [0.2969√x − 0.1260x − 0.3516x² + 0.2843x³ − 0.1036x⁴]
```

where t = 0.12 (12% thickness). The trailing-edge coefficient −0.1036 provides proper closure.

### Mesh Topology

The C-grid wraps around the airfoil and extends downstream:

- **ξ-direction (circumferential)**: lower wake → lower surface → LE → upper surface → upper wake (periodic)
- **η-direction (radial)**: airfoil wall (η = −1) → far-field (η = +1)
- Wall boundary at j = 1 (η = −1), far-field boundary at j = Np1 (η = +1)
- Element ordering: e = (ei − 1) + (ej − 1) × nxi_c + 1

### Grid Convergence Study

The solver includes an automated h-refinement study with three mesh levels:

| Level | n_airfoil | n_wake | n_radial | Total elements |
|-------|-----------|--------|----------|----------------|
| L1 (coarse) | 24 | 8 | 10 | 320 |
| L2 (medium) | 48 | 16 | 20 | 1280 |
| L3 (fine) | 72 | 24 | 30 | 2880 |

Convergence rates for CL and CD are estimated using Richardson extrapolation.

### Validation Targets

Based on Ferrer et al. (2021), DGSEM iLES of NACA 0012:

| α [°] | CL (target) | CD (target) |
|--------|-------------|-------------|
| 0 | ~0.000 | 0.040–0.050 |
| 5 | 0.55–0.65 | 0.045–0.060 |
| 10 | 0.95–1.10 | 0.090–0.130 |

### References

1. Kopriva, D. A. (2009). *Implementing Spectral Methods for Partial Differential Equations*. Springer.
2. Hesthaven, J. S. & Warburton, T. (2008). *Nodal Discontinuous Galerkin Methods*. Springer.
3. Bassi, F. & Rebay, S. (1997). A high-order accurate discontinuous finite element method for the numerical solution of the compressible Navier-Stokes equations. *J. Comput. Phys.*, 131, 267–279.
4. Cockburn, B. & Shu, C.-W. (1998). The local discontinuous Galerkin method for time-dependent convection-diffusion systems. *J. Comput. Phys.*, 141, 199–224.
5. Roe, P. L. (1981). Approximate Riemann solvers, parameter vectors, and difference schemes. *J. Comput. Phys.*, 43, 357–372.
6. Toro, E. F. (2009). *Riemann Solvers and Numerical Methods for Fluid Dynamics*, 3rd ed. Springer.
7. Ferrer, E. et al. (2021). *DGSEM iLES of NACA0012*. In: *High-Order Methods for Incompressible Fluid Flow*, NNFM 143, Springer.
8. Thompson, K. W. (1987). Time dependent boundary conditions for hyperbolic systems. *J. Comput. Phys.*, 68, 1–24.
9. Blazek, J. (2015). *Computational Fluid Dynamics: Principles and Applications*, 3rd ed. Elsevier.
10. Gregory, N. & O'Reilly, C. L. (1970). Low-speed aerodynamic characteristics of NACA 0012 airfoil section. *R&M No. 3726*, ARC.
