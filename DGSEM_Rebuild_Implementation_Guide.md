# DGSEM Code Rebuild: Complete Implementation Guide

## Executive Summary

The original DGSEM/NACA code had **fatal architectural flaws** that made it unsuitable for production CFD. This guide provides a **staged rebuild** with verified components.

---

## Critical Defects Fixed

### 1. Mesh Topology (CRITICAL)
**Original flaw:** Wake cut treated as periodic, connecting last element to first.
**Consequence:** Unphysical flow continuity across wake, wrong vortex shedding.
**Fix:** Closed O-grid baseline; true C-grid requires separate face records.

### 2. Wall Tagging (CRITICAL)
**Original flaw:** All inner boundary faces marked as wall (including wake).
**Consequence:** Wake cut had no-slip BC, completely wrong physics.
**Fix:** Only airfoil elements tagged as wall; wake excluded.

### 3. h_min Computation (CRITICAL)
**Original flaw:** h_min = min(J)*2 (Jacobian is area, not length).
**Consequence:** Wrong CFL, wrong timestep, wrong stability.
**Fix:** Compute h_min from actual edge lengths: sqrt(dx_dxi^2 + dy_dxi^2).

### 4. Derivative Orientation (CRITICAL)
**Original flaw:** du_deta = (u*D')' (wrong transpose).
**Consequence:** Eta derivative computed as xi-oriented operation.
**Fix:** du_deta = u*D' (correct orientation).

### 5. GCL Test Setup (CRITICAL)
**Original flaw:** Initialized alpha=0 while mesh.alpha_deg=5.
**Consequence:** Far-field BC used alpha=5, not a free-stream test.
**Fix:** Set mesh.alpha_deg=0 during GCL test.

### 6. Viscous Energy Flux (MAJOR)
**Original flaw:** Omitted at interfaces.
**Consequence:** Incomplete energy equation.
**Fix:** Add theta_x = tau_xx*u + tau_xy*v + k*dT/dx.

### 7. BR2 Implementation (MAJOR)
**Original flaw:** Hand-written pseudo-lifting without proper variational form.
**Consequence:** Not a valid BR2 scheme.
**Fix:** Implement BR1 first for verification, then BR2 from published weak form.

---

## Staged Rebuild Plan

### Stage 1: Baseline Mesh and State (COMPLETE)
- Closed O-grid mesh generator
- Input validation
- Physical state checks
- Pressure-only aero coefficients (explicitly labeled incomplete)
- Diagnostic without false GCL claims

**Files:**
- `generate_mesh_NACA0012.m`
- `compute_timestep.m`
- `initialize_solution.m`
- `wall_bc.m` (Euler slip only)
- `farfield_bc.m` (prescribed freestream)
- `compute_aero_coefficients.m`
- `compute_pressure_coefficient.m`
- `diagnostic_check.m`

### Stage 2: Minimal Euler DGSEM (COMPLETE)
- Face-loop residual structure
- RK4 time integration
- LGL and D matrix
- Periodic test placeholder

**Files:**
- `rhs_DGSEM_euler.m`
- `rk4_step.m`
- `main_euler_test.m`
- `LGL_quadrature.m`
- `derivative_matrix.m`
- `test_euler_periodic.m`

### Stage 3: Euler Verification (TODO)
- Periodic manufactured solution
- Isentropic vortex convergence
- Conservation tests

### Stage 4: Viscous DG (TODO)
- Implement BR1 or SIPG from published weak form
- Viscous manufactured solution test
- Flat-plate boundary layer

### Stage 5: True C-Grid and Wake (TODO)
- Multiblock mesh with face records
- Proper wake-cut treatment
- Total traction (pressure + viscous)

---

## How to Use This Rebuild

### Week 1-2: Study Foundations
1. Read Kopriva Ch. 8, Hesthaven & Warburton Ch. 6
2. Derive LGL nodes, D matrix, SBP property
3. Implement 1D advection with numerical flux

### Week 3-4: Mesh and Euler
1. Run Stage 1 diagnostic
2. Implement periodic manufactured solution
3. Verify formal order (N+1 for DGSEM)

### Week 5-6: Viscous Terms
1. Implement BR1 or SIPG
2. Viscous MMS test
3. Flat-plate case

### Week 7-8: Airfoil
1. True C-grid mesh
2. Total traction
3. NACA 0012 validation

---

## Verification Checklist

Before claiming any result:
- [ ] LGL nodes exact to machine precision
- [ ] D matrix passes polynomial exactness test
- [ ] Mesh has positive Jacobian at oversampled points
- [ ] Periodic constant state produces zero residual
- [ ] Isentropic vortex shows N+1 convergence
- [ ] Viscous MMS shows N+1 convergence
- [ ] Flat-plate Cf matches theory
- [ ] NACA Cp matches experiment

---

## Professor's Questions (Prepare Answers)

1. **Why is a Jacobian not a cell length?**
   - J measures area scaling for 2D map; h must be a length for CFL.

2. **What is the difference between a wake cut and a periodic interface?**
   - Wake cut: physical discontinuity, separate boundaries.
   - Periodic: mathematically identified boundaries with matching coordinates.

3. **Write the full traction vector for a Newtonian fluid.**
   - t = -p*n + tau*n, where tau = mu*(grad u + grad u^T - 2/3 div u I).

4. **What discrete metric identity ensures free-stream preservation?**
   - Sum of contravariant metric derivatives must vanish: d/dxi(xi_x) + d/deta(eta_x) = 0.

5. **Why is the original GCL test invalid?**
   - Initialized alpha=0 but mesh.alpha_deg=5, so far-field BC used wrong data.

6. **What evidence proves a BR2 implementation is correct?**
   - Viscous MMS with N+1 convergence, flat-plate Cf, and comparison to published BR2 results.

---

## Next Actions

1. **Do NOT run NACA production cases** until Euler and viscous verification pass.
2. **Implement periodic manufactured solution** first.
3. **Add face-record mesh API** for true C-grid.
4. **Document every assumption** in a research notebook.

---

## References

1. Kopriva, D. A. (2009). Implementing Spectral Methods for PDEs. Springer.
2. Hesthaven, J. S. & Warburton, T. (2008). Nodal DG Methods. Springer.
3. Bassi, F. & Rebay, S. (1997). J. Comput. Phys. 131, 267-279.
4. Arnold, D. N. et al. (2002). SIAM J. Numer. Anal. 40, 1-30.
5. NASA Turbulence Modeling Resource: https://www.nasa.gov/nasa-turbulence-modeling-resource/

---

**Generated:** Stage 1 and Stage 2 code packages available in `output/DGSEM_rebuild_stage1/` and `output/DGSEM_rebuild_stage2/`.

**Status:** Baseline repair complete. Euler verification and viscous implementation pending.
