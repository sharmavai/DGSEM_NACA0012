# DGSEM Code Rebuild: Complete Deliverables Index

## 📁 Directory Structure

```
output/
├── DGSEM_rebuild_stage1/          # Baseline mesh and state (9 files)
│   ├── generate_mesh_NACA0012.m
│   ├── compute_timestep.m
│   ├── initialize_solution.m
│   ├── wall_bc.m
│   ├── farfield_bc.m
│   ├── compute_aero_coefficients.m
│   ├── compute_pressure_coefficient.m
│   ├── diagnostic_check.m
│   └── README_STAGE1.md
│
├── DGSEM_rebuild_stage2/          # Minimal Euler DGSEM (8 files)
│   ├── rhs_DGSEM_euler.m
│   ├── rk4_step.m
│   ├── main_euler_test.m
│   ├── LGL_quadrature.m
│   ├── derivative_matrix.m
│   ├── compute_total_energy.m
│   ├── test_euler_periodic.m
│   └── README_STAGE2.md
│
├── DGSEM_Rebuild_Implementation_Guide.md    # Complete rebuild plan
├── DGSEM_Rebuild_Summary.md                 # Executive summary
├── DGSEM_Critical_Fixes_Checklist.txt       # Line-by-line defect list
└── DGSEM_Rebuild_Index.md                   # This file
```

---

## 📦 Package Contents

### Stage 1: Baseline Mesh and State

**Purpose:** Provide a geometrically valid mesh and physical state checks.

**Key Fixes:**
- Closed O-grid (no periodic wake)
- Proper wall tagging (airfoil only)
- h_min from edge lengths
- Input validation
- Physical state checks

**Files:** 9 (1 README + 8 MATLAB)

**Lines of Code:** ~450

**Status:** ✅ COMPLETE

---

### Stage 2: Minimal Euler DGSEM

**Purpose:** Provide correct Euler residual structure.

**Key Features:**
- Face-loop residual
- Correct derivative orientation
- RK4 time integration
- LGL and D matrix
- Conservation check

**Files:** 8 (1 README + 7 MATLAB)

**Lines of Code:** ~350

**Status:** ✅ COMPLETE

---

### Documentation

**Purpose:** Guide the rebuild process and document scope.

**Files:**
1. `DGSEM_Rebuild_Implementation_Guide.md` (5.9 KB)
   - Staged rebuild plan
   - Verification checklist
   - Professor's questions
   - References

2. `DGSEM_Rebuild_Summary.md` (6.0 KB)
   - Executive summary
   - Critical defects fixed
   - Next steps
   - Code statistics

3. `DGSEM_Critical_Fixes_Checklist.txt` (3.0 KB)
   - Line-by-line defect list
   - Status tracking
   - Critical path

**Total Documentation:** ~15 KB

---

## 🎯 Usage Instructions

### For Students

1. **Start with Stage 1**
   - Run `diagnostic_check.m`
   - Verify mesh quality
   - Study mesh structure

2. **Proceed to Stage 2**
   - Run `main_euler_test.m`
   - Implement periodic MMS
   - Verify convergence

3. **Use Documentation**
   - Read Implementation Guide
   - Answer Professor's questions
   - Complete verification checklist

### For Researchers

1. **Baseline First**
   - Do NOT skip Stage 1/2
   - Verify Euler before adding viscosity
   - Document every assumption

2. **Verification Required**
   - Periodic MMS (Stage 3)
   - Isentropic vortex (Stage 3)
   - Viscous MMS (Stage 4)
   - Flat-plate (Stage 4)

3. **Production Cases**
   - Only after all verification passes
   - NACA 0012 validation (Stage 5)
   - Document mesh, BCs, timestep, averaging

---

## 📊 Code Statistics

| Metric | Stage 1 | Stage 2 | Documentation | Total |
|--------|---------|---------|---------------|-------|
| Files | 9 | 8 | 3 | 20 |
| MATLAB files | 8 | 7 | 0 | 15 |
| Documentation files | 1 | 1 | 3 | 5 |
| Lines of code | ~450 | ~350 | ~500 | ~1300 |
| Functions | 12 | 8 | 0 | 20 |

---

## ✅ Verification Status

| Test | Stage 1 | Stage 2 | Stage 3 | Stage 4 | Stage 5 |
|------|---------|---------|---------|---------|---------|
| Mesh quality | ✅ | ✅ | TODO | TODO | TODO |
| State validation | ✅ | ✅ | TODO | TODO | TODO |
| Euler residual | N/A | ✅ | TODO | TODO | TODO |
| Periodic MMS | N/A | N/A | TODO | TODO | TODO |
| Isentropic vortex | N/A | N/A | TODO | TODO | TODO |
| Viscous MMS | N/A | N/A | N/A | TODO | TODO |
| Flat-plate | N/A | N/A | N/A | TODO | TODO |
| NACA validation | N/A | N/A | N/A | N/A | TODO |

**Current Stage:** 2 (Euler baseline complete)

**Next Milestone:** Stage 3 (Euler verification)

---

## 📚 References

### Core Textbooks
1. Kopriva, D. A. (2009). Implementing Spectral Methods for PDEs. Springer.
2. Hesthaven, J. S. & Warburton, T. (2008). Nodal DG Methods. Springer.
3. LeVeque, R. J. (2002). Finite Volume Methods for Hyperbolic Problems. Cambridge.

### Seminal Papers
4. Cockburn, B. & Shu, C.-W. (1989). TVB RKDG for conservation laws. J. Comput. Phys.
5. Bassi, F. & Rebay, S. (1997). BR2 scheme. J. Comput. Phys.
6. Arnold, D. N. et al. (2002). Unified DG analysis. SIAM J. Numer. Anal.

### Resources
7. NASA Turbulence Modeling Resource: https://www.nasa.gov/nasa-turbulence-modeling-resource/

---

## 🎓 Professor's Questions

Prepare answers to these before proceeding:

1. Why is a Jacobian not a cell length?
2. What is the difference between a wake cut and a periodic interface?
3. Write the full traction vector for a Newtonian fluid.
4. What discrete metric identity ensures free-stream preservation?
5. Why is the original GCL test invalid?
6. What evidence proves a BR2 implementation is correct?

---

## ⚠️ Important Warnings

1. **Do NOT run NACA production cases** until Euler and viscous verification pass.
2. **Stage 1 aero coefficients are pressure-only** — do not report CD as total drag.
3. **Stage 1 wall BC is Euler slip** — no-slip requires verified viscous DG.
4. **Stage 1 far-field is prescribed** — characteristic NS BC requires verification.
5. **Document every assumption** in a research notebook.

---

## 📅 Timeline

| Week | Milestone | Deliverable |
|------|-----------|-------------|
| 1-2 | Stage 1/2 study | Run diagnostics, verify mesh |
| 3-4 | Stage 3 | Periodic MMS, isentropic vortex |
| 5-6 | Stage 4 | Viscous MMS, flat-plate |
| 7-8 | Stage 5 | True C-grid, NACA validation |

---

**Generated:** 2026-08-21
**Status:** Stage 1 and Stage 2 complete
**Next:** Stage 3 - Euler verification
**Owner:** DGSEM Rebuild Team
