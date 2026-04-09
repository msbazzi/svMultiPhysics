This case inflates an idealized left ventricle mesh using the **Holzapfel–Gasser–Ogden (HGO)** anisotropic hyperelastic model with two fiber families (longitudinal and sheet directions) read from VTK cell data. Units are cgs; geometry is scaled with `Mesh_scale_factor` like in the other LV struct cases.

Parameters `a4`, `b4`, `a6`, `b6`, and `kappa` map to the solver’s HGO implementation (dispersed fiber structure via `kappa`). Endocardial pressure follows `endo_pressure.dat`.

**Formula cross-check:** `solver_formula.xml` uses `Constitutive_model type="Formula"` with an isochoric strain energy equal to the built-in HGO model:

- `C10 = E/(4(1+nu))` (same as `read_mat_model`: shear modulus `mu = E/(2(1+nu))`, then `C10 = mu/2`).
- `Eff = κ I1 + (1−3κ) I4 − 1`, `Ess = κ I1 + (1−3κ) I8 − 1`, with `I1`, `I4`, `I8` the **modified** invariants used by the CANN helper (same as in `ArtificialNeuralNetMaterial::computeInvariantsAndDerivatives`).
- `Ψ_iso = C10(I1−3) + (a4/(2b4))(exp(b4 Eff²)−1) + (a6/(2b6))(exp(b6 Ess²)−1)`.

Volumetric response is still ST91 with the same `Penalty_parameter` as `solver.xml`. Running both inputs and comparing `Displacement` at step 5 should agree to near machine precision for this mesh.

Example:

```bash
mpirun -np 1 /path/to/svmultiphysics solver.xml
mpirun -np 1 /path/to/svmultiphysics solver_formula.xml
```
