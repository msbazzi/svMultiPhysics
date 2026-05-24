This case exercises external cell-based kinematic growth for `struct` mechanics.

The mesh is a short solid cylinder made of tetrahedral elements. The same VTU file
also stores the external cell fields used by the solver:

- `Growth_Fg`: 9-component row-major growth tensor per cell.
- `Material_C10`: per-cell NeoHookean isochoric parameter.
- `Material_Kpen`: per-cell volumetric penalty parameter used to soften the
  incompressibility response.

The bottom face is fixed and no external load is applied. The mild isotropic
growth field creates a simple mechanics-only test of the `F = Fe Fg` split.
The VTK result writes both `Def_grad` and `Growth_Fg` for inspection.
This case also enables `Save_deformed_geometry_to_VTK`, so the saved VTU point
coordinates include the computed displacement.
The 8% growth field is ramped over eight solver steps with
`Growth_ramp_steps` to avoid applying a discontinuous cell-wise growth field in
one nonlinear solve.

Create or refresh the external cell data on an existing mesh with:

```bash
python3 create_external_cell_data.py
```

The script writes `mesh/growth_properties.vtu`, which is the file referenced by
`solver.xml`. By default it varies `Growth_Fg` smoothly along the longest mesh
direction and uses a softened volumetric penalty. Override the direction or
stretch range with, for example:

```bash
python3 create_external_cell_data.py --axis z --theta-start 1.0 --theta-end 1.08
```

Run the small smoke test with one MPI rank:

```bash
mpirun -np 1 ../../../../build/svMultiPhysics-build/bin/svmultiphysics solver.xml
```

Regenerate the original simple mesh and cell data with:

```bash
python3 generate_case.py
```
