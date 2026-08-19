# Cylinder CANN Material Case

Small cylinder regression case for the legacy CANN `<Add_row>` material input
path.

The committed `solver.xml` is intentionally minimal and self-contained:

- the mesh and surface files are tracked in `mesh/`;
- the pressure history is tracked as `inner_pressure.dat`;
- four tracked fiber VTU files are read from `mesh/`;
- generated result folders, plot images, and scratch solver variants are not
  part of the case.

## Material Rows

Each CANN term is one `<Add_row>` block:

```xml
<Add_row row_name="iso_I1">
  <Invariant_num> 1 </Invariant_num>
  <Activation_functions> (1, 1, 1) </Activation_functions>
  <Weights> (1.0, 1.0, 5.0e3) </Weights>
</Add_row>
```

This case uses one isotropic `I1` row and four recruited fiber-family rows.
The mesh block reads four `FIB_DIR` cell-data files, following the same
`Fiber_direction_file_path` pattern used by the LV CANN cases:

```xml
<Fiber_direction_file_path> mesh/fibersD1Cells.vtu </Fiber_direction_file_path>
<Fiber_direction_file_path> mesh/fibersD2Cells.vtu </Fiber_direction_file_path>
<Fiber_direction_file_path> mesh/fibersZCells.vtu </Fiber_direction_file_path>
<Fiber_direction_file_path> mesh/fibersThetaCells.vtu </Fiber_direction_file_path>
```

The legacy CANN path maps those four directions as `D1`, `D2`, `z`, and
`theta`:

```text
Invariant_num 4   D1
Invariant_num 8   D2
Invariant_num 11  z
Invariant_num 10  theta
```

The fiber files store local cylindrical cell directions. `theta` is the
circumferential direction, `z` is axial, and `D1`/`D2` are a symmetric
opposite-handed helical pair:

```text
D1 = cos(28.807 deg) * theta + sin(28.807 deg) * z
D2 = cos(28.807 deg) * theta - sin(28.807 deg) * z
```

The checked-in rows include beta recruitment:

```xml
<Add_row row_name="recruited_D1">
  <Invariant_num> 4 </Invariant_num>
  <Activation_functions> (1, 2, 2) </Activation_functions>
  <Weights> (1.0, 1.0, 2.5e3) </Weights>
  <Dispersion> 0.0231 </Dispersion>
  <Recruitment distribution="beta">
    <Lower_stretch> 1.00 </Lower_stretch>
    <Upper_stretch> 1.30 </Upper_stretch>
    <Tau> 0.02 </Tau>
    <Alpha> 2.0 </Alpha>
    <Beta> 5.0 </Beta>
    <Quadrature_points> 24 </Quadrature_points>
  </Recruitment>
</Add_row>
<Add_row row_name="recruited_D2">
  <Invariant_num> 8 </Invariant_num>
  <Activation_functions> (1, 2, 2) </Activation_functions>
  <Weights> (1.0, 1.0, 2.5e3) </Weights>
  <Dispersion> 0.0382 </Dispersion>
  <Recruitment distribution="beta">
    <Lower_stretch> 1.00 </Lower_stretch>
    <Upper_stretch> 1.30 </Upper_stretch>
    <Tau> 0.02 </Tau>
    <Alpha> 2.0 </Alpha>
    <Beta> 5.0 </Beta>
    <Quadrature_points> 24 </Quadrature_points>
  </Recruitment>
</Add_row>
```

`Dispersion` is optional and defaults to `0.0`. For supported fiber rows it
applies the GOH-style dispersed feature:

```text
kappa * (I1 - 3) + (1 - 3*kappa) * (I4 - 1)
```

`Recruitment` is optional on `Invariant_num 4`, `8`, `10`, and `11`. The
current implementation supports `distribution="beta"` with explicit lower and
upper stretch bounds, smoothing `Tau`, beta shape parameters, and quadrature
point count.

## Supported Legacy Invariants

```text
1  I1
2  I2
3  J
4  I4_D1
5  I5_D1
6  I8_D1D2
7  I9_D1D2
8  I4_D2
9  I5_D2
10 Etheta
11 Ez
12 E11
13 E22
14 E33
15 E12
16 E23
17 E13
18 EthetaEz
```

Activation functions are stored as `(h0, h1, h2)`:

```text
h0 = 1  identity
h0 = 2  positive Macaulay bracket
h0 = 3  absolute value

h1 = 1  linear
h1 = 2  square
h1 = 3  cubic

h2 = 1  linear
h2 = 2  exponential, exp(w*x)-1
h2 = 3  logarithmic, -log(1-w*x)
```

`solver_NH.xml` remains as a hardcoded Neo-Hookean comparison case.

## HGO Comparison

`solver_HGO.xml` is a hardcoded `HolzapfelOgden` comparison case using the
same mesh, pressure load, and `D1`/`D2` fiber files. It omits recruitment and
matches the non-recruited CANN rows as:

```text
CANN iso_I1             -> HGO a = 5.0e3, b = 0.0
CANN recruited_D1 base  -> HGO a4f = 2.5e3, b4f = 1.0
CANN recruited_D2 base  -> HGO a4s = 2.5e3, b4s = 1.0
```

The hardcoded HGO material does not include the extra recruited `z` and `theta`
rows from `solver.xml`; it is intended as a no-recruitment comparison for the
two HGO-like fiber rows.

`solver_HGO_dispersion.xml` uses the hardcoded `HGO` material with
`kappa = 0.0231`. `solver_CANN_HGO_dispersion.xml` repeats the CANN-HGO
comparison with the same GOH-style dispersion value enabled on both fiber rows.
