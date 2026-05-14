#!/usr/bin/env python3
"""Create external cell data for kinematic-growth mechanics.

The script reads a VTU volume mesh, computes each cell centroid, and writes a
VTU with three CellData arrays:

- Growth_Fg: 9-component row-major growth tensor.
- Material_C10: scalar NeoHookean isochoric parameter.
- Material_Kpen: scalar volumetric penalty parameter.

By default the growth stretch varies linearly along the mesh's longest bounding
box direction. Use --axis to force x, y, or z.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import vtk


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", default="mesh/mesh-complete.mesh.vtu", help="Input VTU volume mesh.")
    parser.add_argument("--output", default="mesh/growth_properties.vtu", help="Output VTU with CellData arrays.")
    parser.add_argument(
        "--axis",
        choices=("auto", "x", "y", "z"),
        default="auto",
        help="Axial direction used for the growth gradient.",
    )
    parser.add_argument("--theta-start", type=float, default=1.00, help="Growth stretch at the minimum axial coordinate.")
    parser.add_argument("--theta-end", type=float, default=1.08, help="Growth stretch at the maximum axial coordinate.")
    parser.add_argument("--c10", type=float, default=5.0e4, help="Constant Material_C10 value.")
    parser.add_argument("--kpen", type=float, default=5.0e5, help="Constant Material_Kpen value.")
    parser.add_argument(
        "--growth-mode",
        choices=("isotropic", "axial"),
        default="isotropic",
        help="Use theta*I or grow only along the axial tensor component.",
    )
    return parser.parse_args()


def read_grid(path: Path) -> vtk.vtkUnstructuredGrid:
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(str(path))
    reader.Update()

    grid = reader.GetOutput()
    if grid is None or grid.GetNumberOfCells() == 0 or grid.GetNumberOfPoints() == 0:
        raise RuntimeError(
            f"Could not read a non-empty VTU mesh from '{path}'. "
            "Check that the file is readable by vtkXMLUnstructuredGridReader."
        )
    return grid


def write_grid(path: Path, grid: vtk.vtkUnstructuredGrid) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetFileName(str(path))
    writer.SetInputData(grid)
    writer.SetDataModeToAscii()
    writer.Write()


def choose_axis(grid: vtk.vtkUnstructuredGrid, axis_name: str) -> int:
    if axis_name != "auto":
        return {"x": 0, "y": 1, "z": 2}[axis_name]

    bounds = grid.GetBounds()
    extents = [bounds[1] - bounds[0], bounds[3] - bounds[2], bounds[5] - bounds[4]]
    return max(range(3), key=lambda i: extents[i])


def cell_centroid(grid: vtk.vtkUnstructuredGrid, cell_id: int) -> tuple[float, float, float]:
    cell = grid.GetCell(cell_id)
    point_ids = cell.GetPointIds()
    centroid = [0.0, 0.0, 0.0]
    npts = point_ids.GetNumberOfIds()

    for i in range(npts):
        point = grid.GetPoint(point_ids.GetId(i))
        centroid[0] += point[0]
        centroid[1] += point[1]
        centroid[2] += point[2]

    return (centroid[0] / npts, centroid[1] / npts, centroid[2] / npts)


def normalized_coordinate(values: list[float], value: float) -> float:
    lo = min(values)
    hi = max(values)
    if abs(hi - lo) < 1.0e-14:
        return 0.0
    return (value - lo) / (hi - lo)


def add_cell_data(grid: vtk.vtkUnstructuredGrid, axis: int, args: argparse.Namespace) -> None:
    ncell = grid.GetNumberOfCells()
    axial_values = [cell_centroid(grid, cell_id)[axis] for cell_id in range(ncell)]

    fg = vtk.vtkDoubleArray()
    fg.SetName("Growth_Fg")
    fg.SetNumberOfComponents(9)

    c10 = vtk.vtkDoubleArray()
    c10.SetName("Material_C10")

    kpen = vtk.vtkDoubleArray()
    kpen.SetName("Material_Kpen")

    for axial_value in axial_values:
        s = normalized_coordinate(axial_values, axial_value)
        theta = args.theta_start + s * (args.theta_end - args.theta_start)

        if args.growth_mode == "isotropic":
            fg_tuple = (theta, 0.0, 0.0, 0.0, theta, 0.0, 0.0, 0.0, theta)
        else:
            diag = [1.0, 1.0, 1.0]
            diag[axis] = theta
            fg_tuple = (diag[0], 0.0, 0.0, 0.0, diag[1], 0.0, 0.0, 0.0, diag[2])

        fg.InsertNextTuple(fg_tuple)
        c10.InsertNextValue(args.c10)
        kpen.InsertNextValue(args.kpen)

    cell_data = grid.GetCellData()
    for name in ("Growth_Fg", "Material_C10", "Material_Kpen"):
        existing = cell_data.GetArray(name)
        if existing is not None:
            cell_data.RemoveArray(name)

    cell_data.AddArray(fg)
    cell_data.AddArray(c10)
    cell_data.AddArray(kpen)


def main() -> None:
    args = parse_args()
    input_path = Path(args.input)
    output_path = Path(args.output)

    grid = read_grid(input_path)
    axis = choose_axis(grid, args.axis)
    add_cell_data(grid, axis, args)
    write_grid(output_path, grid)

    axis_name = ("x", "y", "z")[axis]
    print(f"Wrote {output_path}")
    print(f"Cells: {grid.GetNumberOfCells()}")
    print(f"Axial direction: {axis_name}")
    print(f"Growth stretch: {args.theta_start} -> {args.theta_end}")


if __name__ == "__main__":
    main()
