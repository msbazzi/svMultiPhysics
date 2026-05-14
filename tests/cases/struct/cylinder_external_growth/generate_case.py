import math
from pathlib import Path

import vtk


BASE = Path(__file__).resolve().parent
MESH_DIR = BASE / "mesh"
SURF_DIR = MESH_DIR / "mesh-surfaces"

RADIUS = 1.0
LENGTH = 2.0
N_SEG = 12
Z_LEVELS = [0.0, LENGTH]


def make_points():
    points = vtk.vtkPoints()
    center_ids = []
    ring_ids = []

    for z in Z_LEVELS:
        center_ids.append(points.InsertNextPoint(0.0, 0.0, z))
        level_ring = []
        for i in range(N_SEG):
            theta = 2.0 * math.pi * i / N_SEG
            level_ring.append(points.InsertNextPoint(RADIUS * math.cos(theta), RADIUS * math.sin(theta), z))
        ring_ids.append(level_ring)

    return points, center_ids, ring_ids


def add_tet(grid, ids):
    tet = vtk.vtkTetra()
    for i, point_id in enumerate(ids):
        tet.GetPointIds().SetId(i, point_id)
    grid.InsertNextCell(tet.GetCellType(), tet.GetPointIds())


def add_tri(poly, ids):
    id_list = vtk.vtkIdList()
    id_list.SetNumberOfIds(3)
    for i, point_id in enumerate(ids):
        id_list.SetId(i, point_id)
    poly.InsertNextCell(vtk.VTK_TRIANGLE, id_list)


def add_cell_arrays(grid):
    fg = vtk.vtkDoubleArray()
    fg.SetName("Growth_Fg")
    fg.SetNumberOfComponents(9)

    c10 = vtk.vtkDoubleArray()
    c10.SetName("Material_C10")

    kpen = vtk.vtkDoubleArray()
    kpen.SetName("Material_Kpen")

    # Mild prescribed isotropic growth. C10 and Kpen are cell data to exercise
    # both external material mechanics and softened volumetric response.
    theta = 1.03
    fg_tuple = (theta, 0.0, 0.0, 0.0, theta, 0.0, 0.0, 0.0, theta)
    for _ in range(grid.GetNumberOfCells()):
        fg.InsertNextTuple(fg_tuple)
        c10.InsertNextValue(5.0e4)
        kpen.InsertNextValue(5.0e5)

    grid.GetCellData().AddArray(fg)
    grid.GetCellData().AddArray(c10)
    grid.GetCellData().AddArray(kpen)


def add_node_ids(data_set):
    node_ids = vtk.vtkIntArray()
    node_ids.SetName("GlobalNodeID")
    for i in range(data_set.GetNumberOfPoints()):
        node_ids.InsertNextValue(i + 1)
    data_set.GetPointData().AddArray(node_ids)


def add_element_ids(data_set):
    elem_ids = vtk.vtkIntArray()
    elem_ids.SetName("GlobalElementID")
    for i in range(data_set.GetNumberOfCells()):
        elem_ids.InsertNextValue(i + 1)
    data_set.GetCellData().AddArray(elem_ids)


def write_vtu(path, grid):
    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetFileName(str(path))
    writer.SetInputData(grid)
    writer.SetDataModeToAscii()
    writer.Write()


def write_vtp(path, poly):
    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName(str(path))
    writer.SetInputData(poly)
    writer.SetDataModeToAscii()
    writer.Write()


def make_face_poly(volume_points, global_tris):
    poly = vtk.vtkPolyData()
    points = vtk.vtkPoints()
    poly.Allocate(len(global_tris))

    global_to_local = {}
    local_to_global = []

    for tri in global_tris:
        local_tri = []
        for global_id in tri:
            if global_id not in global_to_local:
                point = volume_points.GetPoint(global_id)
                global_to_local[global_id] = points.InsertNextPoint(point)
                local_to_global.append(global_id)
            local_tri.append(global_to_local[global_id])
        add_tri(poly, local_tri)

    poly.SetPoints(points)

    node_ids = vtk.vtkIntArray()
    node_ids.SetName("GlobalNodeID")
    for global_id in local_to_global:
        node_ids.InsertNextValue(global_id + 1)
    poly.GetPointData().AddArray(node_ids)

    add_element_ids(poly)
    return poly


def main():
    SURF_DIR.mkdir(parents=True, exist_ok=True)

    points, center_ids, ring_ids = make_points()
    grid = vtk.vtkUnstructuredGrid()
    grid.SetPoints(points)

    bottom_tris = []
    top_tris = []
    wall_tris = []

    z0, z1 = 0, 1
    c0 = center_ids[z0]
    c1 = center_ids[z1]

    for i in range(N_SEG):
        j = (i + 1) % N_SEG
        r0i = ring_ids[z0][i]
        r0j = ring_ids[z0][j]
        r1i = ring_ids[z1][i]
        r1j = ring_ids[z1][j]

        # Split the triangular prism into three tetrahedra.
        add_tet(grid, [c0, r0i, r0j, c1])
        add_tet(grid, [r0i, r1i, r0j, c1])
        add_tet(grid, [r0j, r1i, r1j, c1])

        bottom_tris.append([c0, r0j, r0i])
        top_tris.append([c1, r1i, r1j])
        wall_tris.append([r0i, r1i, r1j])
        wall_tris.append([r0i, r1j, r0j])

    add_node_ids(grid)
    add_element_ids(grid)
    add_cell_arrays(grid)

    bottom = make_face_poly(points, bottom_tris)
    top = make_face_poly(points, top_tris)
    wall = make_face_poly(points, wall_tris)

    write_vtu(MESH_DIR / "mesh-complete.mesh.vtu", grid)
    write_vtp(SURF_DIR / "bottom.vtp", bottom)
    write_vtp(SURF_DIR / "top.vtp", top)
    write_vtp(SURF_DIR / "wall.vtp", wall)


if __name__ == "__main__":
    main()
