from pathlib import Path
import typing
import json

import dolfin
import meshio
import pulse
import cardiac_geometries
import gmsh
import numpy as np


def create_cylinder_mesh(
    mesh_folder: Path, char_length: float = 50.0, L: float = 4000, r: float = 300
) -> None:
    Path(mesh_folder).mkdir(exist_ok=True, parents=True)
    gmsh.initialize()
    gmsh.model.add("Cell")

    cylinder = gmsh.model.occ.addCylinder(-L / 2, 0.0, 0.0, L, 0, 0, r)
    gmsh.model.occ.synchronize()

    surfaces = gmsh.model.occ.getEntities(dim=2)
    # inlet_marker, outlet_marker, wall_marker, obstacle_marker = 1, 3, 5, 7
    right = 1
    left = 2
    sides = 3

    for surface in surfaces:
        com = gmsh.model.occ.getCenterOfMass(surface[0], surface[1])
        print(surface, com)
        if np.allclose(com, [L / 2, 0, 0]):
            # Right boundary
            gmsh.model.addPhysicalGroup(surface[0], [surface[1]], right)
            gmsh.model.setPhysicalName(surface[0], right, "Right")
        elif np.allclose(com, [-L / 2, 0, 0]):
            # Left boundary
            gmsh.model.addPhysicalGroup(surface[0], [surface[1]], left)
            gmsh.model.setPhysicalName(surface[0], left, "Left")
        else:
            # Sides
            gmsh.model.addPhysicalGroup(surface[0], [surface[1]], sides)
            gmsh.model.setPhysicalName(surface[0], sides, "Sides")

    gmsh.model.addPhysicalGroup(3, [cylinder], sides)
    gmsh.model.setPhysicalName(3, cylinder, "Volume")

    gmsh.option.setNumber("Mesh.CharacteristicLengthMin", char_length)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", char_length)

    gmsh.model.mesh.generate(3)
    gmsh.model.mesh.optimize("Netgen")
    gmsh.write((Path(mesh_folder) / "cell.msh").as_posix())
    gmsh.finalize()


def preprocess_lv(
    mesh_folder: Path,
    r_short_endo=4.0,
    r_long_endo=8.0,
    r_short_epi=5.5,
    r_long_epi=9.5,
    psize_ref=0.5,
    create_fibers=True,
    fiber_space="Quadrature_6",
    **kwargs,
):
    cardiac_geometries.create_lv_ellipsoid(
        mesh_folder,
        r_short_endo=r_short_endo,
        r_long_endo=r_long_endo,
        r_short_epi=r_short_epi,
        r_long_epi=r_long_epi,
        psize_ref=psize_ref,
        create_fibers=create_fibers,
        fiber_space=fiber_space,
    )


def get_lv_geometry(mesh_folder: Path = Path("meshes/lv")):
    if not mesh_folder.is_dir():
        raise FileNotFoundError(f"Folder {mesh_folder} does not exist")

    return cardiac_geometries.geometry.Geometry.from_folder(mesh_folder)


class Geometry(typing.NamedTuple):
    mesh: dolfin.Mesh
    markers: dict[str, list[int]]
    marker_functions: pulse.MarkerFunctions
    rad0: dolfin.Function
    circ0: dolfin.Function
    long0: dolfin.Function

    def ds(self, marker: int) -> dolfin.Measure:
        return dolfin.Measure(
            "ds",
            domain=self.mesh,
            subdomain_data=self.marker_functions.vfun,
            subdomain_id=marker,
        )

    def dx(self, marker: int) -> dolfin.Measure:
        return dolfin.Measure(
            "dx",
            domain=self.mesh,
            subdomain_id=marker,
        )


def create_mesh(mesh, cell_type):
    # From http://jsdokken.com/converted_files/tutorial_pygmsh.html
    cells = mesh.get_cells_type(cell_type)
    if cells.size == 0:
        return None

    cell_data = mesh.get_cell_data("gmsh:physical", cell_type)
    out_mesh = meshio.Mesh(
        points=mesh.points,
        cells={cell_type: cells},
        cell_data={"name_to_read": [cell_data]},
    )
    return out_mesh


def read_meshfunction(fname, obj):
    try:
        with dolfin.XDMFFile(Path(fname).as_posix()) as f:
            f.read(obj, "name_to_read")
    except RuntimeError:
        pass


def preprocess_cylinder(
    mesh_folder: Path,
    char_length: float = 50.0,
    length: float = 4000,
    radius: float = 300,
    unlink: bool = False,
) -> Geometry:
    outdir = Path(mesh_folder)
    msh_file = outdir / "cell.msh"
    if not msh_file.exists():
        create_cylinder_mesh(
            mesh_folder,
            char_length=char_length,
            L=length,
            r=radius,
        )

    triangle_mesh_name = outdir / "triangle_mesh.xdmf"
    tetra_mesh_name = outdir / "mesh.xdmf"
    marker_name = outdir / "markers.json"

    if dolfin.MPI.comm_world.size == 1:
        # Only use gmsh when running in serial
        msh = meshio.gmsh.read(msh_file)

        outdir.mkdir(exist_ok=True, parents=True)
        triangle_mesh = create_mesh(msh, "triangle")
        tetra_mesh = create_mesh(msh, "tetra")

        if triangle_mesh is not None:
            meshio.write(triangle_mesh_name, triangle_mesh)

        if tetra_mesh is None:
            raise RuntimeError("Unable to create mesh")

        meshio.write(
            tetra_mesh_name,
            tetra_mesh,
        )
        markers = {k: [int(vi) for vi in v] for k, v in msh.field_data.items()}
        marker_name.write_text(json.dumps(markers))

    markers = json.loads(marker_name.read_text())
    mesh = dolfin.Mesh()

    with dolfin.XDMFFile(tetra_mesh_name.as_posix()) as infile:
        infile.read(mesh)

    ffun_val = dolfin.MeshValueCollection("size_t", mesh, 2)
    read_meshfunction(triangle_mesh_name, ffun_val)
    ffun = dolfin.MeshFunction("size_t", mesh, ffun_val)
    ffun.array()[ffun.array() == max(ffun.array())] = 0
    if unlink:
        triangle_mesh_name.unlink(missing_ok=True)
        triangle_mesh_name.with_suffix(".h5").unlink(missing_ok=True)
    else:
        ffun_path = outdir / "ffun.xdmf"
        with dolfin.XDMFFile(ffun_path.as_posix()) as infile:
            infile.write(ffun)

    Q = pulse.QuadratureSpace(mesh, degree=6, dim=3)
    rad0 = dolfin.Function(Q)
    rad0.interpolate(
        dolfin.Expression(
            (
                "0",
                "x[1]/sqrt(x[1]*x[1] + x[2]*x[2])",
                "x[2]/sqrt(x[1]*x[1] + x[2]*x[2])",
            ),
            degree=1,
        )
    )
    circ0 = dolfin.Function(Q)
    circ0.interpolate(
        dolfin.Expression(
            (
                "0",
                "x[2]/sqrt(x[1]*x[1] + x[2]*x[2])",
                "-x[1]/sqrt(x[1]*x[1] + x[2]*x[2])",
            ),
            degree=1,
        )
    )
    long0 = dolfin.Function(Q)
    long0.interpolate(dolfin.Expression(("1", "0", "0"), degree=1))

    import ldrb

    ldrb.fun_to_xdmf(rad0, outdir / "rad0")
    ldrb.fun_to_xdmf(circ0, outdir / "circ0")
    ldrb.fun_to_xdmf(long0, outdir / "long0")

    return Geometry(
        mesh=mesh,
        markers=markers,
        marker_functions=pulse.MarkerFunctions(ffun=ffun),
        rad0=rad0,
        circ0=circ0,
        long0=long0,
    )


def get_cylinder_geometry(mesh_folder: Path) -> Geometry:
    triangle_mesh_name = mesh_folder / "triangle_mesh.xdmf"
    tetra_mesh_name = mesh_folder / "mesh.xdmf"
    marker_name = mesh_folder / "markers.json"
    ffun_path = mesh_folder / "ffun.xdmf"

    markers = json.loads(marker_name.read_text())
    mesh = dolfin.Mesh()

    with dolfin.XDMFFile(tetra_mesh_name.as_posix()) as infile:
        infile.read(mesh)

    ffun_val = dolfin.MeshValueCollection("size_t", mesh, 2)
    read_meshfunction(triangle_mesh_name, ffun_val)
    ffun = dolfin.MeshFunction("size_t", mesh, ffun_val)
    ffun.array()[ffun.array() == max(ffun.array())] = 0
    ffun_path = mesh_folder / "ffun.xdmf"
    with dolfin.XDMFFile(ffun_path.as_posix()) as infile:
        infile.write(ffun)

    Q = pulse.QuadratureSpace(mesh, degree=6, dim=3)
    rad0 = dolfin.Function(Q)
    rad0.interpolate(
        dolfin.Expression(
            (
                "0",
                "x[1]/sqrt(x[1]*x[1] + x[2]*x[2])",
                "x[2]/sqrt(x[1]*x[1] + x[2]*x[2])",
            ),
            degree=1,
        )
    )
    circ0 = dolfin.Function(Q)
    circ0.interpolate(
        dolfin.Expression(
            (
                "0",
                "x[2]/sqrt(x[1]*x[1] + x[2]*x[2])",
                "-x[1]/sqrt(x[1]*x[1] + x[2]*x[2])",
            ),
            degree=1,
        )
    )
    long0 = dolfin.Function(Q)
    long0.interpolate(dolfin.Expression(("1", "0", "0"), degree=1))
    return Geometry(
        mesh=mesh,
        markers=markers,
        marker_functions=pulse.MarkerFunctions(ffun=ffun),
        rad0=rad0,
        circ0=circ0,
        long0=long0,
    )
