from pathlib import Path
from mpi4py import MPI
import adios4dolfinx
import dolfinx
import matplotlib.pyplot as plt
import pyvista as pv


COLORMAP = "cividis"
CRINKLE = True
SHOW_EDGES = True


def plot_native_ES():
    resultsdir = Path("results") / "native"
    imagefile = Path("../data") / "native.png"
    figname = Path("figures") / "native_mesh.svg"

    comm = MPI.COMM_WORLD
    mesh = adios4dolfinx.read_mesh_from_legacy_h5(
        resultsdir / "results_current.xdmf", comm, "/sigma_ff/sigma_ff_2/mesh"
    )

    V = dolfinx.fem.functionspace(mesh, ("DG", 0))

    cells, types, x = dolfinx.plot.vtk_mesh(mesh)
    grid = pv.UnstructuredGrid(cells, types, x)

    sigma = dolfinx.fem.Function(V, dtype=dolfinx.default_real_type())

    adios4dolfinx.read_function_from_legacy_h5(
        resultsdir / "results_current.xdmf", mesh.comm, sigma, group="sigma_ff", step=2
    )
    grid.cell_data["sigma_ff"] = sigma.x.array
    # grid.point_data["sigma_ff"] = sigma.x.array
    grid.set_active_scalars("sigma_ff")

    s = 125.0
    transform = pv.Transform().scale(s, s, s)
    grid_trans = transform.apply(grid)
    grid_trans.translate((910.0, 1000, 0), inplace=True)

    image = pv.read(imagefile)
    bounds = [200, 1800, 400, 1500, 0, 0]
    image = image.clip_box(bounds, invert=False)

    # # Create a PyVista plotter
    pv.start_xvfb()
    plotter = pv.Plotter()
    plotter.add_mesh(image, cmap="gray", show_scalar_bar=False)
    cmap = plt.get_cmap(COLORMAP)
    plotter.add_mesh_clip_plane(
        grid_trans,
        normal=(0, -1, 0),
        crinkle=CRINKLE,
        # tubing=False,
        outline_opacity=False,
        show_scalar_bar=True,
        cmap=cmap,
        show_edges=SHOW_EDGES,
        lighting=False,
        # opacity=0.5,
        scalar_bar_args={
            "color": "black",
            "title": "Fiber stress [kPa]",
            "title_font_size": 30,
            "label_font_size": 30,
            "height": 0.1,
            "width": 0.8,
            "vertical": False,
            "position_x": 0.1,
            "position_y": 0.85,
        },
        clim=[0, 80.0],
    )
    widget = plotter.plane_widgets[0]
    widget.SetEnabled(not widget.GetEnabled())
    # plotter.camera_position = [
    #     (2873.8774814991825, 2647.8774814991825, 1706.8096934058667),
    #     (1167.5, 941.5, 0.4322119066891048),
    #     (0.0, 0.0, 1.0),
    # ]
    plotter.camera_position = [
        (-27.478591656723232, 2730.8560998239436, 620.045130319616),
        (1000.0, 950.0, -0.4971771240234375),
        (0.24374756816294021, -0.19071166164668962, 0.9509028263322237),
    ]

    plotter.save_graphic(figname)


def plot_transplanted_ES():
    resultsdir = Path("results") / "transplanted"
    imagefile = Path("../data") / "transplanted.png"
    figname = Path("figures") / "transplanted_mesh.svg"

    comm = MPI.COMM_WORLD
    mesh = adios4dolfinx.read_mesh_from_legacy_h5(
        resultsdir / "results_current.xdmf", comm, "/sigma_ff/sigma_ff_2/mesh"
    )

    V = dolfinx.fem.functionspace(mesh, ("DG", 0))

    cells, types, x = dolfinx.plot.vtk_mesh(mesh)
    grid = pv.UnstructuredGrid(cells, types, x)

    sigma = dolfinx.fem.Function(V, dtype=dolfinx.default_real_type())

    adios4dolfinx.read_function_from_legacy_h5(
        resultsdir / "results_current.xdmf", mesh.comm, sigma, group="sigma_ff", step=2
    )
    grid.cell_data["sigma_ff"] = sigma.x.array
    grid.set_active_scalars("sigma_ff")

    s = 100.0
    transform = pv.Transform().scale(s, s, s)
    grid_trans = transform.apply(grid)
    grid_trans.translate((670.0, 670.0, 0), inplace=True)

    image = pv.read(imagefile)
    bounds = [100, 1200, 300, 920, 0, 0]
    image = image.clip_box(bounds, invert=False)

    # # Create a PyVista plotter
    pv.start_xvfb()
    plotter = pv.Plotter()
    plotter.add_mesh(image, cmap="gray", show_scalar_bar=False)
    cmap = plt.get_cmap(COLORMAP)
    plotter.add_mesh_clip_plane(
        grid_trans,
        normal=(0, -1, 0),
        crinkle=CRINKLE,
        show_edges=SHOW_EDGES,
        # tubing=False,
        outline_opacity=False,
        show_scalar_bar=True,
        cmap=cmap,
        lighting=False,
        # opacity=0.5,
        scalar_bar_args={
            "color": "black",
            "title": "Fiber stress [kPa]",
            "title_font_size": 30,
            "label_font_size": 30,
            "height": 0.1,
            "width": 0.8,
            "vertical": False,
            "position_x": 0.1,
            "position_y": 0.85,
        },
        clim=[0, 80.0],
    )
    widget = plotter.plane_widgets[0]
    widget.SetEnabled(not widget.GetEnabled())
    # plotter.camera_position = [
    #     (1965.226254107329, 1727.2262541073258, 1188.739149458118),
    #     (776.5, 538.5, 0.012895350782002879),
    #     (-0.4082482904638644, -0.4082482904638633, 0.8164965809277254),
    # ]

    plotter.camera_position = [
        (-97.50252570755485, 1647.5530835985949, 401.1000330530449),
        (776.5, 538.5, 0.012895350782002879),
        (0.14421878471400296, -0.23405925157590624, 0.9614661766735961),
    ]
    plotter.save_graphic(figname)


if __name__ == "__main__":
    plot_native_ES()
    plot_transplanted_ES()
