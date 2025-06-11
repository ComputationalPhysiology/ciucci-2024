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
    grid_trans.translate((840.0, 1000, 0), inplace=True)

    image = pv.read(imagefile)
    bounds = [150, 1880, 380, 1680, 0, 0]
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
        (-391.72899557008725, 3105.815945921396, 482.1905991953405),
        (1000.0, 950.0, -0.4971771240234375),
        (0.11554103254627512, -0.14541487259854116, 0.9826010302382604),
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

    s = 125.0

    x_offset = 70
    y_offset = 70

    transform = pv.Transform().scale(s, s, s)
    grid_trans = transform.apply(grid)
    grid_trans.translate((930.0 + x_offset, 970.0 + y_offset, 0), inplace=True)

    image = pv.read(imagefile)
    bounds = [80, 1240, 50, 930, 0, 0]
    image = image.clip_box(bounds, invert=False)
    transform_img = pv.Transform().scale(1.55, 1.55, 1.55)
    image_trans = transform_img.apply(image)
    image_trans.translate((x_offset, y_offset, 0), inplace=True)

    # # Create a PyVista plotter
    pv.start_xvfb()
    plotter = pv.Plotter()
    plotter.add_mesh(image_trans, cmap="gray", show_scalar_bar=False)
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

    # plotter.camera_position = [
    #     (40.96161473111009, 2162.8016595842296, 593.1096355196673),
    #     (660.0, 490.0, 0.1694488525390625),
    #     (0.14390229106881514, -0.28288076991887, 0.948298792910364),
    # ]
    plotter.camera_position = [
        (-391.72899557008725, 3105.815945921396, 482.1905991953405),
        (1000.0, 950.0, -0.4971771240234375),
        (0.11554103254627512, -0.14541487259854116, 0.9826010302382604),
    ]
    plotter.save_graphic(figname)


if __name__ == "__main__":
    plot_native_ES()
    # plot_transplanted_ES()
