from pathlib import Path
from mpi4py import MPI
import adios4dolfinx
import dolfinx
import matplotlib.pyplot as plt
import pyvista as pv


def plot_native_ES():
    resultsdir = Path("results") / "native"
    imagefile = Path("../data") / "native.png"
    figname = Path("figures") / "native_slice.svg"

    comm = MPI.COMM_WORLD
    mesh = adios4dolfinx.read_mesh_from_legacy_h5(
        resultsdir / "results_current_smooth.xdmf", comm, "/sigma_ff/sigma_ff_2/mesh"
    )

    V = dolfinx.fem.functionspace(mesh, ("DG", 1))

    cells, types, x = dolfinx.plot.vtk_mesh(V)
    grid = pv.UnstructuredGrid(cells, types, x)

    sigma = dolfinx.fem.Function(V, dtype=dolfinx.default_real_type())

    adios4dolfinx.read_function_from_legacy_h5(
        resultsdir / "results_current_smooth.xdmf", mesh.comm, sigma, group="sigma_ff", step=2
    )
    grid.point_data["sigma_ff"] = sigma.x.array
    grid.set_active_scalars("sigma_ff")

    single_slice = grid.slice(normal=[0, 0, 1])

    # # Create a PyVista plotter
    pv.start_xvfb()
    plotter = pv.Plotter()
    plotter.add_background_image(imagefile)
    # cmap = plt.get_cmap("inferno")
    # cmap = plt.get_cmap("plasma")
    cmap = plt.get_cmap("jet")
    single_slice.translate((-2.4, 0.6, 0), inplace=True)
    plotter.add_mesh(
        single_slice,
        show_scalar_bar=True,
        cmap=cmap,
        lighting=False,
        opacity=0.5,
        scalar_bar_args={"color": "white", "title": "Fiber stress", "vertical": True},
        clim=[0, 100],
    )
    plotter.camera_position = [
        (-4.467109849948486, -1.9749349095021953, 26.810133660868157),
        (-0.683612715915249, 0.01947981678538957, 0.003457695253512795),
        (-0.0979259048249759, 0.9933781711560151, 0.06008599033835282),
    ]

    plotter.save_graphic(figname)


def plot_transplanted_ES():
    resultsdir = Path("results") / "transplanted"
    imagefile = Path("../data") / "transplanted.png"
    figname = Path("figures") / "transplanted_slice.svg"

    comm = MPI.COMM_WORLD
    mesh = adios4dolfinx.read_mesh_from_legacy_h5(
        resultsdir / "results_current_smooth.xdmf", comm, "/sigma_ff/sigma_ff_2/mesh"
    )

    V = dolfinx.fem.functionspace(mesh, ("DG", 1))

    cells, types, x = dolfinx.plot.vtk_mesh(V)
    grid = pv.UnstructuredGrid(cells, types, x)

    sigma = dolfinx.fem.Function(V, dtype=dolfinx.default_real_type())

    adios4dolfinx.read_function_from_legacy_h5(
        resultsdir / "results_current_smooth.xdmf", mesh.comm, sigma, group="sigma_ff", step=2
    )
    grid.point_data["sigma_ff"] = sigma.x.array
    grid.set_active_scalars("sigma_ff")

    single_slice = grid.slice(normal=[0, 0, 1])

    # # Create a PyVista plotter
    pv.start_xvfb()
    plotter = pv.Plotter()
    plotter.add_background_image(imagefile)
    # cmap = plt.get_cmap("inferno")
    # cmap = plt.get_cmap("plasma")
    cmap = plt.get_cmap("jet")
    single_slice.translate((-2.3, 1.1, 0), inplace=True)
    plotter.add_mesh(
        single_slice,
        show_scalar_bar=True,
        cmap=cmap,
        lighting=False,
        opacity=0.5,
        scalar_bar_args={"color": "white", "title": "Fiber stress", "vertical": True},
        clim=[0, 1.0],
    )
    plotter.camera_position = [
        (-3.780809292091852, -1.6131616899756525, 21.947585595483087),
        (-0.683612715915249, 0.01947981678538957, 0.003457695253512795),
        (-0.0979259048249759, 0.9933781711560151, 0.06008599033835282),
    ]

    plotter.save_graphic(figname)


if __name__ == "__main__":
    # plot_native_ES()
    plot_transplanted_ES()
