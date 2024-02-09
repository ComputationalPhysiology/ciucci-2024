import dolfin
import cardiac_geometries
from pathlib import Path

import pulse
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

import numpy as np
import ufl_legacy as ufl
import utils


dolfin.parameters["form_compiler"]["quadrature_degree"] = 6
dolfin.parameters["form_compiler"]["cpp_optimize"] = True
dolfin.parameters["form_compiler"]["representation"] = "uflacs"
dolfin.parameters["form_compiler"]["optimize"] = True


def try_except_runtimererror(func, values, args):
    try:
        values.append(func(*args))
    except RuntimeError:
        pass


class SmoothLV(dolfin.UserExpression):
    def __init__(self, f):
        self.f = f
        super().__init__()

    def eval(self, value, x):
        dx = 0.05
        values = [self.f(x[0], x[1], x[2])]

        for i in [-2, -1, 1, 2]:
            try_except_runtimererror(self.f, values, (x[0] + i * dx, x[1], x[2]))
            try_except_runtimererror(self.f, values, (x[0], x[1] + i * dx, x[2]))
            try_except_runtimererror(self.f, values, (x[0], x[1], x[2] + i * dx))

        value[0] = np.mean(values)

    def value_shape(self):
        return ()


def get_geometry(mesh_folder: Path = Path("meshes/lv")):
    if not mesh_folder.is_dir():
        raise FileNotFoundError(f"Folder {mesh_folder} does not exist")

    return cardiac_geometries.geometry.Geometry.from_folder(mesh_folder)


def load_arrs(data_path, output, gammas, pressures, mesh_folder: Path = Path("meshes/lv")):
    geo = get_geometry(mesh_folder=mesh_folder)
    V_DG1 = dolfin.FunctionSpace(geo.mesh, "DG", 1)

    f = dolfin.Function(V_DG1)
    from postprocess import name2latex

    data = []
    with dolfin.XDMFFile(output.as_posix()) as xdmf:
        for ti in range(len(gammas)):
            # xdmf.read_checkpoint(u, "u", ti)
            for name in [
                "sigma_ff",
                "sigma_ss",
                "sigma_nn",
                "sigma_dev_ff",
                "sigma_dev_ss",
                "sigma_dev_nn",
                "E_ff",
                "E_ss",
                "E_nn",
                "p",
            ]:
                xdmf.read_checkpoint(f, name, ti)
                f_arr = f.vector().get_local()

                data.extend(
                    [
                        {
                            "name": name,
                            "value": fi,
                            "gamma": gammas[ti],
                            "pressure": pressures[ti],
                            "latex": name2latex(name),
                        }
                        for fi in f_arr
                    ]
                )
    df = pd.DataFrame(data)
    df.to_csv(data_path)


def postprocess(resultsdir, figdir, print_stats=False):
    output = Path(resultsdir) / "results.xdmf"

    gammas = np.load(resultsdir / "gammas.npy")
    pressures = np.load(resultsdir / "pressures.npy")

    figdir.mkdir(exist_ok=True, parents=True)

    data_path = resultsdir / "results.csv"
    if not data_path.is_file():
        load_arrs(data_path, output, gammas, pressures)

    if print_stats:
        import polars as pl

        df = pl.read_csv(data_path)
        unloaded = df.filter(pl.col("pressure").eq(0.0)).filter(pl.col("gamma").eq(0.2))
        loaded = df.filter(pl.col("pressure").eq(15.0)).filter(pl.col("gamma").eq(0.2))

        print(unloaded.group_by("name").agg(pl.col("*").mean()))
        print(loaded.group_by("name").agg(pl.col("*").mean()))

        return

    df = pd.read_csv(data_path)

    target_gamma = 0.2
    df_unloaded = df[np.isclose(df["pressure"], 0.0) & np.isclose(df["gamma"], target_gamma)]
    df_unloaded = df_unloaded.assign(label="Unloaded systole\nESP = 0 kPa")

    traget_pressure = 15.0
    df_loaded = df[
        np.isclose(df["pressure"], traget_pressure) & np.isclose(df["gamma"], target_gamma)
    ]
    df_loaded = df_loaded.assign(label="Standard systole\nESP = 15 kPa")
    df1 = pd.concat([df_loaded, df_unloaded])

    df1_dev_stress = df1[df1["name"].isin(["sigma_dev_ff", "sigma_dev_ss", "sigma_dev_nn", "p"])]
    plt.rcParams.update({"font.size": 16})
    fig = plt.figure()

    ax = sns.barplot(
        data=df1_dev_stress,
        x="label",
        y="value",
        hue="latex",
        errorbar="ci",
        alpha=0.7,
    )
    ax.get_legend().set_title(None)
    ax.set_xlabel("")
    ax.set_ylabel("Average stress [kPa]")
    ax.grid()
    fig.tight_layout()
    fig.savefig(figdir / "stress_dev.svg")  # type: ignore
    plt.close(fig)

    df1_stress = df1[df1["name"].isin(["sigma_ff", "sigma_ss", "sigma_nn"])]
    plt.rcParams.update({"font.size": 16})
    fig = plt.figure()

    ax = sns.barplot(
        data=df1_stress,
        x="label",
        y="value",
        hue="latex",
        errorbar="ci",
        alpha=0.7,
    )
    ax.get_legend().set_title(None)
    ax.set_xlabel("")
    ax.set_ylabel("Average stress [kPa]")
    ax.grid()
    fig.tight_layout()
    fig.savefig(figdir / "stress.svg")  # type: ignore
    plt.close(fig)

    df1_strain = df1[df1["name"].isin(["E_ff", "E_ss", "E_nn"])]
    fig = plt.figure()
    ax = sns.barplot(
        data=df1_strain,
        x="label",
        y="value",
        hue="latex",
        alpha=0.7,
    )
    sns.move_legend(ax, "lower center", bbox_to_anchor=(0.5, 1), ncol=3, title=None, frameon=False)
    ax.set_xlabel("")
    ax.set_ylabel("Average strain")
    ax.grid()
    fig.savefig(figdir / "strain.svg")  # type: ignore
    plt.close(fig)


def create_paraview_files(resultsdir, mesh_folder: Path = Path("meshes/lv")):
    gammas = np.load(resultsdir / "gammas.npy")
    output = Path(resultsdir) / "results.xdmf"
    pvd_output = Path(resultsdir) / "pvd_files"
    pvd_output.mkdir(exist_ok=True, parents=True)
    geo = get_geometry(mesh_folder=mesh_folder)

    V_DG1 = dolfin.FunctionSpace(geo.mesh, "DG", 1)
    V_CG2 = dolfin.VectorFunctionSpace(geo.mesh, "CG", 2)
    V_CG1 = dolfin.FunctionSpace(geo.mesh, "CG", 1)

    u = dolfin.Function(V_CG2)
    u.rename("u", "")
    f = dolfin.Function(V_DG1)

    with dolfin.XDMFFile(output.as_posix()) as xdmf:
        for ti in range(len(gammas)):
            print(ti)
            xdmf.read_checkpoint(u, "u", ti)

            for i, name in enumerate(
                [
                    "sigma_ff",
                    "sigma_ss",
                    "sigma_nn",
                ],
            ):
                xdmf.read_checkpoint(f, name, ti)
                f_int = dolfin.interpolate(SmoothLV(f), V_CG1)
                f_int.rename(name, "")
                with dolfin.XDMFFile((pvd_output / f"{name}_{ti}.xdmf").as_posix()) as xdmf2:
                    xdmf2.parameters["functions_share_mesh"] = True
                    xdmf2.parameters["flush_output"] = True

                    xdmf2.write(u, ti)
                    xdmf2.write(f_int, ti)


def main(output_folder, mesh_folder: Path = Path("meshes/lv")):
    overwrite = True
    Path(output_folder).mkdir(exist_ok=True, parents=True)
    output = Path(output_folder) / "results.xdmf"
    if output.is_file() and not overwrite:
        print(f"Output {output} already exists")
        return

    geo = get_geometry(mesh_folder=mesh_folder)

    V_DG1 = dolfin.FunctionSpace(geo.mesh, "DG", 2)
    proj = utils.Projector(V_DG1)
    sigma_ff = dolfin.Function(V_DG1)
    sigma_ss = dolfin.Function(V_DG1)
    sigma_nn = dolfin.Function(V_DG1)

    sigma_dev_ff = dolfin.Function(V_DG1)
    sigma_dev_ss = dolfin.Function(V_DG1)
    sigma_dev_nn = dolfin.Function(V_DG1)

    E_ff = dolfin.Function(V_DG1)
    E_ss = dolfin.Function(V_DG1)
    E_nn = dolfin.Function(V_DG1)

    von_Mises = dolfin.Function(V_DG1)

    microstructure = pulse.Microstructure(f0=geo.f0, s0=geo.s0, n0=geo.n0)

    geometry = pulse.HeartGeometry(
        mesh=geo.mesh,
        markers=geo.markers,
        marker_functions=pulse.MarkerFunctions(ffun=geo.ffun),
        microstructure=microstructure,
    )

    m2mm = 1000.0

    matparams = {
        "a": 2.280,
        "b": 9.726,
        "a_f": 1.685,
        "b_f": 15.779,
        "a_s": 0.0,
        "b_s": 0.0,
        "a_fs": 0.0,
        "b_fs": 0.0,
    }
    gamma = dolfin.Constant(0.0)
    activation = gamma

    material = pulse.HolzapfelOgden(
        active_model="active_strain",
        activation=activation,
        parameters=matparams,
        f0=geometry.f0,
        s0=geometry.s0,
        n0=geometry.n0,
    )

    # Pericardium type Robin BC
    spring = dolfin.Constant(500 / m2mm)  # kPa/mm
    robin_bc = [
        pulse.RobinBC(value=dolfin.Constant(spring), marker=geo.markers["EPI"][0]),
    ]

    # LV Pressure
    lvp = dolfin.Constant(0.0)
    lv_marker = geometry.markers["ENDO"][0]
    lv_pressure = pulse.NeumannBC(traction=lvp, marker=lv_marker, name="lv")
    neumann_bc = [lv_pressure]

    # Fix the basal plane in the longitudinal direction
    def fix_basal_plane(W):
        V = W if W.sub(0).num_sub_spaces() == 0 else W.sub(0)
        bc = dolfin.DirichletBC(
            V.sub(0),
            dolfin.Constant(0.0),
            geometry.ffun,
            geometry.markers["BASE"][0],
        )
        return bc

    dirichlet_bc = [fix_basal_plane]

    bcs = pulse.BoundaryConditions(dirichlet=dirichlet_bc, neumann=neumann_bc, robin=robin_bc)

    problem = pulse.MechanicsProblem(geometry, material, bcs)
    gammas = [0.0, 0.2, 0.2]
    pressures = [0.0, 0.0, 15.0]

    np.save(output_folder / "gammas.npy", gammas)
    np.save(output_folder / "pressures.npy", pressures)

    for ti, (lvp_, g) in enumerate(zip(pressures, gammas)):
        print(lvp_, g)
        pulse.iterate.iterate(
            problem,
            (lvp, gamma),
            (lvp_, g),
            initial_number_of_steps=20,
            continuation=False,
        )
        u, p = problem.state.split()

        F = pulse.kinematics.DeformationGradient(u)
        J = ufl.det(F)
        sigma = material.CauchyStress(F, p)
        sigma_dev = sigma - (1 / 3) * ufl.tr(sigma) * ufl.Identity(3)
        E = pulse.kinematics.GreenLagrangeStrain(F)
        f = F * geo.f0
        s = F * geo.s0
        n = F * geo.n0

        proj.project(von_Mises, utils.von_mises(sigma))

        proj.project(sigma_ff, dolfin.inner(f, sigma * f))
        proj.project(sigma_ss, dolfin.inner(s, sigma * s))
        proj.project(sigma_nn, dolfin.inner(n, sigma * n))

        proj.project(sigma_dev_ff, dolfin.inner(f, sigma_dev * f))
        proj.project(sigma_dev_ss, dolfin.inner(s, sigma_dev * s))
        proj.project(sigma_dev_nn, dolfin.inner(n, sigma_dev * n))

        proj.project(E_ff, dolfin.inner(geo.f0, E * geo.f0))
        proj.project(E_ss, dolfin.inner(geo.s0, E * geo.s0))
        proj.project(E_nn, dolfin.inner(geo.n0, E * geo.n0))
        with dolfin.XDMFFile(output.as_posix()) as f:
            f.parameters["functions_share_mesh"] = True
            f.parameters["rewrite_function_mesh"] = False
            f.write_checkpoint(
                u,
                function_name="u",
                time_step=ti,
                encoding=dolfin.XDMFFile.Encoding.HDF5,
                append=True,
            )
            f.write_checkpoint(
                sigma_ff,
                function_name="sigma_ff",
                time_step=ti,
                encoding=dolfin.XDMFFile.Encoding.HDF5,
                append=True,
            )
            f.write_checkpoint(
                sigma_ss,
                function_name="sigma_ss",
                time_step=ti,
                encoding=dolfin.XDMFFile.Encoding.HDF5,
                append=True,
            )
            f.write_checkpoint(
                sigma_nn,
                function_name="sigma_nn",
                time_step=ti,
                encoding=dolfin.XDMFFile.Encoding.HDF5,
                append=True,
            )
            f.write_checkpoint(
                sigma_dev_ff,
                function_name="sigma_dev_ff",
                time_step=ti,
                encoding=dolfin.XDMFFile.Encoding.HDF5,
                append=True,
            )
            f.write_checkpoint(
                sigma_dev_ss,
                function_name="sigma_dev_ss",
                time_step=ti,
                encoding=dolfin.XDMFFile.Encoding.HDF5,
                append=True,
            )
            f.write_checkpoint(
                sigma_dev_nn,
                function_name="sigma_dev_nn",
                time_step=ti,
                encoding=dolfin.XDMFFile.Encoding.HDF5,
                append=True,
            )
            f.write_checkpoint(
                E_ff,
                function_name="E_ff",
                time_step=ti,
                encoding=dolfin.XDMFFile.Encoding.HDF5,
                append=True,
            )
            f.write_checkpoint(
                E_ss,
                function_name="E_ss",
                time_step=ti,
                encoding=dolfin.XDMFFile.Encoding.HDF5,
                append=True,
            )
            f.write_checkpoint(
                E_nn,
                function_name="E_nn",
                time_step=ti,
                encoding=dolfin.XDMFFile.Encoding.HDF5,
                append=True,
            )
            f.write_checkpoint(
                von_Mises,
                function_name="von_Mises",
                time_step=ti,
                encoding=dolfin.XDMFFile.Encoding.HDF5,
                append=True,
            )
            f.write_checkpoint(
                p,
                function_name="p",
                time_step=ti,
                encoding=dolfin.XDMFFile.Encoding.HDF5,
                append=True,
            )


if __name__ == "__main__":
    resultsdir = Path("results/lv_incomp")
    resultsdir.mkdir(exist_ok=True, parents=True)
    # main(resultsdir=resultsdir)
    figdir = Path("figures") / "lv_incomp"
    postprocess(resultsdir=resultsdir, figdir=figdir, print_stats=False)
    # create_paraview_files(resultsdir=resultsdir)
