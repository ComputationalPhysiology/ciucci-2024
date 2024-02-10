from pathlib import Path

import dolfin
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np

from utils import try_except_runtimererror, name2latex
from geometry import get_lv_geometry


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


def load_lv_arrs(data_path, output, gammas, pressures, mesh_folder: Path = Path("meshes/lv")):
    print("Loading LV arrays")
    geo = get_lv_geometry(mesh_folder=mesh_folder)
    V_DG2 = dolfin.FunctionSpace(geo.mesh, "DG", 2)

    f = dolfin.Function(V_DG2)

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


def postprocess_lv(resultsdir, figdir, mesh_folder, print_stats=False, create_paraview=False):
    print("Postprocessing LV")
    output = Path(resultsdir) / "results.xdmf"

    gammas = np.load(resultsdir / "gammas.npy")
    pressures = np.load(resultsdir / "pressures.npy")

    figdir.mkdir(exist_ok=True, parents=True)

    data_path = resultsdir / "results.csv"
    if not data_path.is_file():
        load_lv_arrs(data_path, output, gammas, pressures, mesh_folder=mesh_folder)

    if print_stats:
        try:
            import polars as pl
        except ImportError:
            print("Install polars to print stats (pip install polars)")
            raise SystemExit(1)

        df = pl.read_csv(data_path)
        unloaded = df.filter(pl.col("pressure").eq(0.0)).filter(pl.col("gamma").eq(0.2))
        loaded = df.filter(pl.col("pressure").eq(15.0)).filter(pl.col("gamma").eq(0.2))

        print(unloaded.group_by("name").agg(pl.col("*").mean()))
        print(loaded.group_by("name").agg(pl.col("*").mean()))

        return

    if create_paraview:
        create_paraview_files(resultsdir, figdir=figdir, mesh_folder=mesh_folder)
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


def create_paraview_files(resultsdir, figdir, mesh_folder: Path = Path("meshes/lv")):
    print("Creating Paraview files")
    gammas = np.load(resultsdir / "gammas.npy")
    output = Path(resultsdir) / "results.xdmf"
    pvd_output = Path(figdir) / "pvd_files"
    pvd_output.mkdir(exist_ok=True, parents=True)
    geo = get_lv_geometry(mesh_folder=mesh_folder)

    V_DG2 = dolfin.FunctionSpace(geo.mesh, "DG", 2)
    V_CG2 = dolfin.VectorFunctionSpace(geo.mesh, "CG", 2)
    V_CG1 = dolfin.FunctionSpace(geo.mesh, "CG", 1)

    u = dolfin.Function(V_CG2)
    u.rename("u", "")
    f = dolfin.Function(V_DG2)

    with dolfin.XDMFFile(output.as_posix()) as xdmf:
        for ti in range(len(gammas)):
            # print(ti)
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
