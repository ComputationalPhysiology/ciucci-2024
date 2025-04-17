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


def load_lv_arrs(
    data_path, output, gammas, pressures, volumes, mesh_folder: Path = Path("meshes/lv")
):
    print("Loading LV arrays")
    geo = get_lv_geometry(mesh_folder=mesh_folder)
    V_DG2 = dolfin.FunctionSpace(geo.mesh, "DG", 1)
    V_CG1 = dolfin.FunctionSpace(geo.mesh, "CG", 1)

    f_ = dolfin.Function(V_DG2)
    p = dolfin.Function(V_CG1)

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
                f = p if name == "p" else f_
                xdmf.read_checkpoint(f, name, ti)
                f_arr = f.vector().get_local()

                data.extend(
                    [
                        {
                            "time": ti,
                            "name": name,
                            "value": fi,
                            "gamma": gammas[ti],
                            "pressure": pressures[ti],
                            "volume": volumes[ti],
                            "latex": name2latex(name),
                        }
                        for fi in f_arr
                    ]
                )

    df = pd.DataFrame(data)
    df.to_csv(data_path)


def postprocess_lv(resultsdir, figdir, mesh_folder, print_stats=False):
    print("Postprocessing LV")
    output = Path(resultsdir) / "results_reference.xdmf"

    gammas = np.load(resultsdir / "gammas.npy")
    pressures = np.load(resultsdir / "pressures.npy")
    volumes = np.load(resultsdir / "volumes.npy")
    figdir.mkdir(exist_ok=True, parents=True)

    data_path = resultsdir / "results.csv"
    if data_path.is_file():
        load_lv_arrs(data_path, output, gammas, pressures, volumes, mesh_folder=mesh_folder)

    if print_stats:
        try:
            import polars as pl
        except ImportError:
            print("Install polars to print stats (pip install polars)")
            raise SystemExit(1)

        df = pl.read_csv(data_path)

        ED = df.filter(pl.col("time").eq(1))
        ES = df.filter(pl.col("time").eq(2))
        print(mesh_folder)
        print("ED")
        print(
            ED.group_by("name").agg(pl.col("*").mean())[
                ["name", "value", "pressure", "volume", "gamma"]
            ]
        )
        print("ES")
        print(
            ES.group_by("name").agg(pl.col("*").mean())[
                ["name", "value", "pressure", "volume", "gamma"]
            ]
        )

        return

    df = pd.read_csv(data_path)

    # target_gamma = 0.2
    df_ED = df[np.isclose(df["time"], 1)]
    df_ED = df_ED.assign(label="ED")

    # traget_pressure = 15.0
    df_ES = df[np.isclose(df["time"], 2)]
    df_ES = df_ES.assign(label="ES")
    df1 = pd.concat([df_ED, df_ES])

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


def postprocess_lv_ES(nativedir, transplanteddir, figdir):
    try:
        df_native = pd.read_csv(nativedir / "results.csv")
        df_trans = pd.read_csv(transplanteddir / "results.csv")
    except FileNotFoundError:
        print("No results found. Please run postprocess_lv first.")
        return

    # breakpoint()

    df_native_ES = df_native[np.isclose(df_native["time"], 2)]
    df_native_ES = df_native_ES.assign(label="Native")

    df_trans_ES = df_trans[np.isclose(df_trans["time"], 2)]
    df_trans_ES = df_trans_ES.assign(label="Transplanted")

    df1 = pd.concat([df_trans_ES, df_native_ES])

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
    fig.savefig(figdir / "stress_ES.svg")  # type: ignore
    plt.close(fig)
