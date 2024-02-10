from pathlib import Path
import dolfin
import pandas as pd
import seaborn as sns

import ufl_legacy as ufl
import numpy as np
import matplotlib.pyplot as plt

from geometry import get_cylinder_geometry
from utils import name2latex


def create_arrays(
    resultsdirs,
    data_file,
    data_file_points,
    plot_file,
    gamma,
    geo,
    key2title,
    is_twitch=False,
):
    R = geo.mesh.coordinates().max(0)[-1]
    L = geo.mesh.coordinates().max(0)[0]

    N = len(gamma)

    V_DG1 = dolfin.FunctionSpace(geo.mesh, "DG", 1)

    dr = R / 5.0
    points = np.arange(0, R + dr, dr)

    funcs = {
        name: dolfin.Function(V_DG1)
        for name in [
            "sigma_xx",
            "sigma_r",
            "sigma_c",
            "sigma_dev_xx",
            "sigma_dev_r",
            "sigma_dev_c",
            "E_xx",
            "E_r",
            "E_c",
            "hydrostatic",
        ]
    }

    for v in funcs.values():
        v.set_allow_extrapolation(True)

    data = []
    data_points = []
    plots = {}

    for idx, (key, resultsdir) in enumerate(resultsdirs.items()):
        output = Path(resultsdir) / "results.xdmf"
        print(output)
        plots[key] = {}

        with dolfin.XDMFFile(output.as_posix()) as f:
            for i in range(N):
                print("Timepoint", i + 1, "of", N, "for", key, "with gamma", gamma[i])

                for name, func in funcs.items():
                    print(name)

                    f.read_checkpoint(func, name, i)

                    if name == "sigma_xx":
                        name = "sigma_ff"
                    if name == "sigma_dev_xx":
                        name = "sigma_dev_ff"
                    if name == "hydrostatic":
                        name = "p"

                    if not is_twitch:
                        # Average over the dofs in the center
                        dofs = np.abs(V_DG1.tabulate_dof_coordinates()[:, 0]) < 0.1 * L
                        f_arr = func.vector().get_local()[dofs]

                        data.extend(
                            [
                                {
                                    "name": name,
                                    "key": key,
                                    "value": fi,
                                    "gamma": gamma[i],
                                    "latex": name2latex(name),
                                    "xlabel": key2title(key),
                                }
                                for fi in f_arr
                            ]
                        )

                    if not is_twitch and i == 0:
                        continue

                    for point_start, point_end in zip(points[:-1], points[1:]):
                        for point in np.arange(point_start, point_end, dr / 5):
                            for y_ in [-1, 0, 1, np.sqrt(2) / 2, -np.sqrt(2) / 2]:
                                y = y_ * point
                                for z in [
                                    np.sqrt(point**2 - y**2),
                                    -np.sqrt(point**2 - y**2),
                                ]:
                                    print(point, y, z, np.sqrt(y**2 + z**2))
                                    data_points.append(
                                        {
                                            "name": name,
                                            "key": key,
                                            "value": func(0, y, z),
                                            "gamma": gamma[i] * point**2 / R**2,
                                            "latex": name2latex(name),
                                            "xlabel": key2title(key),
                                            "point": point,
                                            "timepoint": i,
                                            "point_start": point_start,
                                            "point_end": point_end,
                                        }
                                    )

    if data_file is not None:
        df = pd.DataFrame(data)
        df.to_csv(data_file)

    df_points = pd.DataFrame(data_points)
    df_points.to_csv(data_file_points)
    if plot_file is not None:
        np.save(plot_file, plots, allow_pickle=True)


def plot_points_in_slice(data_file, figdir):
    df = pd.read_csv(data_file)

    for i, name in enumerate(
        ["sigma_ff", "sigma_r", "sigma_c", "sigma_dev_ff", "sigma_dev_r", "sigma_dev_c"]
    ):
        df_plot = df[df["name"] == name]
        fig = plt.figure()
        ax = sns.barplot(
            data=df_plot,
            x="xlabel",
            y="value",
            hue="point",
            alpha=0.7,
            errorbar="ci",
        )
        ax.get_legend().set_title("Radius")
        ax.set_xlabel("")
        ax.grid()
        ax.set_ylabel(f"{name2latex(name)} [kPa]")
        fig.tight_layout()
        fig.savefig(figdir / f"{name}_points.png")
        fig.savefig(figdir / f"{name}_points.svg")

    for i, name in enumerate(["E_xx", "E_r", "E_c"]):
        df_plot = df[df["name"] == name]
        fig = plt.figure()
        ax = sns.barplot(
            data=df_plot,
            x="xlabel",
            y="value",
            hue="point",
            alpha=0.7,
            errorbar="ci",
        )
        ax.get_legend().set_title("Radius")
        ax.set_xlabel("")
        ax.grid()
        ax.set_ylabel(f"{name2latex(name)}")
        fig.tight_layout()
        fig.savefig(figdir / f"{name}_points.png")

    # df_plot = df[df["name"] == "hydrostatic2"]
    df_plot = df[df["name"] == "p"]
    fig = plt.figure()
    ax = sns.barplot(
        data=df_plot,
        x="xlabel",
        y="value",
        hue="point",
        alpha=0.7,
        errorbar="ci",
    )
    ax.get_legend().set_title("Radius")
    ax.set_xlabel("")
    ax.grid()
    ax.set_ylabel("$p$ [kPa]")
    fig.tight_layout()
    fig.savefig(figdir / "p_points.png")
    fig.savefig(figdir / "p_points.svg")


def plot_stress(data_file, figdir):
    df = pd.read_csv(data_file)
    max_df = df[np.isclose(df["gamma"], df["gamma"].max())]
    stress_df = max_df[max_df["name"].isin(["sigma_ff", "sigma_r", "sigma_c"])]
    # sns.set(font_scale=2)
    plt.rcParams.update({"font.size": 16})
    fig = plt.figure()
    ax = sns.barplot(
        data=stress_df,
        y="value",
        x="xlabel",
        hue="latex",
        errorbar="ci",
        alpha=0.7,
    )
    ax.get_legend().set_title("")
    ax.set_xlabel("")
    # ax.grid()
    ax.set_ylabel("Average stress [kPa]")
    fig.tight_layout()
    fig.savefig(figdir / "stress.svg")
    fig.savefig(figdir / "stress.png")  # , dpi=500)

    stress_dev_df = max_df[max_df["name"].isin(["sigma_dev_ff", "sigma_dev_r", "sigma_dev_c", "p"])]
    # sns.set(font_scale=2)
    plt.rcParams.update({"font.size": 16})
    fig = plt.figure()
    ax = sns.barplot(
        data=stress_dev_df,
        y="value",
        x="xlabel",
        hue="latex",
        errorbar="ci",
        alpha=0.7,
    )
    ax.get_legend().set_title("")
    ax.set_xlabel("")
    # ax.grid()
    ax.set_ylabel("Average stress [kPa]")
    fig.tight_layout()
    fig.savefig(figdir / "stress_dev.svg")
    fig.savefig(figdir / "stress_dev.png")  # , dpi=500)

    fig = plt.figure()
    ax = sns.barplot(
        data=stress_df,
        x="xlabel",
        y="value",
        hue="latex",
        errorbar="sd",
        alpha=0.7,
    )
    ax.get_legend().set_title("")
    ax.set_xlabel("")
    ax.set_ylabel("Average stress [kPa]")
    fig.savefig(figdir / "stress_yerr.png")

    ax.set_yscale("log")
    fig.tight_layout()
    fig.savefig(figdir / "stress_log.png")
    plt.close(fig)

    fig = plt.figure()
    strain_df = max_df[max_df["name"].isin(["E_xx", "E_r", "E_c"])]
    ax = sns.barplot(
        data=strain_df,
        x="xlabel",
        y="value",
        hue="latex",
        alpha=0.7,
    )
    ax.get_legend().set_title("")
    ax.set_xlabel("")
    ax.set_ylabel("Average strain")
    fig.tight_layout()
    fig.savefig(figdir / "strain.png")
    fig.savefig(figdir / "strain.svg")
    plt.close(fig)


# def postprocess_stress(
#     resultsdirs,
#     key2title,
#     mesh_folder="meshes-cylinder",
#     figdir="figures-cylinder",
#     plot_slice=False,
# ):
#     geo = get_cylinder_geometry(mesh_folder=mesh_folder)
#     figdir = Path(figdir)
#     figdir.mkdir(exist_ok=True, parents=True)

#     fst = next(iter(resultsdirs.values()))
#     t = np.load(Path(fst) / "t.npy")
#     gamma = np.load(Path(fst) / "gamma.npy")

#     R = geo.mesh.coordinates().max(0)[-1]
#     disk = dolfin.UnitDiscMesh.create(dolfin.MPI.comm_world, 50, 1, 2)
#     disk.coordinates()[:] *= R
#     V_disk = dolfin.FunctionSpace(disk, "CG", 1)

#     data_file = figdir / "data.csv"
#     data_file_points = figdir / "data_points.csv"
#     plot_file = figdir / "plots.npy"
#     if not data_file.exists():
#         create_arrays(
#             resultsdirs=resultsdirs,
#             gamma=gamma,
#             figdir=figdir,
#             plot_slice=plot_slice,
#             geo=geo,
#             key2title=key2title,
#             data_file=data_file,
#             data_file_points=data_file_points,
#             plot_file=plot_file,
#             V_disk=V_disk,
#         )

#     # data_arrs = np.load(arr_path, allow_pickle=True).item()
#     # arrs_std = data_arrs["arrs_std"]
#     # arrs_avg = data_arrs["arrs_avg"]
#     plot_stress(data_file, figdir)
#     # plot_slices(V_disk, plot_file, figdir)
#     plot_relevant_slices(V_disk, plot_file, figdir)
#     plot_points_in_slice(data_file_points, figdir)


def load_displacement(geo, resultsdirs, data_path, key2title):
    R = geo.mesh.coordinates().max(0)[-1]
    L = geo.mesh.coordinates().max(0)[0]

    V_CG2 = dolfin.VectorFunctionSpace(geo.mesh, "CG", 2)
    u = dolfin.Function(V_CG2)
    u.set_allow_extrapolation(True)

    num_dirs = len(resultsdirs)
    fst = next(iter(resultsdirs.values()))
    t = np.load(Path(fst) / "t.npy")
    N = len(t)

    vols = np.zeros((num_dirs, N))
    mesh_vol = dolfin.assemble(dolfin.Constant(1) * dolfin.dx(geo.mesh))
    data = []
    for idx, (key, resultsdir) in enumerate(resultsdirs.items()):
        output = Path(resultsdir) / "results.xdmf"
        print(output)

        with dolfin.XDMFFile(output.as_posix()) as f:
            for i in range(N):
                f.read_checkpoint(u, "u", i)

                vols[idx, i] = (
                    dolfin.assemble(ufl.det(ufl.grad(u) + ufl.Identity(3)) * dolfin.dx) / mesh_vol
                )
                u_center = u(0, 0, R)
                u_long = u(L, 0, 0)
                data.append(
                    {
                        "u_center_x": u_center[0],
                        "u_center_y": u_center[1],
                        "u_center_z": u_center[2],
                        "u_long_x": u_long[0],
                        "u_long_y": u_long[1],
                        "u_long_z": u_long[2],
                        "length": 2 * (L + u_long[0]),
                        "fractional_shortening": abs(u_long[0] / L),
                        "diameter": 2 * (R + np.sqrt(u_center[1] ** 2 + u_center[2] ** 2)),
                        "key": key,
                        "volume": vols[idx, i],
                        "gamma": t[i],
                        "title": key2title(key),
                        "time_point": i,
                        "label": "Contaction" if i == 1 else "Relaxation",
                    }
                )
    df = pd.DataFrame(data)
    df.to_csv(data_path)


def postprocess_disp(
    data_path,
    figdir="figures-cylinder",
):
    df = pd.read_csv(data_path)
    df_rc = df[df["time_point"].isin([0, 1])]
    df_c = df_rc[df_rc["time_point"] == 1]
    k3 = df_c[df_c["title"] == "Standard load\n$k$ = 3 Pa/μm"]
    k3 = k3.assign(k=3.0)
    k01 = df_c[df_c["title"] == "Unloaded\n$k$ = 0.1 Pa/μm"]
    k01 = k01.assign(k=0.1)
    k10 = df_c[df_c["title"] == "Increased load\n$k$ = 10 Pa/μm"]
    k10 = k10.assign(k=10.0)
    df_c = pd.concat([k01, k3, k10])

    plt.rcParams.update({"font.size": 16})
    fig = plt.figure(figsize=(4, 4))
    ax = sns.lineplot(
        data=df_c,
        x="k",
        y="fractional_shortening",
        hue="label",
        alpha=0.7,
        marker="o",
        linewidth=3,
        markersize=10,
        legend=False,
    )
    ax.set_xscale("log")
    ax.set_xticks([0.1, 1, 10])
    ax.set_xticklabels(["0.1", "1", "10"])
    ax.set_xlabel("$k$ [Pa/μm]")
    # sns.move_legend(
    #     ax, "lower center", bbox_to_anchor=(0.5, 1), ncol=3, title=None, frameon=False
    # )

    ax.set_ylim(0, 0.25)
    ax.set_ylabel("Fractional shortening")
    ax.grid()

    ax.set_yticks(np.arange(0, 0.3, 0.05))
    ax.set_yticklabels(["0%", "5%", "10%", "15%", "20%", "25%"])

    # ax.yaxis.set_major_formatter(ticker.PercentFormatter())
    fig.tight_layout()
    fig.savefig(figdir / "frac.svg")
    fig.savefig(figdir / "frac.png")
    plt.close(fig)

    fig = plt.figure()
    ax = sns.barplot(
        data=df_rc,
        x="title",
        y="length",
        hue="label",
        alpha=0.7,
    )
    sns.move_legend(ax, "lower center", bbox_to_anchor=(0.5, 1), ncol=3, title=None, frameon=False)
    ax.set_xlabel("")
    ax.set_ylabel("Length [mm]")
    ax.grid()
    fig.tight_layout()
    fig.savefig(figdir / "length.svg")
    fig.savefig(figdir / "length.png")
    plt.close(fig)

    fig = plt.figure()
    ax = sns.barplot(
        data=df_rc,
        x="label",
        y="volume",
        hue="title",
        alpha=0.7,
    )
    sns.move_legend(ax, "lower center", bbox_to_anchor=(0.5, 1), ncol=3, title=None, frameon=False)
    ax.set_xlabel("")
    ax.set_ylabel("Volume change")
    fig.savefig(figdir / "volume.png")
    plt.close(fig)

    fig = plt.figure()
    ax = sns.barplot(
        data=df_rc,
        x="title",
        hue="label",
        y="diameter",
        alpha=0.7,
    )
    ax.grid()
    sns.move_legend(ax, "lower center", bbox_to_anchor=(0.5, 1), ncol=3, title=None, frameon=False)
    ax.set_xlabel("")
    ax.set_ylabel("Diameter [mm]")
    fig.tight_layout()
    fig.savefig(figdir / "diameter.svg")
    fig.savefig(figdir / "diameter.png")
    plt.close(fig)


def _plot_twitch_single(
    df,
    y,
    ylabel,
    path,
    xticks,
    xticklabels,
    x="time_point",
    xlabel="Time",
    hue=None,
    lengend=False,
    marker="o",
    palette=None,
    errorbar=None,
):
    fig = plt.figure(figsize=(8, 8))
    ax = sns.lineplot(
        data=df,
        x=x,
        y=y,
        hue=hue,
        marker=marker,
        linewidth=3,
        markersize=10,
        legend=lengend,
        errorbar=errorbar,
        palette=palette,
    )

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid()
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticklabels)

    fig.tight_layout()
    fig.savefig(path.with_suffix(".svg"))
    fig.savefig(path.with_suffix(".png"))
    plt.close(fig)


def postprocess_twitch_disp(
    data_path,
    t: np.ndarray,
    figdir="figures-cylinder",
):
    figdir = Path(figdir)
    figdir.mkdir(exist_ok=True, parents=True)

    xticks = np.array([0, len(t) // 4, len(t) // 2, 3 * len(t) // 4, len(t) - 1])
    xticklabels = [f"{ti:.2f}" for ti in t[xticks]]

    df = pd.read_csv(data_path)

    plt.rcParams.update({"font.size": 16})

    # Fractional shortening
    _plot_twitch_single(
        df=df,
        y="fractional_shortening",
        ylabel="Fractional shortening",
        path=figdir / "frac.svg",
        xticks=xticks,
        xticklabels=xticklabels,
    )

    # Length
    _plot_twitch_single(
        df=df,
        y="length",
        ylabel="Length [mm]",
        path=figdir / "length.svg",
        xticks=xticks,
        xticklabels=xticklabels,
    )

    # Diameter
    _plot_twitch_single(
        df=df,
        y="diameter",
        ylabel="Diameter [mm]",
        path=figdir / "diameter.svg",
        xticks=xticks,
        xticklabels=xticklabels,
    )


def postprocess_twitch_stress(
    data_file,
    figdir="figures-cylinder",
):
    figdir = Path(figdir)
    figdir.mkdir(exist_ok=True, parents=True)
    df = pd.read_csv(data_file)

    points = df.point_end.unique()
    points_start = df.point_start.unique()
    points_str = [f"{p_start:.0f}-{p:.0f}" for p_start, p in zip(points_start, points)]
    # Get one color for each point using tab10
    colors = plt.get_cmap("viridis")(np.linspace(0, 1, len(points)))

    fig, axs = plt.subplots(2, 4, figsize=(16, 8), sharex=True)
    for i, name in enumerate(
        [
            "sigma_ff",
            "sigma_r",
            "sigma_c",
            "p",
            "sigma_dev_ff",
            "sigma_dev_r",
            "sigma_dev_c",
        ]
    ):
        df_plot = df[df["name"] == name]
        lines = []
        ax = axs.flatten()[i]
        for j, point in enumerate(points):
            df_point = df_plot[df_plot["point_end"] == point]
            value = df_point[["timepoint", "value"]].groupby("timepoint")
            value_mean = value.agg(lambda x: x.mean())
            value_std = value.agg(lambda x: x.std())

            (l,) = ax.plot(
                df_point["timepoint"].unique(),
                value_mean,
                label=point,
                color=colors[j],
            )
            lines.append(l)
            ax.fill_between(
                value_mean.index,
                value_mean["value"] - value_std["value"],
                value_mean["value"] + value_std["value"],
                alpha=0.2,
                color=colors[j],
            )

        ax.set_xlabel("Time")
        ax.set_ylabel("Average stress [kPa]")
        ax.set_title(name2latex(name))
        ax.grid()

    ax_legend = axs.flatten()[-1]
    ax_legend.axis("off")
    fig.legend(
        lines,
        points_str,
        loc="center",
        bbox_to_anchor=(0.9, 0.25),
        title="Radius",
        fontsize="large",
        title_fontsize="x-large",
    )
    fig.tight_layout()
    fig.savefig(figdir / "stress.svg")
    fig.savefig(figdir / "stress.png")
    plt.close(fig)

    fig = plt.figure()
    ax = fig.add_axes(
        [0, 0, 1, 1],
        polar=True,
    )

    for i, color in enumerate(colors):
        ax.bar(
            i * 2 * np.pi,
            60,
            width=2 * np.pi,
            bottom=60 * i,
            color=color,
            edgecolor=color,
        )

    ax.set_xticks([])
    ax.set_yticks(points)
    fig.savefig(figdir / "stress_polar_label.png")
    fig.savefig(figdir / "stress_polar_label.svg")


def postprocess_cylinder(
    resultsdir,
    mesh_folder="meshes-cylinder",
    figdir="figures-cylinder",
):
    springs = [3000, 100, 10000]

    resultsdirs = {spring: f"{resultsdir.as_posix()}/spring{spring}" for spring in springs}

    def key2title(k):
        if np.isclose(k, 100.0):
            return "Unloaded\n$k$ = 0.1 Pa/\u03BCm"
        elif np.isclose(k, 3000.0):
            return "Standard load\n$k$ = 3 Pa/\u03BCm"
        elif np.isclose(k, 10000.0):
            return "Increased load\n$k$ = 10 Pa/\u03BCm"

    geo = get_cylinder_geometry(mesh_folder=mesh_folder)
    figdir = Path(figdir)
    figdir.mkdir(exist_ok=True, parents=True)
    gamma = np.load(Path(next(iter(resultsdirs.values()))) / "gamma.npy")

    data_file = figdir / "data.csv"
    data_file_points = figdir / "data_points.csv"
    plot_file = figdir / "plots.npy"
    if not data_file.exists():
        create_arrays(
            resultsdirs=resultsdirs,
            gamma=gamma,
            geo=geo,
            key2title=key2title,
            data_file=data_file,
            data_file_points=data_file_points,
            plot_file=plot_file,
            is_twitch=False,
        )

    data_path_u = figdir / "data_u.csv"
    if not data_path_u.is_file():
        load_displacement(
            geo=geo,
            resultsdirs=resultsdirs,
            data_path=data_path_u,
            key2title=key2title,
        )

    postprocess_disp(data_path_u, figdir)
    plot_stress(data_file, figdir)
    plot_points_in_slice(data_file_points, figdir)


def postprocess_cylinder_twitch(
    resultsdir,
    mesh_folder="meshes-cylinder",
    figdir="figures-cylinder",
):
    geo = get_cylinder_geometry(mesh_folder=mesh_folder)
    figdir = Path(figdir)
    figdir.mkdir(exist_ok=True, parents=True)
    t = np.load(Path(resultsdir) / "t.npy")
    gamma = np.load(Path(resultsdir) / "gamma.npy")

    data_file = figdir / "data.csv"
    data_file_points = figdir / "data_points.csv"
    plot_file = figdir / "plots.npy"
    if not data_file.exists():
        create_arrays(
            resultsdirs={0: resultsdir},
            gamma=gamma,
            geo=geo,
            key2title=lambda x: str(x),
            data_file=data_file,
            data_file_points=data_file_points,
            plot_file=plot_file,
            is_twitch=True,
        )

    data_path_u = figdir / "data_u.csv"
    if not data_path_u.is_file():
        load_displacement(
            geo=geo,
            resultsdirs={0: resultsdir},
            data_path=data_path_u,
            key2title=lambda x: str(x),
        )

    postprocess_twitch_disp(data_path_u, t, figdir)
    postprocess_twitch_stress(data_file=data_file_points, figdir=figdir)
