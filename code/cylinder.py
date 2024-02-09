from pathlib import Path
import dolfin
import pulse
import numpy as np
import ufl_legacy as ufl


from geometry import load_cylinder_geometry
from utils import Projector

# import postprocess

dolfin.parameters["form_compiler"]["quadrature_degree"] = 6
dolfin.parameters["form_compiler"]["cpp_optimize"] = True
dolfin.parameters["form_compiler"]["representation"] = "uflacs"
dolfin.parameters["form_compiler"]["optimize"] = True

m2um = 1e6


def ca_transient(t, tstart=0.05, ca_ampl=0.3):
    tau1 = 0.05
    tau2 = 0.110

    ca_diast = 0.0

    beta = (tau1 / tau2) ** (-1 / (tau1 / tau2 - 1)) - (tau1 / tau2) ** (-1 / (1 - tau2 / tau1))
    ca = np.zeros_like(t)

    ca[t <= tstart] = ca_diast

    ca[t > tstart] = (ca_ampl - ca_diast) / beta * (
        np.exp(-(t[t > tstart] - tstart) / tau1) - np.exp(-(t[t > tstart] - tstart) / tau2)
    ) + ca_diast
    return ca


def main(
    output_folder="results",
    mesh_folder=Path("meshes/cylinder"),
    spring=3000.0,
    target_preload=0.0,
    overwrite: bool = False,
    varying_gamma: bool = False,
    twitch: bool = False,
    amp=0.2,
):
    pulse.set_log_level(10)
    Path(output_folder).mkdir(exist_ok=True, parents=True)
    output = Path(output_folder) / "results.xdmf"
    if output.is_file() and not overwrite:
        print(f"Output {output} already exists")
        return
    output.unlink(missing_ok=True)
    output.with_suffix(".h5").unlink(missing_ok=True)

    geo = load_cylinder_geometry(mesh_folder=mesh_folder)._asdict()

    microstructure = pulse.Microstructure(
        f0=dolfin.as_vector([1.0, 0.0, 0.0]),
        s0=dolfin.as_vector([0.0, 1.0, 0.0]),
        n0=dolfin.as_vector([0.0, 0.0, 1.0]),
    )

    geometry = pulse.Geometry(
        mesh=geo["mesh"],
        markers=geo["markers"],
        marker_functions=geo["marker_functions"],
        microstructure=microstructure,
    )

    gamma = dolfin.Constant(0.0)
    if varying_gamma:
        V_CG1 = dolfin.FunctionSpace(geometry.mesh, "CG", 1)
        R = geometry.mesh.coordinates().max(0)[-1]
        # gamma_expr = dolfin.Expression("(x[1]*x[1] + x[2]*x[2])/(R * R)", R=R, degree=2)
        gamma_expr = dolfin.Expression("(sqrt(x[1]*x[1] + x[2]*x[2]))/R", R=R, degree=2)

        activation = dolfin.Function(V_CG1)
        activation.interpolate(gamma_expr)
        activation *= gamma
    else:
        activation = gamma

    # m2mm = 1000.0

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

    material = pulse.HolzapfelOgden(
        active_model="active_strain",
        activation=activation,
        parameters=matparams,
        f0=geometry.f0,
        s0=geometry.s0,
        n0=geometry.n0,
    )
    robin_bc = [
        pulse.RobinBC(value=dolfin.Constant(spring / m2um), marker=geometry.markers["Right"][0]),
        pulse.RobinBC(value=dolfin.Constant(spring / m2um), marker=geometry.markers["Left"][0]),
    ]

    preload = dolfin.Constant(0.0)

    neumann_bc = [
        pulse.NeumannBC(traction=preload, marker=geometry.markers["Right"][0]),
        pulse.NeumannBC(traction=-preload, marker=geometry.markers["Left"][0]),
    ]

    bcs = pulse.BoundaryConditions(dirichlet=[lambda x: []], neumann=neumann_bc, robin=robin_bc)

    problem = pulse.MechanicsProblem(geometry, material, bcs)

    V_DG1 = dolfin.FunctionSpace(geometry.mesh, "DG", 1)
    proj = Projector(V_DG1)
    sigma_xx = dolfin.Function(V_DG1)
    sigma_r = dolfin.Function(V_DG1)
    sigma_c = dolfin.Function(V_DG1)

    sigma_dev_xx = dolfin.Function(V_DG1)
    sigma_dev_r = dolfin.Function(V_DG1)
    sigma_dev_c = dolfin.Function(V_DG1)

    E_xx = dolfin.Function(V_DG1)
    E_r = dolfin.Function(V_DG1)
    E_c = dolfin.Function(V_DG1)

    def save(U, p, ti):
        F = ufl.variable(pulse.kinematics.DeformationGradient(U))
        J = ufl.variable(ufl.det(F))

        sigma = material.CauchyStress(F, p)
        sigma_iso = sigma - (1 / 3) * ufl.tr(sigma) * ufl.Identity(3)
        E = pulse.kinematics.GreenLagrangeStrain(F)
        rad = F * geo["rad0"]
        circ = F * geo["circ0"]
        long = F * geo["long0"]

        proj.project(sigma_xx, dolfin.inner(long, sigma * long))
        proj.project(sigma_r, dolfin.inner(rad, sigma * rad))
        proj.project(sigma_c, dolfin.inner(circ, sigma * circ))

        proj.project(sigma_dev_xx, dolfin.inner(long, sigma_iso * long))
        proj.project(sigma_dev_r, dolfin.inner(rad, sigma_iso * rad))
        proj.project(sigma_dev_c, dolfin.inner(circ, sigma_iso * circ))

        proj.project(E_xx, dolfin.inner(geo["long0"], E * geo["long0"]))
        proj.project(E_r, dolfin.inner(geo["rad0"], E * geo["rad0"]))
        proj.project(E_c, dolfin.inner(geo["circ0"], E * geo["circ0"]))
        with dolfin.XDMFFile(output.as_posix()) as f:
            f.parameters["functions_share_mesh"] = True
            f.parameters["rewrite_function_mesh"] = False
            f.write_checkpoint(
                U,
                function_name="u",
                time_step=ti,
                encoding=dolfin.XDMFFile.Encoding.HDF5,
                append=True,
            )
            f.write_checkpoint(
                sigma_xx,
                function_name="sigma_xx",
                time_step=ti,
                encoding=dolfin.XDMFFile.Encoding.HDF5,
                append=True,
            )
            f.write_checkpoint(
                sigma_r,
                function_name="sigma_r",
                time_step=ti,
                encoding=dolfin.XDMFFile.Encoding.HDF5,
                append=True,
            )
            f.write_checkpoint(
                sigma_c,
                function_name="sigma_c",
                time_step=ti,
                encoding=dolfin.XDMFFile.Encoding.HDF5,
                append=True,
            )
            f.write_checkpoint(
                sigma_dev_xx,
                function_name="sigma_dev_xx",
                time_step=ti,
                encoding=dolfin.XDMFFile.Encoding.HDF5,
                append=True,
            )
            f.write_checkpoint(
                sigma_dev_r,
                function_name="sigma_dev_r",
                time_step=ti,
                encoding=dolfin.XDMFFile.Encoding.HDF5,
                append=True,
            )
            f.write_checkpoint(
                sigma_dev_c,
                function_name="sigma_dev_c",
                time_step=ti,
                encoding=dolfin.XDMFFile.Encoding.HDF5,
                append=True,
            )
            f.write_checkpoint(
                E_xx,
                function_name="E_xx",
                time_step=ti,
                encoding=dolfin.XDMFFile.Encoding.HDF5,
                append=True,
            )
            f.write_checkpoint(
                E_r,
                function_name="E_r",
                time_step=ti,
                encoding=dolfin.XDMFFile.Encoding.HDF5,
                append=True,
            )
            f.write_checkpoint(
                E_c,
                function_name="E_c",
                time_step=ti,
                encoding=dolfin.XDMFFile.Encoding.HDF5,
                append=True,
            )
            f.write_checkpoint(
                p,
                function_name="hydrostatic",
                time_step=ti,
                encoding=dolfin.XDMFFile.Encoding.HDF5,
                append=True,
            )

    if twitch:
        N = 50
        t = np.linspace(0, 1, N)
        gamma_arr = ca_transient(t, ca_ampl=amp)
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots()
        ax.plot(t, gamma_arr)
        ax.set_xlabel("Time [s]")
        ax.set_ylabel("Activation [-]")
        fig.savefig(output_folder / "activation.png")
    else:
        t = [0.2]
        gamma_arr = [amp]

    np.save(Path(output_folder) / "t.npy", [0.0] + t)
    np.save(Path(output_folder) / "gamma.npy", [0.0] + gamma_arr)

    if not np.isclose(target_preload, 0.0):
        pulse.iterate.iterate(problem, preload, target_preload, initial_number_of_steps=20)

    save(*problem.state.split(), 0.0)
    for i, (ti, g) in enumerate(zip(t, gamma_arr)):
        pulse.iterate.iterate(problem, gamma, g, initial_number_of_steps=20)
        save(*problem.state.split(), ti)


def effect_of_spring(
    main_resultdir,
    main_figdir="figures2",
    datadir="meshes/cylinder_fine2",
    varying_gamma=True,
    amp=0.2,
):
    springs = [3000.0, 100.0, 10000.0]

    subfolder = "cylinder_fine_varying_incomp" if varying_gamma else "cylinder_fine_incomp"
    resultdirs = {
        spring: f"{main_resultdir.as_posix()}/{subfolder}/spring{spring}" for spring in springs
    }
    for spring in springs:
        print("Spring : ", spring)
        if Path(resultdirs[spring]).is_dir():
            continue
        run(
            resultsdir=resultdirs[spring],
            spring=spring,
            datadir=datadir,
            varying_gamma=varying_gamma,
            amp=amp,
        )

    def key2title(k):
        if np.isclose(k, 100.0):
            return "Unloaded\n$k$ = 0.1 Pa/\u03BCm"
        elif np.isclose(k, 3000.0):
            return "Standard load\n$k$ = 1 Pa/\u03BCm"
        elif np.isclose(k, 10000.0):
            return "Increased load\n$k$ = 10 Pa/\u03BCm"

    subfigdir = "cylinder_fine_varying_incomp" if varying_gamma else "cylinder_fine_incomp"
    figdir = f"{main_figdir}/{subfigdir}"

    postprocess.postprocess_stress(
        resultdirs=resultdirs,
        figdir=figdir,
        key2title=key2title,
        datadir=datadir,
        plot_slice=True,
    )
    postprocess.postprocess_disp(
        resultdirs=resultdirs,
        figdir=figdir,
        key2title=key2title,
        datadir=datadir,
    )


def run_compressiblity(main_resultdir):
    kappas = [1, 10, 1e2, 1e3, 1e4]
    resultdirs = {kappa: f"{main_resultdir.as_posix()}/kappa{kappa}" for kappa in kappas}
    for kappa in kappas:
        run(
            resultsdir=resultdirs[kappa],
            overwrite=True,
            datadir="meshes/cylin",
            kappa=kappa,
        )

    key2title = lambda kappa: "incomp" if kappa is None else rf"$\kappa = {kappa:.0f}$"
    postprocess.postprocess_stress(
        resultdirs=resultdirs, figdir="figures/cylinder/incomp", key2title=key2title
    )
    postprocess.postprocess_disp(
        resultdirs=resultdirs, figdir="figures/cylinder/incomp", key2title=key2title
    )


def run_basic_with_preload(main_resultdir):
    preloads = [0.0, 0.1, 0.5, 1.0, 1.5]
    resultdirs = {preload: f"{main_resultdir.as_posix()}/preload_{preload}" for preload in preloads}
    for preload in preloads:
        run(
            resultsdir=resultdirs[preload],
            datadir="meshes/cylin",
            overwrite=True,
            target_preload=preload,
        )

    key2title = lambda preload: rf"$F = {preload:.3f}$"
    postprocess.postprocess_stress(
        resultdirs=resultdirs, figdir="figures/cylinder/preload", key2title=key2title
    )
    postprocess.postprocess_disp(
        resultdirs=resultdirs, figdir="figures/cylinder/preload", key2title=key2title
    )


def run_twitch(main_resultdir, main_figdir="figures2", datadir="meshes/cylinder_fine2", amp=0.2):
    for varying_gamma in [False, True]:
        resultsdir = (
            main_resultdir / "cylinder_fine_varying_incomp" / "twitch"
            if varying_gamma
            else main_resultdir / "cylinder_fine_incomp" / "twitch"
        )

        if not resultsdir.is_dir():
            run(
                resultsdir=resultsdir,
                datadir=datadir,
                varying_gamma=varying_gamma,
                twitch=True,
                amp=amp,
            )

        subfigdir = (
            "cylinder_fine_varying_incomp_twitch"
            if varying_gamma
            else "cylinder_fine_incomp_twitch"
        )
        figdir = f"{main_figdir}/{subfigdir}"

        postprocess.postprocess_twitch_disp(
            resultdirs={0: resultsdir},
            figdir=figdir,
            key2title=lambda x: str(x),
            datadir=datadir,
        )
        postprocess.postprocess_twitch_stress(
            resultdirs={0: resultsdir},
            figdir=figdir,
            key2title=lambda x: str(x),
            datadir=datadir,
        )


# def main():
#     # main_resultdir = Path("results/cylinder")
#     # main_resultdir.mkdir(exist_ok=True, parents=True)
#     # # run_basic_with_preload(main_resultdir)
#     # effect_of_spring(main_resultdir)
#     # # run_compressiblity(main_resultdir)

#     # main_resultdir = Path("results/cylinder_varying_incomp")
#     # main_resultdir.mkdir(exist_ok=True, parents=True)
#     # run_basic_with_preload(main_resultdir)

#     main_resultdir = Path("results4")
#     main_resultdir.mkdir(exist_ok=True, parents=True)
#     effect_of_spring(main_resultdir, main_figdir="figures4", varying_gamma=False)
#     effect_of_spring(main_resultdir, main_figdir="figures4", varying_gamma=True)
#     run_twitch(main_resultdir, main_figdir="figures4")

#     main_resultdir = Path("results5")
#     main_resultdir.mkdir(exist_ok=True, parents=True)
#     effect_of_spring(
#         main_resultdir, main_figdir="figures5", varying_gamma=False, amp=0.25
#     )
#     effect_of_spring(
#         main_resultdir, main_figdir="figures5", varying_gamma=True, amp=0.25
#     )
#     run_twitch(main_resultdir, main_figdir="figures5", amp=0.25)

#     main_resultdir = Path("results6")
#     main_resultdir.mkdir(exist_ok=True, parents=True)
#     effect_of_spring(
#         main_resultdir, main_figdir="figures6", varying_gamma=False, amp=0.3
#     )
#     effect_of_spring(
#         main_resultdir, main_figdir="figures6", varying_gamma=True, amp=0.3
#     )
#     run_twitch(main_resultdir, main_figdir="figures6", amp=0.3)


# if __name__ == "__main__":
#     main()
