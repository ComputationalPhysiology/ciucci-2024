from pathlib import Path
import dolfin
import pulse
import numpy as np
import ufl_legacy as ufl


from geometry import get_cylinder_geometry
from utils import Projector, ca_transient


dolfin.parameters["form_compiler"]["quadrature_degree"] = 6
dolfin.parameters["form_compiler"]["cpp_optimize"] = True
dolfin.parameters["form_compiler"]["representation"] = "uflacs"
dolfin.parameters["form_compiler"]["optimize"] = True

m2um = 1e6


def run(
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

    geo = get_cylinder_geometry(mesh_folder=mesh_folder)._asdict()

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
    u, p = problem.state.split()
    save(u, p, 0.0)
    for i, (ti, g) in enumerate(zip(t, gamma_arr)):
        pulse.iterate.iterate(problem, gamma, g, initial_number_of_steps=20)
        u, p = problem.state.split()
        save(u, p, ti)


def main(
    output_folder,
    mesh_folder="meshes-cylinder",
    varying_gamma=True,
    amp=0.2,
):
    springs = [3000, 100, 10000]

    output_folders = {spring: f"{output_folder.as_posix()}/spring{spring}" for spring in springs}
    for spring in springs:
        print("Spring : ", spring)
        if Path(output_folders[spring]).is_dir():
            continue
        run(
            output_folder=output_folders[spring],
            spring=spring,
            mesh_folder=mesh_folder,
            varying_gamma=varying_gamma,
            amp=amp,
        )

    # def key2title(k):
    #     if np.isclose(k, 100.0):
    #         return "Unloaded\n$k$ = 0.1 Pa/\u03BCm"
    #     elif np.isclose(k, 3000.0):
    #         return "Standard load\n$k$ = 1 Pa/\u03BCm"
    #     elif np.isclose(k, 10000.0):
    #         return "Increased load\n$k$ = 10 Pa/\u03BCm"

    # subfigdir = (
    #     "cylinder_fine_varying_incomp" if varying_gamma else "cylinder_fine_incomp"
    # )
    # figdir = f"{main_figdir}/{subfigdir}"

    # postprocess.postprocess_stress(
    #     resultdirs=resultdirs,
    #     figdir=figdir,
    #     key2title=key2title,
    #     datadir=datadir,
    #     plot_slice=True,
    # )
    # postprocess.postprocess_disp(
    #     resultdirs=resultdirs,
    #     figdir=figdir,
    #     key2title=key2title,
    #     datadir=datadir,
    # )


def main_twitch(
    output_folder: Path,
    mesh_folder="meshes/cylinder_fine2",
    amp=0.2,
    spring: float = 3000.0,
    varying_gamma=False,
):
    if not output_folder.is_dir():
        run(
            output_folder=output_folder,
            mesh_folder=mesh_folder,
            varying_gamma=varying_gamma,
            twitch=True,
            spring=spring,
            amp=amp,
        )

        # subfigdir = (
        #     "cylinder_fine_varying_incomp_twitch"
        #     if varying_gamma
        #     else "cylinder_fine_incomp_twitch"
        # )
        # figdir = f"{main_figdir}/{subfigdir}"

        # postprocess.postprocess_twitch_disp(
        #     resultdirs={0: resultsdir},
        #     figdir=figdir,
        #     key2title=lambda x: str(x),
        #     datadir=datadir,
        # )
        # postprocess.postprocess_twitch_stress(
        #     resultdirs={0: resultsdir},
        #     figdir=figdir,
        #     key2title=lambda x: str(x),
        #     datadir=datadir,
        # )
